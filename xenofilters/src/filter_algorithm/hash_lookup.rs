//! [`HashLookup`] -- two-pass fragment-matching for position-sorted BAMs.
//!
//! **Pass 1** (sequential scan, lightweight `MappedRecord`s):
//! Reads name, flags, ref_id, pos, CIGAR, MD, qualities, virtual_offset.
//! No sequence. Inserts into a `FragmentTable`. At fragment completion, classifies
//! each stream as Early (all primaries perfect, no BED/VCF overlap) or
//! NeedsScoring. Early fragments retain only virtual offsets; Scoring
//! fragments retain `MappedRecord`s for NW scoring.
//!
//! **Pass 2** (selective seek):
//! For each completed fragment, seeks to stored virtual offsets and reads
//! full records for output. Supplementary virtual offsets follow the
//! fragment's decision.
//!
//! Single-threaded only. Output order follows driving-stream (stream 0) order.

pub mod assemble;
pub(crate) mod name_codec;
pub(crate) mod stage;
#[cfg(test)]
pub(crate) mod tests;

use crate::alignment::pre_assess::min_delta_for_early_decision;
use crate::alignment::{mate_slot, segment_id, Fragment, MateKind, MdCigFlags, SimpleRec};
use crate::alignment::{pre_assess_scoring_records, PreAssessResult};
use crate::aln_stream::AlignmentStream;
use crate::config::args::resolve_threshold;
use crate::config::run_config::RunConfig;
use crate::filter_algorithm::line_by_line::COUNTER_STRIDE;
use crate::filter_algorithm::line_by_line::{READ_CT, Scratch, ordering::Decision};
use crate::penalty::Penalty;
use crate::region::tabix_query::{TabixBed, TabixVcf};
use crate::region::{ScoreFn, ScoredRegions};
use crate::variant::FragEvalVec;
use crate::Error;
use assemble::{
    EarlyKind, FragmentTable, MappedRecord, PendingFragment, RecordKind, StreamKind, insert,
    new_fragment_table,
};
use noodles::core::Position;
use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record_buf::{
    Cigar, Data, QualityScores, RecordBuf, Sequence, data::field::Value as BufValue,
};
use smallvec::SmallVec;
use stage::StagedOutput;
use std::cmp::Ordering;

// ---------------------------------------------------------------------------
// ScoredFragment
// ---------------------------------------------------------------------------

pub(crate) struct ScoredFragment {
    /// (stream_nr, virtual_offset) for winning records.
    pub(crate) winner_offsets: SmallVec<[(usize, u64); 2]>,
    /// (stream_nr, virtual_offset) for losing records.
    pub(crate) loser_offsets: SmallVec<[(usize, u64); 2]>,
    /// Supplementary offsets per stream -- follow winner's decision.
    pub(crate) supp_offsets: [SmallVec<[u64; 1]>; 2],
    pub(crate) decision: Option<Decision>,
    pub(crate) winner_nr: usize,
    /// True if result is ambiguous (both streams go to ambiguous output).
    pub(crate) is_ambiguous: bool,
}

// ---------------------------------------------------------------------------
// HashLookup
// ---------------------------------------------------------------------------

pub(crate) struct HashLookup<R: SimpleRec> {
    aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    table: FragmentTable,
    staged: StagedOutput,
    seq_counter: u64,
    /// Per-stream record counters used as virtual offsets for MockStream compat
    /// and correct pass-2 fetch indexing.
    record_counters: [u64; 2],
    penalties: Penalty,
    scratch: Scratch,
    pub(crate) routing_counters: SmallVec<[u64; 8]>,
    add_decision_tag: bool,
    ambiguous_log_threshold: f64,
    bed: [Option<TabixBed>; 2],
    vcf: [Option<TabixVcf>; 2],
    bisulfite: bool,
    positive: [Option<(ScoredRegions, ScoreFn)>; 2],
    codec_prefix: Vec<u8>,
}

impl<R: SimpleRec> HashLookup<R> {
    pub(crate) fn new(
        args: &RunConfig,
        aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>,
        bed: [Option<TabixBed>; 2],
        vcf: [Option<TabixVcf>; 2],
        pos: [Option<(ScoredRegions, ScoreFn)>; 2],
    ) -> Result<HashLookup<RecordBuf>, Error> {
        let ambiguous_log_threshold = resolve_threshold(args.scoring.ambiguous_threshold, false);
        Ok(HashLookup {
            aln,
            table: new_fragment_table(),
            staged: StagedOutput::new(),
            seq_counter: 0,
            record_counters: [0, 0],
            penalties: args.scoring.to_penalty(),
            scratch: Scratch::new(),
            routing_counters: SmallVec::from_elem(0, 2 * COUNTER_STRIDE),
            add_decision_tag: args.io.add_decision_tag,
            ambiguous_log_threshold,
            bed,
            vcf,
            positive: pos,
            codec_prefix: Vec::new(),
            bisulfite: args.scoring.bisulfite,
        })
    }

    pub(crate) fn process(&mut self, config: &RunConfig) -> Result<(), Error> {
        let mut exhausted = [false; 2];
        loop {
            let mut finished = true;
            for (nr, ex) in exhausted.iter_mut().enumerate() {
                if *ex {
                    continue;
                }
                match self.next_scoring_record(nr)? {
                    None => *ex = true,
                    Some((key, rec)) => {
                        finished = false;
                        self.ingest(key, rec, nr)?;
                    }
                }
            }
            if finished {
                break;
            }
            self.staged.flush(
                &mut self.aln,
                &mut self.routing_counters,
                self.add_decision_tag,
            )?;
        }

        let unmatched: Vec<_> = self.table.drain().collect();
        for (_, pending) in unmatched {
            self.emit_unmatched(pending)?;
        }
        self.staged.flush_all(
            &mut self.aln,
            &mut self.routing_counters,
            self.add_decision_tag,
        )?;
        config.print_routing_counters(&self.routing_counters, "Hash-lookup");
        Ok(())
    }

    fn next_scoring_record(&mut self, nr: usize) -> Result<Option<(String, RecordKind)>, Error> {
        let rec = match self.aln[nr].next_rec()? {
            Some(r) => r,
            None => return Ok(None),
        };

        let raw_name: &[u8] = rec
            .name()
            .map(|n| n.as_ref())
            .ok_or(Error::RecordHasNoReadName)?;
        let key = String::from_utf8_lossy(raw_name).to_string();

        let flags = rec.flags()?;
        let virtual_offset = self.record_counters[nr];
        self.record_counters[nr] += 1;

        // Secondary: virtual offset only -- no CIGAR/MD parsing needed.
        if flags.is_secondary() {
            return Ok(Some((
                key,
                RecordKind::Secondary {
                    flags,
                    virtual_offset,
                },
            )));
        }

        // Unmapped primary: never NW-scored -- skip CIGAR/MD/ref_len, keep qualities.
        if flags.is_unmapped() {
            let qualities: Vec<u8> = rec
                .quality_scores()
                .as_ref()
                .iter()
                .collect::<Result<Vec<u8>, std::io::Error>>()?;
            return Ok(Some((
                key,
                RecordKind::UnmappedPrimary {
                    flags,
                    virtual_offset,
                    qualities,
                },
            )));
        }

        // Mapped (primary or supplementary): full parse, as before.
        let ref_id = rec.ref_seq_id().transpose()?.unwrap_or(usize::MAX);
        let pos = rec
            .alignment_start()
            .transpose()?
            .map(|p| p.get())
            .unwrap_or(0);

        let mut ref_len = 0usize;
        let mut cigar_bytes = Vec::new();
        for op_result in rec.cigar().as_ref().iter() {
            let op = op_result?;
            match op.kind() {
                Kind::Match
                | Kind::Deletion
                | Kind::Skip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch => ref_len += op.len(),
                _ => {}
            }
            let code: u32 = match op.kind() {
                Kind::Match => 0,
                Kind::Insertion => 1,
                Kind::Deletion => 2,
                Kind::Skip => 3,
                Kind::SoftClip => 4,
                Kind::HardClip => 5,
                Kind::Pad => 6,
                Kind::SequenceMatch => 7,
                Kind::SequenceMismatch => 8,
            };
            let encoded = ((op.len() as u32) << 4) | code;
            cigar_bytes.extend_from_slice(&encoded.to_le_bytes());
        }

        let md = match rec.data().get(&Tag::MISMATCHED_POSITIONS).transpose()? {
            Some(Value::String(s)) => {
                let bytes: &[u8] = s.as_ref();
                bytes.to_vec()
            }
            _ => Vec::new(),
        };
        let qualities: Vec<u8> = rec
            .quality_scores()
            .iter()
            .collect::<Result<Vec<u8>, std::io::Error>>()?;
        let sequence: Vec<u8> = rec.sequence().iter().collect::<Vec<u8>>();
        let supp_count = match rec.data().get(&Tag::OTHER_ALIGNMENTS).transpose()? {
            Some(Value::String(s)) => {
                let bytes: &[u8] = s.as_ref();
                bytes.iter().filter(|&&c| c == b';').count()
            }
            _ => 0,
        };

        Ok(Some((
            key,
            RecordKind::Mapped(Box::new(MappedRecord {
                flags,
                ref_id,
                pos,
                ref_len,
                cigar_bytes,
                md,
                qualities,
                sequence,
                virtual_offset,
                supp_count,
            })),
        )))
    }

    fn ingest(&mut self, key: String, rec: RecordKind, nr: usize) -> Result<(), Error> {
        let bed = self.bed[nr].as_ref().filter(|b| !b.is_empty());
        let vcf = self.vcf[nr].as_ref().filter(|v| !v.is_empty());
        let (key, complete) = insert(
            &mut self.table,
            rec,
            key,
            nr,
            &mut self.seq_counter,
            bed,
            vcf,
        )?;
        if complete {
            let pending = self.table.remove(&key).unwrap();
            let seq_nr = pending.seq_nr;
            let sf = self.resolve_fragment(pending)?;
            self.staged.push(seq_nr, sf);
        }
        Ok(())
    }
    fn resolve_fragment(&mut self, pending: PendingFragment) -> Result<ScoredFragment, Error> {
        let supp_offsets = pending.supplementary_offsets;
        let driving_offsets = pending.driving.virtual_offsets();
        let lookup_offsets = pending.lookup.virtual_offsets();
        let dk = pending.driving.early_kind();
        let lk = pending.lookup.early_kind();

        // Both streams are Scoring -> run NW scoring on the retained MappedRecords.
        if dk.is_none() && lk.is_none() {
            let (
                StreamKind::Scoring {
                    records: recs_a,
                    mate_kinds: mk_a,
                },
                StreamKind::Scoring {
                    records: recs_b,
                    mate_kinds: mk_b,
                },
            ) = (pending.driving, pending.lookup)
            else {
                return Err(Error::BugExpectedScoringOnBothSides);
            };
            return self.evaluate_scoring_pair(
                *recs_a,
                *recs_b,
                mk_a,
                mk_b,
                driving_offsets,
                lookup_offsets,
                supp_offsets,
            );
        }

        use EarlyKind::*;
        let (winner_nr, decision, is_ambiguous) = match (dk, lk) {
            (Some(d1), Some(d2)) if d1 == d2 => (0, Decision::Ambiguous, true),
            (Some(AllPerfect), _) | (None, Some(AllUnmapped)) => (0, Decision::First, false),
            (Some(AllUnmapped), _) | (None, Some(AllPerfect)) => (1, Decision::Last, false),
            (None, None) => {
                unreachable!("Both streams are Scoring, should have been handled above")
            }
        };

        let drv = || driving_offsets.iter().map(|&o| (0, o));
        let lkp = || lookup_offsets.iter().map(|&o| (1, o));

        let (winner_offsets, loser_offsets) = if is_ambiguous {
            (drv().chain(lkp()).collect(), SmallVec::new())
        } else if winner_nr == 0 {
            (drv().collect(), lkp().collect())
        } else {
            (lkp().collect(), drv().collect())
        };

        Ok(ScoredFragment {
            winner_offsets,
            loser_offsets,
            supp_offsets,
            decision: self.add_decision_tag.then_some(decision),
            winner_nr,
            is_ambiguous,
        })
    }
    /// Run mate-level cancellation, then Tier 2.5 pre-assessment, then (if
    /// necessary) full NW scoring excluding any mates that cancelled.
    fn evaluate_scoring_pair(
        &mut self,
        recs_a: SmallVec<[RecordKind; 2]>,
        recs_b: SmallVec<[RecordKind; 2]>,
        mk_a: [Option<MateKind>; 2],
        mk_b: [Option<MateKind>; 2],
        offsets_a: SmallVec<[u64; 2]>,
        offsets_b: SmallVec<[u64; 2]>,
        supp_offsets: [SmallVec<[u64; 1]>; 2],
    ) -> Result<ScoredFragment, Error> {
        // -- Mate cancellation ----------------------------------------------
        // A mate slot cancels when both streams classify it identically as
        // Unmapped or Perfect: its per-base contribution is provably equal in
        // both streams (same read, same quality string), so it can be excluded
        // from NW scoring entirely without affecting the delta.
        let mut cancel_slot = [false, false];
        let mut all_present_cancel = true;
        let mut any_present = false;
        for i in 0..2 {
            match (mk_a[i], mk_b[i]) {
                (Some(a), Some(b)) if a == b && a != MateKind::Other => {
                    cancel_slot[i] = true;
                    any_present = true;
                }
                (None, None) => {}
                _ => {
                    all_present_cancel = false;
                    any_present = any_present || mk_a[i].is_some() || mk_b[i].is_some();
                }
            }
        }

        // Every present mate cancelled -> guaranteed tie, no NW needed at all.
        if all_present_cancel && any_present {
            let mut both: SmallVec<[(usize, u64); 2]> = offsets_a.iter().map(|&o| (0, o)).collect();
            both.extend(offsets_b.iter().map(|&o| (1, o)));
            return Ok(ScoredFragment {
                winner_offsets: both,
                loser_offsets: SmallVec::new(),
                supp_offsets,
                decision: self.add_decision_tag.then_some(Decision::Ambiguous),
                winner_nr: 0,
                is_ambiguous: true,
            });
        }

        let off_a: SmallVec<[(usize, u64); 2]> = offsets_a.iter().map(|&o| (0, o)).collect();
        let off_b: SmallVec<[(usize, u64); 2]> = offsets_b.iter().map(|&o| (1, o)).collect();
        let (winner_nr, decision, is_ambiguous) = match pre_assess_scoring_records(
            &recs_a,
            &recs_b,
            &self.penalties,
            self.ambiguous_log_threshold,
        ) {
            PreAssessResult::EarlyDecision(ord) => match ord {
                Ordering::Greater => (0, self.add_decision_tag.then_some(Decision::First), false),
                Ordering::Less => (1, self.add_decision_tag.then_some(Decision::Last), false),
                Ordering::Equal => (
                    0,
                    self.add_decision_tag.then_some(Decision::Ambiguous),
                    true,
                ),
            },
            PreAssessResult::FullScoring => {
                // Score only the mates that did NOT cancel, on both streams.
                let delta = self.nw_score_records(&recs_a, 0, cancel_slot)?
                    - self.nw_score_records(&recs_b, 1, cancel_slot)?;
                if delta > self.ambiguous_log_threshold {
                    (0, self.phred_delta(delta), false)
                } else if delta < -self.ambiguous_log_threshold {
                    (1, self.phred_delta(-delta), false)
                } else {
                    (
                        0,
                        self.add_decision_tag.then_some(Decision::Ambiguous),
                        true,
                    )
                }
            }
        };

        let (winner_offsets, loser_offsets) = if is_ambiguous {
            let mut both = off_a;
            both.extend(off_b);
            (both, SmallVec::new())
        } else if winner_nr == 0 {
            (off_a, off_b)
        } else {
            (off_b, off_a)
        };

        Ok(ScoredFragment {
            winner_offsets,
            loser_offsets,
            supp_offsets,
            decision,
            winner_nr,
            is_ambiguous,
        })
    }

    fn phred_delta(&self, abs_delta: f64) -> Option<Decision> {
        self.add_decision_tag.then(|| {
            let p = (10.0 * abs_delta / std::f64::consts::LN_10) as u32;
            Decision::PhredConfidence(p.min(255) as u8)
        })
    }

    fn nw_score_records(
        &mut self,
        records: &SmallVec<[RecordKind; 2]>,
        aln_idx: usize,
        cancel_slot: [bool; 2],
    ) -> Result<f64, Error> {
        let mut penalty = 0.0;
        let mut bufs: SmallVec<[RecordBuf; 2]> = SmallVec::new();

        for rk in records.iter() {
            let flags = rk.flags();
            if flags.is_secondary() {
                continue;
            }

            let slot = mate_slot(segment_id(&flags));
            if cancel_slot[slot] {
                continue;
            } // identical contribution in both streams -- skip

            let RecordKind::Mapped(m) = rk else { continue }; // Unmapped never scored

            if m.flags.is_supplementary() {
                penalty += self.penalties.chimeric_junction_penalty;
            }

            let mut buf = RecordBuf::default();
            *buf.flags_mut() = m.flags;
            *buf.reference_sequence_id_mut() = Some(m.ref_id);
            if m.pos > 0 {
                *buf.alignment_start_mut() =
                    Some(Position::new(m.pos).ok_or(Error::InvalidPosition(m.pos))?);
            }
            let mut cigar_ops = Vec::new();
            for chunk in m.cigar_bytes.chunks_exact(4) {
                let encoded = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
                let len = (encoded >> 4) as usize;
                let kind = match encoded & 0xF {
                    0 => Kind::Match,
                    1 => Kind::Insertion,
                    2 => Kind::Deletion,
                    3 => Kind::Skip,
                    4 => Kind::SoftClip,
                    5 => Kind::HardClip,
                    6 => Kind::Pad,
                    7 => Kind::SequenceMatch,
                    8 => Kind::SequenceMismatch,
                    k => return Err(Error::UnknownCigarOp(k)),
                };
                cigar_ops.push(Op::new(kind, len));
            }
            *buf.cigar_mut() = Cigar::from(cigar_ops);
            *buf.quality_scores_mut() = QualityScores::from_iter(m.qualities.iter().cloned());
            *buf.sequence_mut() = Sequence::from_iter(m.sequence.iter().cloned());
            let md_str = String::from_utf8(m.md.clone())?;
            let data: Data = [(Tag::MISMATCHED_POSITIONS, BufValue::from(md_str))]
                .into_iter()
                .collect();
            *buf.data_mut() = data;
            bufs.push(buf);
        }

        let flags_vec: SmallVec<[_; 2]> = bufs.iter().map(|b| b.flags()).collect();
        let mut mcfs: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();

        for (buf, flags) in bufs.iter().zip(flags_vec.iter()) {
            mcfs.push(MdCigFlags::try_from_record(buf, flags, self.bisulfite)?);
        }

        let seg: SmallVec<[&RecordBuf; READ_CT]> = bufs.iter().collect();

        let store = self.aln[aln_idx].variant_store();

        let mut dvnt: FragEvalVec<'_> = SmallVec::new();

        for buf in bufs.iter() {
            if buf.flags().is_unmapped() {
                dvnt.push(SmallVec::new());
            } else {
                let tid = buf.reference_sequence_id().ok_or(Error::NoRefSeqId)?;
                let start = buf.alignment_start().ok_or(Error::NoAlnStart)?.get();
                let end = start + buf.cigar().len();

                // Borrow from the `store` we fetched outside the loop
                let vars = if let Some(ref s) = store {
                    s.overlapping_multi(tid, start, end)
                } else {
                    SmallVec::new()
                };
                dvnt.push(vars);
            }
        }

        penalty += Fragment::new(&self.penalties, seg, mcfs)?
            .score(&mut self.scratch, &mut dvnt)
            .map_err(|e| Error::ScoreStreamError(aln_idx, e.to_string()))?;
        Ok(penalty)
    }

    fn emit_unmatched(&mut self, pending: PendingFragment) -> Result<(), Error> {
        let seq_nr = pending.seq_nr;
        let supp = pending.supplementary_offsets;
        let (driving_empty, _lookup_empty) =
            (pending.driving.is_empty(), pending.lookup.is_empty());
        let driving_offsets = pending.driving.virtual_offsets();
        let lookup_offsets = pending.lookup.virtual_offsets();

        let (winner_offsets, winner_nr): (SmallVec<[(usize, u64); 2]>, usize) = if !driving_empty {
            (driving_offsets.iter().map(|&o| (0, o)).collect(), 0)
        } else {
            (lookup_offsets.iter().map(|&o| (1, o)).collect(), 1)
        };

        self.staged.push(
            seq_nr,
            ScoredFragment {
                winner_offsets,
                loser_offsets: SmallVec::new(),
                supp_offsets: supp,
                decision: None,
                winner_nr,
                is_ambiguous: false,
            },
        );
        Ok(())
    }
}
