//! [`HashLookup`] — two-pass fragment-matching for position-sorted BAMs.
//!
//! **Pass 1** (sequential scan, lightweight `ScoringRecord`s):
//! Reads name, flags, ref_id, pos, CIGAR, MD, qualities, virtual_offset.
//! No sequence. Inserts into a `FragmentTable`. At fragment completion, classifies
//! each stream as Early (all primaries perfect, no BED/VCF overlap) or
//! NeedsScoring. Early fragments retain only virtual offsets; Scoring
//! fragments retain `ScoringRecord`s for NW scoring.
//!
//! **Pass 2** (selective seek):
//! For each completed fragment, seeks to stored virtual offsets and reads
//! full records for output. Supplementary virtual offsets follow the
//! fragment's decision.
//!
//! Single-threaded only. Output order follows driving-stream (stream 0) order.

pub(crate) mod assemble;
pub(crate) mod stage;
#[cfg(test)]
pub(crate) mod tests;

use crate::alignment::{pre_assess_scoring_records, PreAssessResult};
use crate::alignment::{Fragment, MdCigFlags, SimpleRec};
use crate::aln_stream::AlignmentStream;
use crate::config::{Config, StripReadSuffix};
use crate::filter_algorithm::collated::reader::canonical_name;
use crate::filter_algorithm::line_by_line::{ordering::Decision, Scratch, READ_CT};
use crate::penalty::Penalty;
use crate::region::{AmbiguousRegions, DiagnosticVariants};
use crate::variant::FragEvalVec;
use assemble::{insert, EarlyKind, FragmentTable, PendingFragment, ScoringRecord, StreamKind};
use noodles::sam::alignment::record::Cigar;
use smallvec::SmallVec;
use stage::StagedOutput;
use std::cmp::Ordering;
use crate::Error;

// ---------------------------------------------------------------------------
// ScoredFragment
// ---------------------------------------------------------------------------

pub(crate) struct ScoredFragment {
    /// (stream_nr, virtual_offset) for winning records.
    pub(crate) winner_offsets: SmallVec<[(usize, u64); 2]>,
    /// (stream_nr, virtual_offset) for losing records.
    pub(crate) loser_offsets: SmallVec<[(usize, u64); 2]>,
    /// Supplementary offsets per stream — follow winner's decision.
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
    strip: StripReadSuffix,
    bed: [Option<AmbiguousRegions>; 2],
    vcf: [Option<DiagnosticVariants>; 2],
}

impl<R: SimpleRec> HashLookup<R> {
    pub(crate) fn new(
        config: Config,
        mut aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
        bed: [Option<AmbiguousRegions>; 2],
        vcf: [Option<DiagnosticVariants>; 2],
    ) -> Result<Self, Error> {
        if aln.len() != 2 {
            return Err(Error::AlgoRequiresTwoStreams);
        }
        let ambiguous_log_threshold = match config.ambiguous_threshold {
            0 => 0.0,
            t => (t as f64) * std::f64::consts::LN_10 / 10.0,
        };
        let aln_len = aln.len();
        for (i, a) in aln.iter_mut().enumerate() {
            a.init_writers(&config, i)?;
        }
        Ok(Self {
            aln,
            table: FragmentTable::new(),
            staged: StagedOutput::new(),
            seq_counter: 0,
            record_counters: [0, 0],
            penalties: config.to_penalties(),
            scratch: Scratch::new(),
            routing_counters: SmallVec::from_elem(0, aln_len * 4),
            add_decision_tag: config.add_decision_tag,
            ambiguous_log_threshold,
            strip: config.strip_read_suffix,
            bed,
            vcf,
        })
    }

    pub(crate) fn process(&mut self) -> Result<(), Error> {
        let mut exhausted = [false; 2];
        loop {
            let mut progress = false;
            for (nr, ex) in exhausted.iter_mut().enumerate() {
                if *ex {
                    continue;
                }
                match self.next_scoring_record(nr)? {
                    None => {
                        *ex = true;
                    }
                    Some((key, rec)) => {
                        progress = true;
                        self.ingest(key, rec, nr)?;
                    }
                }
            }
            if !progress {
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
        self.print_counters();
        Ok(())
    }

    fn next_scoring_record(&mut self, nr: usize) -> Result<Option<(Box<[u8]>, ScoringRecord)>, Error> {
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::data::field::{Tag, Value};

        let rec = match self.aln[nr].next_rec()? {
            Some(r) => r,
            None => return Ok(None),
        };

        let raw_name: Vec<u8> = rec
            .name()
            .map(|n| {
                let b: &[u8] = n.as_ref();
                b.to_vec()
            })
            .unwrap_or_default();
        let key = canonical_name(&raw_name, self.strip);

        let flags = rec.flags()?;
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
                let b: &[u8] = s.as_ref();
                b.to_vec()
            }
            _ => Vec::new(),
        };

        // Quality scores: noodles quality scores iterate as Result<u8, Error>.
        let qualities: Vec<u8> = rec
            .quality_scores()
            .as_ref()
            .iter()
            .collect::<Result<Vec<u8>, std::io::Error>>()?;

        // SA:Z: pending supplementary alignments for this read.
        // Each supplementary entry ends with ';', so semicolon count == count.
        let supp_count = match rec.data().get(&Tag::OTHER_ALIGNMENTS).transpose()? {
            Some(Value::String(s)) => {
                let b: &[u8] = s.as_ref();
                b.iter().filter(|&&c| c == b';').count()
            }
            _ => 0,
        };

        // Use a per-stream counter so fetch_by_virtual_offset(offset) correctly
        // indexes into stream[nr].original_reads[offset] for both real and mock streams.
        let virtual_offset = self.record_counters[nr];
        self.record_counters[nr] += 1;

        Ok(Some((
            key,
            ScoringRecord {
                flags,
                ref_id,
                pos,
                ref_len,
                cigar_bytes,
                md,
                qualities,
                virtual_offset,
                supp_count,
            },
        )))
    }

    fn ingest(&mut self, key: Box<[u8]>, rec: ScoringRecord, nr: usize) -> Result<(), Error> {
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
        );
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
        let dk = pending.driving.early_kind();
        let lk = pending.lookup.early_kind();
        let driving_offsets = pending.driving.virtual_offsets();
        let lookup_offsets = pending.lookup.virtual_offsets();

        // Consume scoring records before match (moving out of pending).
        let (driving_records, lookup_records) = match (pending.driving, pending.lookup) {
            (StreamKind::Scoring { records: ra }, StreamKind::Scoring { records: rb }) => {
                (Some(*ra), Some(*rb))
            }
            (StreamKind::Scoring { records: ra }, _) => (Some(*ra), None),
            (_, StreamKind::Scoring { records: rb }) => (None, Some(*rb)),
            _ => (None, None),
        };

        match (dk, lk) {
            // Both early-classified.
            (Some(d), Some(l)) => {
                use EarlyKind::*;
                match (d, l) {
                    (AllPerfect, AllPerfect) | (AllUnmapped, AllUnmapped) => {
                        let dec = self.add_decision_tag.then_some(Decision::Ambiguous);
                        let winner_offsets = driving_offsets
                            .iter()
                            .map(|&o| (0, o))
                            .chain(lookup_offsets.iter().map(|&o| (1, o)))
                            .collect();
                        Ok(ScoredFragment {
                            winner_offsets,
                            loser_offsets: SmallVec::new(),
                            supp_offsets,
                            decision: dec,
                            winner_nr: 0,
                            is_ambiguous: true,
                        })
                    }
                    (AllPerfect, AllUnmapped) => {
                        let dec = self.add_decision_tag.then_some(Decision::First);
                        Ok(ScoredFragment {
                            winner_offsets: driving_offsets.iter().map(|&o| (0, o)).collect(),
                            loser_offsets: lookup_offsets.iter().map(|&o| (1, o)).collect(),
                            supp_offsets,
                            decision: dec,
                            winner_nr: 0,
                            is_ambiguous: false,
                        })
                    }
                    (AllUnmapped, AllPerfect) => {
                        let dec = self.add_decision_tag.then_some(Decision::Last);
                        Ok(ScoredFragment {
                            winner_offsets: lookup_offsets.iter().map(|&o| (1, o)).collect(),
                            loser_offsets: driving_offsets.iter().map(|&o| (0, o)).collect(),
                            supp_offsets,
                            decision: dec,
                            winner_nr: 1,
                            is_ambiguous: false,
                        })
                    }
                }
            }
            // Driving early, lookup needs scoring.
            (Some(EarlyKind::AllPerfect), None) => {
                let dec = self.add_decision_tag.then_some(Decision::First);
                Ok(ScoredFragment {
                    winner_offsets: driving_offsets.iter().map(|&o| (0, o)).collect(),
                    loser_offsets: lookup_offsets.iter().map(|&o| (1, o)).collect(),
                    supp_offsets,
                    decision: dec,
                    winner_nr: 0,
                    is_ambiguous: false,
                })
            }
            (Some(EarlyKind::AllUnmapped), None) => {
                // Driving unmapped; Scoring stream has ≥1 mapped primary (else it would be AllUnmapped).
                let dec = self.add_decision_tag.then_some(Decision::Last);
                Ok(ScoredFragment {
                    winner_offsets: lookup_offsets.iter().map(|&o| (1, o)).collect(),
                    loser_offsets: driving_offsets.iter().map(|&o| (0, o)).collect(),
                    supp_offsets,
                    decision: dec,
                    winner_nr: 1,
                    is_ambiguous: false,
                })
            }
            // Lookup early, driving needs scoring.
            (None, Some(EarlyKind::AllPerfect)) => {
                let dec = self.add_decision_tag.then_some(Decision::Last);
                Ok(ScoredFragment {
                    winner_offsets: lookup_offsets.iter().map(|&o| (1, o)).collect(),
                    loser_offsets: driving_offsets.iter().map(|&o| (0, o)).collect(),
                    supp_offsets,
                    decision: dec,
                    winner_nr: 1,
                    is_ambiguous: false,
                })
            }
            (None, Some(EarlyKind::AllUnmapped)) => {
                let dec = self.add_decision_tag.then_some(Decision::First);
                Ok(ScoredFragment {
                    winner_offsets: driving_offsets.iter().map(|&o| (0, o)).collect(),
                    loser_offsets: lookup_offsets.iter().map(|&o| (1, o)).collect(),
                    supp_offsets,
                    decision: dec,
                    winner_nr: 0,
                    is_ambiguous: false,
                })
            }
            // Both need full scoring.
            (None, None) => {
                let recs_a = driving_records
                    .ok_or(Error::MissingDrivingRecords)?;
                let recs_b = lookup_records
                    .ok_or(Error::MissingLookupRecords)?;
                self.evaluate_scoring_pair(
                    recs_a,
                    recs_b,
                    driving_offsets,
                    lookup_offsets,
                    supp_offsets,
                )
            }
        }
    }

    /// Run Tier 2.5 pre-assessment then, if necessary, full NW scoring for a
    /// pair of `ScoringRecord` sets, and return a fully resolved `ScoredFragment`.
    ///
    /// Decision path:
    ///   1. `pre_assess_scoring_records` → `EarlyDecision` → return immediately.
    ///   2. Otherwise: `nw_score_records` on both sets → compare delta → route.
    fn evaluate_scoring_pair(
        &mut self,
        recs_a: SmallVec<[ScoringRecord; 2]>,
        recs_b: SmallVec<[ScoringRecord; 2]>,
        offsets_a: SmallVec<[u64; 2]>,
        offsets_b: SmallVec<[u64; 2]>,
        supp_offsets: [SmallVec<[u64; 1]>; 2],
    ) -> Result<ScoredFragment, Error> {
        // Tier 2.5: CIGAR/MD structural subsumption.
        // Borrows recs_a/recs_b immutably — both are still available for nw_score_records below.
        match pre_assess_scoring_records(&recs_a, &recs_b) {
            PreAssessResult::EarlyDecision(ord) => {
                let off_a: SmallVec<[(usize, u64); 2]> =
                    offsets_a.iter().map(|&o| (0, o)).collect();
                let off_b: SmallVec<[(usize, u64); 2]> =
                    offsets_b.iter().map(|&o| (1, o)).collect();
                return match ord {
                    Ordering::Greater => Ok(ScoredFragment {
                        winner_offsets: off_a,
                        loser_offsets: off_b,
                        supp_offsets,
                        decision: self.add_decision_tag.then_some(Decision::First),
                        winner_nr: 0,
                        is_ambiguous: false,
                    }),
                    Ordering::Less => Ok(ScoredFragment {
                        winner_offsets: off_b,
                        loser_offsets: off_a,
                        supp_offsets,
                        decision: self.add_decision_tag.then_some(Decision::Last),
                        winner_nr: 1,
                        is_ambiguous: false,
                    }),
                    Ordering::Equal => {
                        let mut both = off_a;
                        both.extend(off_b);
                        Ok(ScoredFragment {
                            winner_offsets: both,
                            loser_offsets: SmallVec::new(),
                            supp_offsets,
                            decision: self.add_decision_tag.then_some(Decision::Ambiguous),
                            winner_nr: 0,
                            is_ambiguous: true,
                        })
                    }
                };
            }
            PreAssessResult::FullScoring => {}
        }

        // Full NW scoring.
        let score_a = self.nw_score_records(&recs_a, 0)?;
        let score_b = self.nw_score_records(&recs_b, 1)?;
        let delta = score_a - score_b;

        let off_a: SmallVec<[(usize, u64); 2]> = offsets_a.iter().map(|&o| (0, o)).collect();
        let off_b: SmallVec<[(usize, u64); 2]> = offsets_b.iter().map(|&o| (1, o)).collect();

        if delta > self.ambiguous_log_threshold {
            let dec = self.phred_delta(delta);
            Ok(ScoredFragment {
                winner_offsets: off_a,
                loser_offsets: off_b,
                supp_offsets,
                decision: dec,
                winner_nr: 0,
                is_ambiguous: false,
            })
        } else if delta < -self.ambiguous_log_threshold {
            let dec = self.phred_delta(-delta);
            Ok(ScoredFragment {
                winner_offsets: off_b,
                loser_offsets: off_a,
                supp_offsets,
                decision: dec,
                winner_nr: 1,
                is_ambiguous: false,
            })
        } else {
            let dec = self.add_decision_tag.then_some(Decision::Ambiguous);
            let mut both = off_a;
            both.extend(off_b);
            Ok(ScoredFragment {
                winner_offsets: both,
                loser_offsets: SmallVec::new(),
                supp_offsets,
                decision: dec,
                winner_nr: 0,
                is_ambiguous: true,
            })
        }
    }

    fn phred_delta(&self, abs_delta: f64) -> Option<Decision> {
        self.add_decision_tag.then(|| {
            let p = (10.0 * abs_delta / std::f64::consts::LN_10) as u32;
            Decision::PhredConfidence(p.min(255) as u8)
        })
    }

    fn nw_score_records(
        &mut self,
        records: &SmallVec<[ScoringRecord; 2]>,
        aln_idx: usize,
    ) -> Result<f64, Error> {
        use noodles::core::Position;
        use noodles::sam::alignment::record::cigar::op::{Kind, Op};
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::{
            data::field::Value as BufValue, Cigar, Data, QualityScores, RecordBuf, Sequence,
        };

        let mut penalty = 0.0;

        let mut bufs: SmallVec<[RecordBuf; 2]> = SmallVec::new();
        for sr in records.iter() {
            if sr.flags.is_secondary() {
                continue;
            }
            if sr.is_supplementary() {
                penalty += self.penalties.chimeric_junction_penalty;
            }
            let mut buf = RecordBuf::default();
            *buf.flags_mut() = sr.flags;
            *buf.reference_sequence_id_mut() = Some(sr.ref_id);
            if sr.pos > 0 {
                *buf.alignment_start_mut() = Some(
                    Position::new(sr.pos).ok_or(Error::InvalidPosition(sr.pos))?,
                );
            }
            let mut cigar_ops = Vec::new();
            for chunk in sr.cigar_bytes.chunks_exact(4) {
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
            *buf.quality_scores_mut() = QualityScores::from_iter(sr.qualities.iter().cloned());
            let seq_len = sr.qualities.len();
            *buf.sequence_mut() = Sequence::from(vec![b'N'; seq_len]);
            let md_str = String::from_utf8(sr.md.clone())?;
            let data: Data = [(Tag::MISMATCHED_POSITIONS, BufValue::from(md_str))]
                .into_iter()
                .collect();
            *buf.data_mut() = data;
            bufs.push(buf);
        }
        let flags_vec: SmallVec<[_; 2]> = bufs.iter().map(|b| b.flags()).collect();
        let mut mcfs: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();

        for (buf, flags) in bufs.iter().zip(flags_vec.iter()) {
            mcfs.push(MdCigFlags::try_from_record(buf, flags)?);
        }

        let seg: SmallVec<[&RecordBuf; READ_CT]> = bufs.iter().collect();

        let store = self.aln[aln_idx].variant_store();

        let mut dvnt: FragEvalVec<'_> = SmallVec::new();

        for buf in bufs.iter() {
            if buf.flags().is_unmapped() {
                dvnt.push(SmallVec::new());
            } else {
                let tid = buf
                    .reference_sequence_id()
                    .ok_or(Error::NoRefSeqId)?;
                let start = buf
                    .alignment_start()
                    .ok_or(Error::NoAlignmentStart)?
                    .get();
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
    pub(crate) fn print_counters(&self) {
        let len = self.routing_counters.len();
        for nr in 0..(len / 4) {
            for (i, set) in ["discard", "out", "ambig"].iter().enumerate() {
                eprintln!(
                    "collated[{set}:{i}]: {}",
                    self.routing_counters[i + (nr * 4)]
                );
            }
        }
    }
}
