//! [`HashLookup`] — two-pass fragment-matching for position-sorted, indexed BAMs.
//!
//! **Pass 1** (sequential scan, lightweight):
//! - Reads `ScoringRecord`s (name, flags, ref_id, pos, CIGAR, MD, qualities,
//!   virtual_offset). No sequence.
//! - Inserts into a `NameTable` keyed by canonical read name.
//! - At fragment completion, classifies each stream as Early or NeedsScoring.
//! - Early fragments: virtual offsets stored, ScoringRecords dropped.
//! - NeedsScoring fragments: scored immediately using retained records.
//! - Completed fragments staged in `StagedOutput` for ordered emission.
//!
//! **Pass 2** (selective seek):
//! - For Early-assigned fragments: seeks to stored virtual offsets in the BAM
//!   file handle, reads full records, emits to output.
//! - For fully-scored fragments: full records were not retained; seeks via
//!   virtual offsets stored per-record during pass 1.
//! - Supplementary records: virtual offsets stored separately; retrieved in
//!   pass 2 following the fragment's decision.
//!
//! Single-threaded only. Output order follows driving-stream (stream 0) order.

pub(crate) mod assemble;
pub(crate) mod stage;
#[cfg(test)]
pub(crate) mod tests;

use crate::alignment::{stringify_record, Fragment, FragmentState, MdCigFlags, SimpleRec};
use crate::aln_stream::AlignmentStream;
use crate::config::{Config, StripReadSuffix};
use crate::filter_algorithm::collated::reader::canonical_name;
use crate::filter_algorithm::line_by_line::{ordering::Decision, Scratch, READ_CT};
use crate::penalty::Penalty;
use crate::region::{AmbiguousRegions, DiagnosticVariants};
use crate::variant::FragEvalVec;
use anyhow::{anyhow, Result};
use assemble::{insert, NameTable, PendingFragment, ScoringRecord, StreamKind};
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::SmallVec;
use stage::StagedOutput;
use std::cmp::Ordering;

// ---------------------------------------------------------------------------
// ScoredFragment — result of scoring, waiting for pass-2 emission
// ---------------------------------------------------------------------------

pub(crate) struct ScoredFragment {
    /// Virtual offsets of winning records (stream index + offset).
    pub(crate) winner_offsets: SmallVec<[(usize, u64); 2]>,
    /// Virtual offsets of losing records.
    pub(crate) loser_offsets: SmallVec<[(usize, u64); 2]>,
    /// Supplementary offsets — follow winner's decision.
    pub(crate) supp_offsets: [SmallVec<[u64; 1]>; 2],
    pub(crate) decision: Option<Decision>,
    /// Which stream index is the winner (for branch counter tracking).
    pub(crate) winner_nr: usize,
    /// True if result is ambiguous (both go to ambiguous output).
    pub(crate) is_ambiguous: bool,
}

// ---------------------------------------------------------------------------
// HashLookup
// ---------------------------------------------------------------------------

pub(crate) struct HashLookup<R: SimpleRec> {
    aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    table: NameTable,
    staged: StagedOutput,
    seq_counter: u64,
    penalties: Penalty,
    scratch: Scratch,
    pub(crate) branch_counters: [u64; 32],
    add_decision_tag: bool,
    ambiguous_log_threshold: f64,
    strip: StripReadSuffix,
    /// Per-stream ambiguous BED regions.
    bed: [Option<AmbiguousRegions>; 2],
    /// Per-stream diagnostic VCF sites.
    vcf: [Option<DiagnosticVariants>; 2],
}

impl<R: SimpleRec> HashLookup<R> {
    pub(crate) fn new(
        config: Config,
        mut aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
        bed: [Option<AmbiguousRegions>; 2],
        vcf: [Option<DiagnosticVariants>; 2],
    ) -> Result<Self> {
        assert_eq!(
            aln.len(),
            2,
            "HashLookup requires exactly 2 alignment streams"
        );
        let ambiguous_log_threshold = match config.ambiguous_threshold {
            0 => 0.0,
            t => (t as f64) * std::f64::consts::LN_10 / 10.0,
        };
        for (i, a) in aln.iter_mut().enumerate() {
            a.init_writers(&config, i)?;
        }
        Ok(Self {
            aln,
            table: NameTable::new(),
            staged: StagedOutput::new(),
            seq_counter: 0,
            penalties: config.to_penalties(),
            scratch: Scratch::new(),
            branch_counters: [0u64; 32],
            add_decision_tag: config.add_decision_tag,
            ambiguous_log_threshold,
            strip: config.strip_read_suffix,
            bed,
            vcf,
        })
    }

    pub(crate) fn process(&mut self) -> Result<()> {
        // Pass 1: interleave both streams.
        let mut exhausted = [false; 2];
        loop {
            let mut progress = false;
            for nr in 0..2usize {
                if exhausted[nr] {
                    continue;
                }
                // Read a lightweight scoring record from stream nr.
                match self.read_scoring_record(nr)? {
                    None => {
                        exhausted[nr] = true;
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
                &mut self.branch_counters,
                self.add_decision_tag,
            )?;
        }

        // Drain unmatched fragments.
        let unmatched: Vec<_> = self.table.drain().collect();
        for (_, pending) in unmatched {
            self.handle_unmatched(pending)?;
        }

        self.staged.flush_all(
            &mut self.aln,
            &mut self.branch_counters,
            self.add_decision_tag,
        )?;

        self.print_counters();
        Ok(())
    }

    /// Read one record from stream `nr` as a lightweight `ScoringRecord`.
    /// Returns the canonical name key and the record, or `None` at EOF.
    fn read_scoring_record(&mut self, nr: usize) -> Result<Option<(Box<[u8]>, ScoringRecord)>> {
        let rec = match self.aln[nr].next_rec()? {
            Some(r) => r,
            None => return Ok(None),
        };
        let raw_name = rec.name().map(|n| n.as_ref().to_vec()).unwrap_or_default();
        let key = canonical_name(&raw_name, self.strip);

        let flags = rec.flags()?;
        let ref_id = rec.ref_seq_id().transpose()?.unwrap_or(usize::MAX);
        let pos = rec
            .alignment_start()
            .transpose()?
            .map(|p| p.get())
            .unwrap_or(0);

        // Compute ref_len from CIGAR (count ops that consume reference).
        let ref_len = {
            use noodles::sam::alignment::record::cigar::op::Kind;
            let mut len = 0usize;
            for op in rec.cigar().as_ref().iter() {
                let op = op?;
                match op.kind() {
                    Kind::Match
                    | Kind::Deletion
                    | Kind::Skip
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch => len += op.len(),
                    _ => {}
                }
            }
            len
        };

        // Extract MD tag.
        let md = {
            use noodles::sam::alignment::record::data::field::{Tag, Value};
            match rec.data().get(&Tag::MISMATCHED_POSITIONS).transpose()? {
                Some(Value::String(s)) => s.as_ref().to_vec(),
                _ => Vec::new(),
            }
        };

        // Encode CIGAR as BAM bytes.
        let cigar_bytes = {
            use noodles::sam::alignment::record::cigar::op::Kind;
            let mut bytes = Vec::new();
            for op in rec.cigar().as_ref().iter() {
                let op = op?;
                let code: u8 = match op.kind() {
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
                let encoded = ((op.len() as u32) << 4) | (code as u32);
                bytes.extend_from_slice(&encoded.to_le_bytes());
            }
            bytes
        };

        let qualities: Vec<u8> = rec.quality_scores().as_ref().iter().copied().collect();

        // Virtual offset: not directly available from AlignmentStream trait.
        // We store a monotonic counter as a proxy; actual BGZF seek in pass 2
        // requires the BAM reader to expose virtual_position(). This is
        // implementation-dependent on the noodles BAM reader.
        // For now: store seq_counter as placeholder; pass-2 seek uses this.
        // TODO: expose virtual_position() from noodles BamReader when available.
        let virtual_offset = self.seq_counter;

        let srec = ScoringRecord {
            flags,
            ref_id,
            pos,
            ref_len,
            cigar_bytes,
            md,
            qualities,
            virtual_offset,
        };
        Ok(Some((key, srec)))
    }

    fn ingest(&mut self, key: Box<[u8]>, rec: ScoringRecord, nr: usize) -> Result<()> {
        let bed = self.bed[nr].as_ref();
        let vcf = self.vcf[nr].as_ref();
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

    fn resolve_fragment(&mut self, pending: PendingFragment) -> Result<ScoredFragment> {
        let supp_offsets = pending.supplementary_offsets;

        match (pending.driving, pending.lookup) {
            // Both early: ambiguous without scoring.
            (
                StreamKind::Early {
                    virtual_offsets: vo_a,
                    ..
                },
                StreamKind::Early {
                    virtual_offsets: vo_b,
                    ..
                },
            ) => {
                let dec = self.add_decision_tag.then_some(Decision::Ambiguous);
                let winner_offsets = vo_a.iter().map(|&o| (0, o)).collect();
                let loser_offsets = vo_b.iter().map(|&o| (1, o)).collect();
                Ok(ScoredFragment {
                    winner_offsets,
                    loser_offsets,
                    supp_offsets,
                    decision: dec,
                    winner_nr: 0,
                    is_ambiguous: true,
                })
            }
            // Driving early, lookup empty or worse → driving wins.
            (
                StreamKind::Early {
                    virtual_offsets: vo_a,
                    ..
                },
                lookup,
            ) => {
                let loser_offsets = match lookup {
                    StreamKind::Scoring { ref records, .. } => {
                        records.iter().map(|r| (1usize, r.virtual_offset)).collect()
                    }
                    StreamKind::Early {
                        ref virtual_offsets,
                        ..
                    } => virtual_offsets.iter().map(|&o| (1, o)).collect(),
                    StreamKind::Empty => SmallVec::new(),
                };
                let dec = self.add_decision_tag.then_some(Decision::First);
                Ok(ScoredFragment {
                    winner_offsets: vo_a.iter().map(|&o| (0, o)).collect(),
                    loser_offsets,
                    supp_offsets,
                    decision: dec,
                    winner_nr: 0,
                    is_ambiguous: false,
                })
            }
            // Lookup early, driving worse → lookup wins.
            (
                driving,
                StreamKind::Early {
                    virtual_offsets: vo_b,
                    ..
                },
            ) => {
                let loser_offsets = match driving {
                    StreamKind::Scoring { ref records, .. } => {
                        records.iter().map(|r| (0usize, r.virtual_offset)).collect()
                    }
                    StreamKind::Early {
                        ref virtual_offsets,
                        ..
                    } => virtual_offsets.iter().map(|&o| (0, o)).collect(),
                    StreamKind::Empty => SmallVec::new(),
                };
                let dec = self.add_decision_tag.then_some(Decision::Last);
                Ok(ScoredFragment {
                    winner_offsets: vo_b.iter().map(|&o| (1, o)).collect(),
                    loser_offsets,
                    supp_offsets,
                    decision: dec,
                    winner_nr: 1,
                    is_ambiguous: false,
                })
            }
            // Both need scoring.
            (
                StreamKind::Scoring {
                    records: recs_a, ..
                },
                StreamKind::Scoring {
                    records: recs_b, ..
                },
            ) => self.score_and_build(recs_a, recs_b, supp_offsets),
            // One empty (unmatched) — handled in handle_unmatched, not here.
            _ => Err(anyhow!("resolve_fragment called on incomplete fragment")),
        }
    }

    fn score_and_build(
        &mut self,
        recs_a: SmallVec<[ScoringRecord; 2]>,
        recs_b: SmallVec<[ScoringRecord; 2]>,
        supp_offsets: [SmallVec<[u64; 1]>; 2],
    ) -> Result<ScoredFragment> {
        // Build FragmentState equivalents from ScoringRecords for the scorer.
        // The scorer needs MdCigFlags; we reconstruct RecordBuf from scoring fields.
        // This is the full-scoring path — allocations here are expected.
        let offsets_a: SmallVec<[(usize, u64); 2]> =
            recs_a.iter().map(|r| (0, r.virtual_offset)).collect();
        let offsets_b: SmallVec<[(usize, u64); 2]> =
            recs_b.iter().map(|r| (1, r.virtual_offset)).collect();

        let score_a = self.score_records(&recs_a, 0)?;
        let score_b = self.score_records(&recs_b, 1)?;
        let delta = score_a - score_b;

        let (winner_offsets, loser_offsets, decision, winner_nr, is_ambiguous) =
            self.apply_delta_offsets(delta, offsets_a, offsets_b);

        Ok(ScoredFragment {
            winner_offsets,
            loser_offsets,
            supp_offsets,
            decision,
            winner_nr,
            is_ambiguous,
        })
    }

    fn apply_delta_offsets(
        &self,
        delta: f64,
        offsets_a: SmallVec<[(usize, u64); 2]>,
        offsets_b: SmallVec<[(usize, u64); 2]>,
    ) -> (
        SmallVec<[(usize, u64); 2]>,
        SmallVec<[(usize, u64); 2]>,
        Option<Decision>,
        usize,
        bool,
    ) {
        if delta > self.ambiguous_log_threshold {
            let dec = self.phred_delta(delta);
            (offsets_a, offsets_b, dec, 0, false)
        } else if delta < -self.ambiguous_log_threshold {
            let dec = self.phred_delta(-delta);
            (offsets_b, offsets_a, dec, 1, false)
        } else {
            let dec = self.add_decision_tag.then_some(Decision::Ambiguous);
            // Ambiguous: both are "winners", loser is empty.
            let mut combined = offsets_a;
            combined.extend(offsets_b);
            (combined, SmallVec::new(), dec, 0, true)
        }
    }

    fn phred_delta(&self, abs_delta: f64) -> Option<Decision> {
        self.add_decision_tag.then(|| {
            let p = (10.0 * abs_delta / std::f64::consts::LN_10) as u32;
            Decision::ConfDelta(p.min(255) as u8)
        })
    }

    /// Score a set of `ScoringRecord`s for one stream.
    /// Reconstructs minimal `RecordBuf`s sufficient for `Fragment::score`.
    fn score_records(
        &mut self,
        records: &SmallVec<[ScoringRecord; 2]>,
        aln_idx: usize,
    ) -> Result<f64> {
        use noodles::core::Position;
        use noodles::sam::alignment::record::data::field::{Tag, Value};
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};

        let mut bufs: SmallVec<[RecordBuf; 2]> = SmallVec::new();
        for sr in records.iter().filter(|r| r.is_primary()) {
            let mut buf = RecordBuf::default();
            *buf.flags_mut() = sr.flags;
            *buf.reference_sequence_id_mut() = Some(sr.ref_id);
            if sr.pos > 0 {
                *buf.alignment_start_mut() = Some(
                    Position::new(sr.pos).ok_or_else(|| anyhow!("Invalid position {}", sr.pos))?,
                );
            }
            // Decode CIGAR from BAM bytes.
            let mut cigar_ops = Vec::new();
            use noodles::sam::alignment::record::cigar::op::{Kind, Op};
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
                    k => return Err(anyhow!("Unknown CIGAR op {k}")),
                };
                cigar_ops.push(Op::new(kind, len));
            }
            *buf.cigar_mut() = Cigar::from(cigar_ops);
            *buf.quality_scores_mut() = QualityScores::from_iter(sr.qualities.iter().cloned());
            // Sequence: not stored — fill with N for the length implied by qualities.
            let seq_len = sr.qualities.len();
            *buf.sequence_mut() = Sequence::from(vec![b'N'; seq_len]);
            // MD tag.
            let data: noodles::sam::alignment::record_buf::Data =
                [(Tag::MISMATCHED_POSITIONS, Value::from(sr.md.as_slice()))]
                    .into_iter()
                    .collect();
            *buf.data_mut() = data;
            bufs.push(buf);
        }

        if bufs.is_empty() {
            return Ok(f64::NEG_INFINITY);
        }

        // Build MdCigFlags for each buf.
        let mut mcfs: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();
        for buf in bufs.iter() {
            let flags = buf.flags()?;
            mcfs.push(MdCigFlags::try_from_record(buf, &flags)?);
        }

        let seg: SmallVec<[&RecordBuf; READ_CT]> = bufs.iter().collect();
        let mut dvnt: FragEvalVec<'_> = SmallVec::new();
        for buf in bufs.iter() {
            if buf.flags()?.is_unmapped() {
                dvnt.push(SmallVec::new());
            } else {
                let tid = buf
                    .reference_sequence_id()
                    .ok_or_else(|| anyhow!("No ref seq id"))?;
                let start = buf
                    .alignment_start()
                    .transpose()?
                    .ok_or_else(|| anyhow!("No alignment start"))?
                    .get();
                let end = start + buf.cigar().len();
                let vars = if let Some(store) = self.aln[aln_idx].variant_store() {
                    store.overlapping_multi(tid, start, end)
                } else {
                    SmallVec::new()
                };
                dvnt.push(vars);
            }
        }

        Fragment::new(&self.penalties, seg, mcfs)?
            .score(&mut self.scratch, &mut dvnt)
            .map_err(|e| anyhow!("Score error stream {aln_idx}: {e}"))
    }

    fn handle_unmatched(&mut self, pending: PendingFragment) -> Result<()> {
        let seq_nr = pending.seq_nr;
        let supp = pending.supplementary_offsets;
        let (winner_offsets, winner_nr) = match (pending.driving, pending.lookup) {
            (
                StreamKind::Early {
                    virtual_offsets, ..
                },
                _,
            )
            | (StreamKind::Scoring { records: _, .. }, StreamKind::Empty) => {
                // driving has something, lookup empty
                // get offsets from driving
                match pending.driving {
                    StreamKind::Early {
                        virtual_offsets, ..
                    } => (virtual_offsets.iter().map(|&o| (0, o)).collect(), 0),
                    StreamKind::Scoring { records, .. } => {
                        (records.iter().map(|r| (0, r.virtual_offset)).collect(), 0)
                    }
                    _ => (SmallVec::new(), 0),
                }
            }
            (
                StreamKind::Empty,
                StreamKind::Early {
                    virtual_offsets, ..
                },
            ) => (virtual_offsets.iter().map(|&o| (1, o)).collect(), 1),
            (StreamKind::Empty, StreamKind::Scoring { records, .. }) => {
                (records.iter().map(|r| (1, r.virtual_offset)).collect(), 1)
            }
            _ => (SmallVec::new(), 0),
        };
        let sf = ScoredFragment {
            winner_offsets,
            loser_offsets: SmallVec::new(),
            supp_offsets: supp,
            decision: None,
            winner_nr,
            is_ambiguous: false,
        };
        self.staged.push(seq_nr, sf);
        Ok(())
    }

    pub(crate) fn print_counters(&self) {
        for i in 0..2 {
            eprintln!("hashlookup[filter:{}]: {}", i, self.branch_counters[i << 1]);
            eprintln!(
                "hashlookup[out:{}]: {}",
                i,
                self.branch_counters[1 + (i << 1)]
            );
            eprintln!("hashlookup[ambig:{}]: {}", i, self.branch_counters[16 + i]);
        }
    }
}
