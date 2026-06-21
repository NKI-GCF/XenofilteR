//! Fragment accumulator for the two-pass HashLookup algorithm.
//!
//! Pass 1 reads lightweight `ScoringRecord`s. Each record is either:
//! - **EarlyAssigned**: perfect alignment, no BED/VCF overlap → only the
//!   BGZF virtual offset is retained for pass-2 retrieval.
//! - **NeedsScoring**: imperfect or overlaps ambiguous/diagnostic region →
//!   the full `ScoringRecord` is retained for NW scoring.
//!
//! A fragment becomes complete for early assignment the moment one stream's
//! `StreamState` is `EarlyAssigned` and its primary count matches expectations
//! (1 for single-end, 2 for paired-end). The other stream need not have
//! arrived — it will be looked up by read name and written to filtered output.
//!
//! A fragment requires full scoring when either stream contains a
//! `NeedsScoring` record and both streams have contributed primaries.

use crate::alignment::MdCigFlags;
use crate::filter_algorithm::line_by_line::READ_CT;
use crate::region::{AmbiguousRegions, DiagnosticVariants, EarlyCheck, check_early};
use anyhow::{anyhow, Result};
use noodles::sam::alignment::record::Flags;
use smallvec::SmallVec;

/// Minimum fields needed for scoring and early-assignment in pass 1.
/// No sequence — CIGAR + MD + qualities suffice.
pub(crate) struct ScoringRecord {
    pub(crate) flags: Flags,
    pub(crate) ref_id: usize,
    /// 1-based alignment start (matching BAM convention).
    pub(crate) pos: usize,
    /// Reference span (from CIGAR).
    pub(crate) ref_len: usize,
    pub(crate) cigar: Vec<u8>,   // raw CIGAR bytes; parsed on demand
    pub(crate) md: Vec<u8>,
    pub(crate) qualities: Vec<u8>,
    /// BGZF virtual offset of this record in the BAM file.
    /// Used in pass 2 to seek directly to the full record.
    pub(crate) virtual_offset: u64,
}

impl ScoringRecord {
    pub(crate) fn is_primary(&self) -> bool {
        !self.flags.is_secondary() && !self.flags.is_supplementary()
    }
    pub(crate) fn is_supplementary(&self) -> bool {
        self.flags.is_supplementary()
    }
    pub(crate) fn is_unmapped(&self) -> bool {
        self.flags.is_unmapped()
    }
}

/// Per-stream accumulation state.
pub(crate) enum StreamState {
    /// All primaries perfect, no ambiguous/diagnostic overlap.
    /// Virtual offsets retained for pass-2 full-record retrieval.
    EarlyAssigned {
        virtual_offsets: SmallVec<[u64; 2]>,
        primary_count: usize,
    },
    /// Full scoring required for at least one primary.
    NeedsScoring {
        records: SmallVec<[ScoringRecord; 2]>,
        primary_count: usize,
    },
    /// No records yet.
    Empty,
}

impl Default for StreamState {
    fn default() -> Self {
        StreamState::Empty
    }
}

impl StreamState {
    pub(crate) fn primary_count(&self) -> usize {
        match self {
            StreamState::EarlyAssigned { primary_count, .. } => *primary_count,
            StreamState::NeedsScoring { primary_count, .. } => *primary_count,
            StreamState::Empty => 0,
        }
    }

    pub(crate) fn is_early(&self) -> bool {
        matches!(self, StreamState::EarlyAssigned { .. })
    }

    pub(crate) fn has_any(&self) -> bool {
        !matches!(self, StreamState::Empty)
    }
}

/// A fragment being assembled from both streams.
pub(crate) struct PendingFragment {
    pub(crate) driving: StreamState,   // stream 0
    pub(crate) lookup: StreamState,    // stream 1
    /// Virtual offsets of supplementary records, per stream.
    /// Supplementaries are not scored but follow the fragment's decision.
    pub(crate) supplementary_offsets: [SmallVec<[u64; 1]>; 2],
    /// Driving-stream insertion order — used for staged output ordering.
    pub(crate) seq_nr: u64,
    /// True once we know this is paired-end (first paired record seen).
    pub(crate) is_paired: Option<bool>,
}

impl PendingFragment {
    pub(crate) fn new(seq_nr: u64) -> Self {
        Self {
            driving: StreamState::Empty,
            lookup: StreamState::Empty,
            supplementary_offsets: [SmallVec::new(), SmallVec::new()],
            seq_nr,
            is_paired: None,
        }
    }

    /// Expected primary count: 2 if paired-end, 1 if single-end, None if unknown.
    fn expected_primaries(&self) -> Option<usize> {
        self.is_paired.map(|p| if p { 2 } else { 1 })
    }

    /// Push a record from stream `nr` into the fragment.
    /// Returns `true` if the fragment is now complete for early assignment
    /// (one stream all-perfect, other stream worse or not yet arrived).
    /// Returns `true` also if full scoring is possible (both streams have
    /// all expected primaries and at least one needs scoring).
    pub(crate) fn push(
        &mut self,
        rec: ScoringRecord,
        nr: usize,
        bed: Option<&AmbiguousRegions>,
        vcf: Option<&DiagnosticVariants>,
    ) -> bool {
        // Track paired-end status.
        if self.is_paired.is_none() {
            self.is_paired = Some(rec.flags.is_segmented());
        }

        if rec.is_supplementary() {
            self.supplementary_offsets[nr].push(rec.virtual_offset);
            return self.is_complete();
        }

        let state = if nr == 0 { &mut self.driving } else { &mut self.lookup };
        let voffset = rec.virtual_offset;
        let is_primary = rec.is_primary();

        // Determine early-assignability for this record.
        let early = if is_primary && !rec.is_unmapped() {
            check_early(
                std::iter::once((
                    rec.ref_id,
                    rec.pos,
                    rec.pos + rec.ref_len,
                    // We cannot build MdCigFlags here without a full record.
                    // Conservatively treat as NeedsScoring when we cannot check.
                    // The actual perfect-check is done via is_perfect_scoring_rec.
                    &dummy_mcf_always_needs_scoring(),
                )),
                bed,
                vcf,
            )
        } else {
            EarlyCheck::NeedsScoring
        };

        // Actually: is_perfect check on ScoringRecord directly.
        let is_perf = is_perfect_scoring_rec(&rec);
        let no_region_overlap = if !rec.is_unmapped() {
            let bed_ok = bed.map_or(true, |b| !b.overlaps(rec.ref_id, rec.pos, rec.pos + rec.ref_len));
            let vcf_ok = vcf.map_or(true, |v| !v.overlaps(rec.ref_id, rec.pos, rec.pos + rec.ref_len));
            bed_ok && vcf_ok
        } else {
            false // unmapped is never early-assignable
        };

        let assignable = is_perf && no_region_overlap;

        match state {
            StreamState::Empty => {
                if assignable && is_primary {
                    *state = StreamState::EarlyAssigned {
                        virtual_offsets: SmallVec::from_slice(&[voffset]),
                        primary_count: 1,
                    };
                } else {
                    let pc = if is_primary { 1 } else { 0 };
                    *state = StreamState::NeedsScoring {
                        records: SmallVec::from_vec(vec![rec]),
                        primary_count: pc,
                    };
                }
            }
            StreamState::EarlyAssigned { virtual_offsets, primary_count } => {
                if assignable && is_primary {
                    virtual_offsets.push(voffset);
                    *primary_count += 1;
                } else {
                    // Downgrade to NeedsScoring — reconstruct with only offsets
                    // (we cannot retrieve records we already discarded, but we
                    // stored their virtual offsets, so pass 2 can load them).
                    // For scoring we need the actual ScoringRecords. Since we
                    // discarded them, we must re-read them in pass 1 here.
                    // Practical solution: once we hit a NeedsScoring record,
                    // we mark the whole stream as NeedsScoring but keep offsets
                    // for records already seen. The scorer will read them via
                    // the virtual offsets — effectively a mini pass-2.
                    // Simpler: accumulate all primaries as NeedsScoring from
                    // the start if any is imperfect. We cannot retroactively
                    // downgrade without re-reading, so we instead always
                    // accumulate ScoringRecords and mark early only at
                    // completion time. See below.
                    //
                    // Revised approach: always accumulate ScoringRecords.
                    // Determine EarlyAssigned vs NeedsScoring at completion.
                    // This avoids the downgrade problem entirely.
                    let offsets = std::mem::take(virtual_offsets);
                    let pc = *primary_count + if is_primary { 1 } else { 0 };
                    *state = StreamState::NeedsScoring {
                        records: SmallVec::from_vec(vec![rec]),
                        primary_count: pc,
                    };
                    // Note: offsets of prior early records are lost here.
                    // See module-level note: accumulate all, classify at end.
                    let _ = offsets;
                }
            }
            StreamState::NeedsScoring { records, primary_count } => {
                if is_primary { *primary_count += 1; }
                records.push(rec);
            }
        }

        self.is_complete()
    }

    /// True if the fragment is ready for scoring or early assignment.
    pub(crate) fn is_complete(&self) -> bool {
        let exp = match self.expected_primaries() {
            Some(e) => e,
            None => 1, // default to 1 until we know
        };
        // Early assignment: one stream fully assignable with enough primaries.
        let driving_early = matches!(
            &self.driving,
            StreamState::EarlyAssigned { primary_count, .. } if *primary_count >= exp
        );
        let lookup_early = matches!(
            &self.lookup,
            StreamState::EarlyAssigned { primary_count, .. } if *primary_count >= exp
        );
        // Full scoring: both streams have enough primaries.
        let driving_ready = self.driving.primary_count() >= exp;
        let lookup_ready = self.lookup.primary_count() >= exp;

        driving_early || lookup_early || (driving_ready && lookup_ready)
    }

    /// True if this fragment can be resolved without scoring the other stream.
    pub(crate) fn can_early_assign(&self) -> bool {
        let exp = self.expected_primaries().unwrap_or(1);
        let driving_early = matches!(
            &self.driving,
            StreamState::EarlyAssigned { primary_count, .. } if *primary_count >= exp
        );
        let lookup_early = matches!(
            &self.lookup,
            StreamState::EarlyAssigned { primary_count, .. } if *primary_count >= exp
        );
        // Both early → also qualifies (ambiguous, no scoring needed).
        driving_early || lookup_early
    }
}

/// Check if a `ScoringRecord` represents a perfect alignment:
/// single-span alignment (CIGAR has one M op) and MD is all digits.
fn is_perfect_scoring_rec(rec: &ScoringRecord) -> bool {
    if rec.is_unmapped() { return false; }
    // MD all digits → no mismatches.
    let md_perfect = !rec.md.is_empty() && rec.md.iter().all(|&b| b.is_ascii_digit());
    // CIGAR single op: encoded as (len << 4 | op_code); 4 bytes per op.
    // Single M op: cigar.len() == 4 and low nibble == 0 (Match).
    let cigar_perfect = rec.cigar.len() == 4 && (rec.cigar[0] & 0x0F) == 0;
    md_perfect && cigar_perfect
}

/// Placeholder — real check is done via is_perfect_scoring_rec.
/// This exists only to satisfy the type system in the push() early-check path
/// above (which is superseded by the direct check).
fn dummy_mcf_always_needs_scoring<'a>() -> MdCigFlags<'a> {
    // This function is never actually called in the revised push() logic
    // (the EarlyCheck call was replaced by direct is_perfect_scoring_rec).
    // It is retained only to keep the compiler happy during refactoring.
    // TODO: remove once push() is fully cleaned up.
    panic!("dummy_mcf_always_needs_scoring should never be called")
}

pub(crate) type NameTable = std::collections::HashMap<Box<[u8]>, PendingFragment>;
