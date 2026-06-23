//!
//! Design principle: always accumulate `ScoringRecord`s during pass 1.
//! Early-assignability is evaluated at fragment *completion* time, not on
//! each individual record arrival. This avoids the downgrade problem (a
//! stream that looked assignable until its second primary arrived imperfect).
//!
//! At completion:
//! - If one stream's primaries are all perfect with no BED/VCF overlap →
//!   `StreamKind::Early`: only virtual offsets are kept; `ScoringRecord`s
//!   are dropped to free memory.
//! - Otherwise → `StreamKind::Scoring`: records retained for NW scoring.
//!
//! Supplementary records never affect assignment; their virtual offsets are
//! always stored separately for pass-2 retrieval.

use crate::region::{AmbiguousRegions, DiagnosticVariants};
use noodles::sam::alignment::record::Flags;
use smallvec::SmallVec;

// ---------------------------------------------------------------------------
// ScoringRecord
// ---------------------------------------------------------------------------

/// Minimum fields needed for pass-1 scoring and early-assignment checks.
/// No sequence field — CIGAR + MD + qualities suffice for NW scoring
/// and for deriving the reference/read base at diagnostic positions.
#[derive(Debug)]
pub(crate) struct ScoringRecord {
    pub(crate) flags: Flags,
    pub(crate) ref_id: usize,
    /// 1-based alignment start (BAM convention).
    pub(crate) pos: usize,
    /// Reference span derived from CIGAR.
    pub(crate) ref_len: usize,
    /// Raw BAM-encoded CIGAR bytes (4 bytes per op, little-endian u32).
    pub(crate) cigar_bytes: Vec<u8>,
    pub(crate) md: Vec<u8>,
    pub(crate) qualities: Vec<u8>,
    /// BGZF virtual offset — used in pass 2 for direct seek.
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
    /// True iff single-op CIGAR (all-Match) and MD is all digits (no mismatches).
    pub(crate) fn is_perfect(&self) -> bool {
        if self.is_unmapped() {
            return false;
        }
        // BAM CIGAR: each op is a u32le; low 4 bits = op code (0 = Match).
        // Single perfect op: exactly 4 bytes, op code 0.
        let cigar_ok = self.cigar_bytes.len() == 4 && (self.cigar_bytes[0] & 0x0F) == 0;
        let md_ok = !self.md.is_empty() && self.md.iter().all(|&b| b.is_ascii_digit());
        cigar_ok && md_ok
    }
}

// ---------------------------------------------------------------------------
// StreamBuf — accumulates records for one stream before classification
// ---------------------------------------------------------------------------

#[derive(Default)]
struct StreamBuf {
    records: SmallVec<[ScoringRecord; 2]>,
    primary_count: usize,
}

impl StreamBuf {
    fn push(&mut self, rec: ScoringRecord) {
        if rec.is_primary() {
            self.primary_count += 1;
        }
        self.records.push(rec);
    }

    /// Evaluate early-assignability once all primaries have arrived.
    /// Returns `StreamKind::Early` only if every primary is perfect and
    /// none overlaps an ambiguous BED or diagnostic VCF region.
    fn classify(self, bed: Option<&AmbiguousRegions>, vcf: Option<&DiagnosticVariants>) -> StreamKind {
        let all_assignable = self.records.iter().filter(|r| r.is_primary()).all(|r| {
            if !r.is_perfect() { return false; }
            if let Some(b) = bed
                && b.overlaps(r.ref_id, r.pos, r.pos + r.ref_len) {
                    return false;
                }
            if let Some(v) = vcf
                && v.overlaps(r.ref_id, r.pos, r.pos + r.ref_len) {
                    return false;
                }
            true
        });

        if all_assignable && self.primary_count > 0 {
            let offsets = self.records.iter().map(|r| r.virtual_offset).collect();
            StreamKind::Early { virtual_offsets: offsets, kind: EarlyKind::AllPerfect }
        } else {
            StreamKind::Scoring { records: Box::new(self.records) }
        }
    }
}

// ---------------------------------------------------------------------------
// StreamKind — post-classification state
// ---------------------------------------------------------------------------

pub(crate) enum EarlyKind {
    AllPerfect,
    AllUnmapped,
}

#[derive(Default)]
pub(crate) enum StreamKind {
    Early {
        virtual_offsets: SmallVec<[u64; 2]>,
        kind: EarlyKind,
    },
    Scoring {
        records: Box<SmallVec<[ScoringRecord; 2]>>,
    },
    #[default]
    Empty,
}

impl StreamKind {
    pub(crate) fn is_early(&self) -> bool {
        matches!(self, StreamKind::Early { .. })
    }
    pub(crate) fn is_empty(&self) -> bool {
        matches!(self, StreamKind::Empty)
    }
    pub(crate) fn virtual_offsets(&self) -> SmallVec<[u64; 2]> {
        match self {
            StreamKind::Early { virtual_offsets, .. } => virtual_offsets.clone(),
            StreamKind::Scoring { records, .. } => {
                records.iter().map(|r| r.virtual_offset).collect()
            }
            StreamKind::Empty => SmallVec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// PendingFragment
// ---------------------------------------------------------------------------

/// A fragment being assembled across both streams.
pub(crate) struct PendingFragment {
    /// Accumulation buffers — used until `classify()` is called.
    driving_buf: StreamBuf,
    lookup_buf: StreamBuf,
    /// Post-classification states — set after `classify()`.
    pub(crate) driving: StreamKind,
    pub(crate) lookup: StreamKind,
    /// Virtual offsets of supplementary records, per stream.
    pub(crate) supplementary_offsets: [SmallVec<[u64; 1]>; 2],
    /// Driving-stream insertion-order sequence number for staged output.
    pub(crate) seq_nr: u64,
    pub(crate) is_paired: Option<bool>,
}

impl PendingFragment {
    pub(crate) fn new(seq_nr: u64) -> Self {
        Self {
            driving_buf: StreamBuf::default(),
            lookup_buf: StreamBuf::default(),
            driving: StreamKind::Empty,
            lookup: StreamKind::Empty,
            supplementary_offsets: [SmallVec::new(), SmallVec::new()],
            seq_nr,
            is_paired: None,
        }
    }

    /// Push a record from stream `nr`. Returns `true` when the fragment is
    /// complete and ready for scoring or early assignment.
    pub(crate) fn push(
        &mut self,
        rec: ScoringRecord,
        nr: usize,
        bed: Option<&AmbiguousRegions>,
        vcf: Option<&DiagnosticVariants>,
    ) -> bool {
        if self.is_paired.is_none() {
            self.is_paired = Some(rec.flags.is_segmented());
        }
        if rec.is_supplementary() {
            self.supplementary_offsets[nr].push(rec.virtual_offset);
            return self.check_complete(bed, vcf);
        }
        if nr == 0 {
            self.driving_buf.push(rec);
        } else {
            self.lookup_buf.push(rec);
        }
        self.check_complete(bed, vcf)
    }

    fn expected_primaries(&self) -> usize {
        self.is_paired.map_or(1, |p| p as usize + 1)
    }

    /// Check completion and classify streams that have enough primaries.
    fn check_complete(
        &mut self,
        bed: Option<&AmbiguousRegions>,
        vcf: Option<&DiagnosticVariants>,
    ) -> bool {
        let exp = self.expected_primaries();

        // Classify driving stream if it has enough primaries and not yet classified.
        if self.driving.is_empty() && self.driving_buf.primary_count >= exp {
            let buf = std::mem::take(&mut self.driving_buf);
            self.driving = buf.classify(bed, vcf);
        }
        // Classify lookup stream similarly.
        if self.lookup.is_empty() && self.lookup_buf.primary_count >= exp {
            let buf = std::mem::take(&mut self.lookup_buf);
            self.lookup = buf.classify(bed, vcf);
        }

        self.is_complete_inner()
    }

    fn is_complete_inner(&self) -> bool {
        !self.driving.is_empty() && !self.lookup.is_empty()
    }
}

// ---------------------------------------------------------------------------
// NameTable
// ---------------------------------------------------------------------------

pub(crate) type NameTable = std::collections::HashMap<Box<[u8]>, PendingFragment>;

/// Insert `rec` from stream `nr` into `table`.
/// Returns the canonical key and `true` if the fragment is now complete.
pub(crate) fn insert(
    table: &mut NameTable,
    rec: ScoringRecord,
    canonical_name: Box<[u8]>,
    nr: usize,
    seq_counter: &mut u64,
    bed: Option<&AmbiguousRegions>,
    vcf: Option<&DiagnosticVariants>,
) -> (Box<[u8]>, bool) {
    let entry = table.entry(canonical_name.clone()).or_insert_with(|| {
        let sn = *seq_counter;
        *seq_counter += 1;
        PendingFragment::new(sn)
    });
    let complete = entry.push(rec, nr, bed, vcf);
    (canonical_name, complete)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::Flags;

    fn rec(flags_bits: u16, perfect: bool, ref_id: usize, pos: usize) -> ScoringRecord {
        // Perfect: single 10M op (BAM encoding: len=10, op=0 → 10<<4|0 = 0xA0 as u32le)
        // Imperfect: two ops (5M5S)
        let cigar_bytes = if perfect {
            // 10M: u32le = (10 << 4) | 0 = 160 = 0x000000A0
            vec![0xA0u8, 0x00, 0x00, 0x00]
        } else {
            // 5M5S: two ops
            vec![0x50u8, 0x00, 0x00, 0x00, 0x54u8, 0x00, 0x00, 0x00]
        };
        let md = if perfect {
            b"10".to_vec()
        } else {
            b"5".to_vec()
        };
        ScoringRecord {
            flags: Flags::from_bits(flags_bits).unwrap(),
            ref_id,
            pos,
            ref_len: 10,
            cigar_bytes,
            md,
            qualities: vec![30; 10],
            virtual_offset: pos as u64 * 1000,
        }
    }
    #[test]
    fn test_single_end_perfect_both_streams_classified() {
        let mut frag = PendingFragment::new(0);
        let complete = frag.push(rec(0, true, 0, 100), 0, None, None);
        assert!(!complete); // driving classified as Early but lookup still Empty
        let complete = frag.push(rec(0, true, 0, 200), 1, None, None);
        assert!(complete); // both Early → complete
        assert!(matches!(frag.driving, StreamKind::Early { .. }));
        assert!(matches!(frag.lookup, StreamKind::Early { .. }));
    }

    #[test]
    fn test_single_end_driving_perfect_lookup_imperfect() {
        let mut frag = PendingFragment::new(0);
        let complete = frag.push(rec(0, true, 0, 100), 0, None, None);
        assert!(!complete);
        let complete = frag.push(rec(0, false, 0, 200), 1, None, None);
        assert!(complete);
        assert!(matches!(frag.driving, StreamKind::Early { .. }));
        assert!(matches!(frag.lookup, StreamKind::Scoring { .. }));
    }

    #[test]
    fn test_driving_early_without_lookup() {
        // Both-stream requirement: driving alone is classified but fragment is not complete.
        let mut frag = PendingFragment::new(0);
        let complete = frag.push(rec(0, true, 0, 100), 0, None, None);
        assert!(!complete);
        assert!(frag.lookup.is_empty());
        assert!(matches!(&frag.driving, StreamKind::Early { .. }));
    }

    #[test]
    fn test_paired_end_requires_two_primaries() {
        let mut frag = PendingFragment::new(0);
        let complete = frag.push(rec(0x41, true, 0, 100), 0, None, None);
        assert!(!complete); // only 1 of 2 primaries for stream 0
        let complete = frag.push(rec(0x81, true, 0, 200), 0, None, None);
        assert!(!complete); // stream 0 classified (Early) but lookup still Empty
        let complete = frag.push(rec(0x41, true, 0, 300), 1, None, None);
        assert!(!complete); // stream 1 has only 1 of 2 expected primaries
        let complete = frag.push(rec(0x81, true, 0, 400), 1, None, None);
        assert!(complete); // both streams classified → complete
    }

    #[test]
    fn test_supplementary_does_not_count_as_primary() {
        let mut frag = PendingFragment::new(0);
        // Supplementary: flags 0x800
        let complete = frag.push(rec(0x800, true, 0, 500), 0, None, None);
        assert!(!complete);
        assert_eq!(frag.driving_buf.primary_count, 0);
    }

    #[test]
    fn test_ambiguous_region_forces_scoring() {
        use crate::region::ambiguous::{AmbiguousRegions, Region};
        use std::collections::HashMap;

        let mut per_ref = HashMap::new();
        per_ref.insert(0usize, vec![Region { start: 90, end: 110 }]);
        let bed = AmbiguousRegions::from_test(per_ref);

        let mut frag = PendingFragment::new(0);
        // pos=100, ref_len=10 → [100,110) overlaps [90,110)
        let complete = frag.push(rec(0, true, 0, 100), 0, Some(&bed), None);
        assert!(!complete);
        let complete = frag.push(rec(0, false, 0, 200), 1, Some(&bed), None);
        assert!(complete);
        // Driving overlaps ambiguous region → NeedsScoring despite being perfect.
        assert!(matches!(frag.driving, StreamKind::Scoring { .. }));
    }
}
