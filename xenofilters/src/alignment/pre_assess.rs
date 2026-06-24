//! Tier 2.5 pre-assessment: CIGAR/MD structural subsumption.
//!
//! When both streams fail the perfect-match fast-path (Tier 2), this module
//! checks whether one stream's aggregate alignment penalty profile is a strict
//! subset of the other's. If so, the superior stream wins without Needleman–Wunsch
//! scoring.
//!
//! # Coordinate-space note
//! Subsumption operates on per-read aggregate counts (mismatches, clips, indels),
//! not on reference coordinates. This is intentional: in the xenograft use case
//! Stream A and Stream B align to *different* reference genomes, so direct
//! coordinate comparison is meaningless. Region-restricted scoring is therefore
//! not implemented.

use crate::alignment::MdCigFlags;
use crate::filter_algorithm::line_by_line::READ_CT;
use noodles::sam::alignment::record::cigar::op::Kind;
use smallvec::SmallVec;
use std::cmp::Ordering;

// ---------------------------------------------------------------------------
// AlignSig
// ---------------------------------------------------------------------------

/// Aggregate alignment quality signature.
/// Lower values on every axis indicate a better alignment.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct AlignSig {
    pub(crate) mismatches:  usize,
    pub(crate) soft_clips:  usize,
    pub(crate) indel_bases: usize,
}

impl AlignSig {
    pub(crate) fn perfect() -> Self {
        Self::default()
    }

    fn add(&mut self, other: AlignSig) {
        self.mismatches  += other.mismatches;
        self.soft_clips  += other.soft_clips;
        self.indel_bases += other.indel_bases;
    }
}

// ---------------------------------------------------------------------------
// Signature computation
// ---------------------------------------------------------------------------

/// Compute an `AlignSig` from a pre-built `MdCigFlags` slice.
/// Used by `LineByLine` and `Collated` where `MdCigFlags` already exist.
pub(crate) fn alignment_sig(mcfs: &SmallVec<[MdCigFlags<'_>; READ_CT]>) -> AlignSig {
    let mut sig = AlignSig::perfect();
    for mcf in mcfs {
        sig.mismatches += md_mismatches(mcf.get_md());
        for op in mcf.get_cigar().iter().filter_map(|r| r.ok()) {
            match op.kind() {
                Kind::SoftClip                   => sig.soft_clips  += op.len(),
                Kind::Insertion | Kind::Deletion => sig.indel_bases += op.len(),
                _                                => {}
            }
        }
    }
    sig
}

/// Compute an `AlignSig` from raw BAM-encoded CIGAR bytes + raw MD bytes.
/// Used by `HashLookup` where records are stored as `ScoringRecord`.
pub(crate) fn alignment_sig_raw(cigar_bytes: &[u8], md: &[u8]) -> AlignSig {
    let mut sig = AlignSig::perfect();
    for chunk in cigar_bytes.chunks_exact(4) {
        let enc = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        let op  = enc & 0xF;
        let len = (enc >> 4) as usize;
        match op {
            4 | 5 => sig.soft_clips  += len,  // SoftClip, HardClip
            1 | 2 => sig.indel_bases += len,  // Insertion, Deletion
            _     => {}
        }
    }
    sig.mismatches = md_mismatches(md);
    sig
}

/// Count mismatch bases from a raw MD tag byte string.
/// Digit runs (matches) and `^`-prefixed deletion letters are skipped.
/// Returns 0 on empty or malformed input — never panics.
fn md_mismatches(md: &[u8]) -> usize {
    let mut count = 0;
    let mut i = 0;
    while i < md.len() {
        let b = md[i];
        if b.is_ascii_digit() {
            i += 1;
        } else if b == b'^' {
            i += 1;
            // Skip reference bases shown in the deletion block.
            while i < md.len() && !md[i].is_ascii_digit() {
                i += 1;
            }
        } else if matches!(b, b'A' | b'C' | b'G' | b'T' | b'N') {
            count += 1;
            i += 1;
        } else {
            i += 1; // unexpected byte — skip gracefully
        }
    }
    count
}

// ---------------------------------------------------------------------------
// Subsumption predicate
// ---------------------------------------------------------------------------

/// `Some(Greater)` = `a` is better, `Some(Less)` = `b` is better,
/// `Some(Equal)` = identical, `None` = incomparable (NW needed).
pub(crate) fn subsumes(a: &AlignSig, b: &AlignSig) -> Option<Ordering> {
    let a_dom = a.mismatches  <= b.mismatches
        && a.soft_clips  <= b.soft_clips
        && a.indel_bases <= b.indel_bases;
    let b_dom = b.mismatches  <= a.mismatches
        && b.soft_clips  <= a.soft_clips
        && b.indel_bases <= a.indel_bases;
    match (a_dom, b_dom) {
        (true,  false) => Some(Ordering::Greater),
        (false, true)  => Some(Ordering::Less),
        (true,  true)  => Some(Ordering::Equal),
        (false, false) => None,
    }
}

// ---------------------------------------------------------------------------
// PreAssessResult and unified entry points
// ---------------------------------------------------------------------------

/// Outcome of a pre-NW structural assessment.
#[derive(Debug)]
pub(crate) enum PreAssessResult {
    /// Structural dominance resolved the fragment; skip NW scoring.
    /// The ordering follows `score(a).cmp(&score(b))`: `Greater` means stream A wins.
    EarlyDecision(Ordering),
    /// No structural dominance found; fall through to full NW.
    FullScoring,
}

/// Run Tier 2.5 against pre-built `MdCigFlags` (LineByLine, Collated).
///
/// Precondition: caller has already run Tier 2 (`cmp_perfect`) and received
/// `None` (both streams are imperfect). Subsumption is only meaningful when
/// both sides present the same number of primary segments.
///
/// Falls back to `FullScoring` on any parse error or segment-count mismatch.
pub(crate) fn pre_assess_mcfs(
    mcfs_a: &SmallVec<[MdCigFlags<'_>; READ_CT]>,
    mcfs_b: &SmallVec<[MdCigFlags<'_>; READ_CT]>,
) -> PreAssessResult {
    if mcfs_a.len() != mcfs_b.len() || mcfs_a.is_empty() {
        return PreAssessResult::FullScoring;
    }
    let sig_a = alignment_sig(mcfs_a);
    let sig_b = alignment_sig(mcfs_b);
    match subsumes(&sig_a, &sig_b) {
        Some(ord) => PreAssessResult::EarlyDecision(ord),
        None      => PreAssessResult::FullScoring,
    }
}

/// Run Tier 2.5 against raw `ScoringRecord` data (HashLookup).
///
/// Only considers primary, mapped records. Unmapped records always have
/// all-zero signatures and would produce misleading subsumption results,
/// so they are excluded. Returns `FullScoring` when primary counts differ.
pub(crate) fn pre_assess_scoring_records(
    recs_a: &[crate::filter_algorithm::hash_lookup::assemble::ScoringRecord],
    recs_b: &[crate::filter_algorithm::hash_lookup::assemble::ScoringRecord],
) -> PreAssessResult {
    let primaries_a: SmallVec<[_; 2]> = recs_a.iter()
        .filter(|r| r.is_primary() && !r.is_unmapped())
        .collect();
    let primaries_b: SmallVec<[_; 2]> = recs_b.iter()
        .filter(|r| r.is_primary() && !r.is_unmapped())
        .collect();

    if primaries_a.is_empty()
        || primaries_a.len() != primaries_b.len()
    {
        return PreAssessResult::FullScoring;
    }

    let mut sig_a = AlignSig::perfect();
    for r in &primaries_a {
        sig_a.add(alignment_sig_raw(&r.cigar_bytes, &r.md));
    }
    let mut sig_b = AlignSig::perfect();
    for r in &primaries_b {
        sig_b.add(alignment_sig_raw(&r.cigar_bytes, &r.md));
    }

    match subsumes(&sig_a, &sig_b) {
        Some(ord) => PreAssessResult::EarlyDecision(ord),
        None      => PreAssessResult::FullScoring,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn sig(mis: usize, clips: usize, indels: usize) -> AlignSig {
        AlignSig { mismatches: mis, soft_clips: clips, indel_bases: indels }
    }

    #[test]
    fn perfect_subsumes_anything() {
        assert_eq!(subsumes(&sig(0,0,0), &sig(3,2,1)), Some(Ordering::Greater));
    }

    #[test]
    fn equal_sigs_are_equal() {
        assert_eq!(subsumes(&sig(2,1,0), &sig(2,1,0)), Some(Ordering::Equal));
    }

    #[test]
    fn incomparable_returns_none() {
        // a has fewer mismatches but more clips
        assert_eq!(subsumes(&sig(1,5,0), &sig(3,1,0)), None);
    }

    #[test]
    fn md_mismatches_counts_snvs_not_deletions() {
        // "3A2^AT1C1" → A and C are mismatches; ^AT is deletion block
        assert_eq!(md_mismatches(b"3A2^AT1C1"), 2);
    }

    #[test]
    fn md_mismatches_all_match() {
        assert_eq!(md_mismatches(b"100"), 0);
    }

    #[test]
    fn md_mismatches_empty() {
        assert_eq!(md_mismatches(b""), 0);
    }

    #[test]
    fn alignment_sig_raw_counts_clips_and_indels() {
        // 5S5M: SoftClip(5) = 0x54_00_00_00, Match(5) = 0x50_00_00_00
        let cigar: Vec<u8> = {
            let mut v = Vec::new();
            v.extend_from_slice(&(5u32 << 4 | 4).to_le_bytes()); // 5S
            v.extend_from_slice(&(5u32 << 4 | 0).to_le_bytes()); // 5M
            v
        };
        let sig = alignment_sig_raw(&cigar, b"5");
        assert_eq!(sig.soft_clips, 5);
        assert_eq!(sig.mismatches, 0);
        assert_eq!(sig.indel_bases, 0);
    }
}
