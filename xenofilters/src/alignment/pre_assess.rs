//! Tier 2.5 pre-assessment: unified structural + read-space alignment comparison.
//!
//! Two public entry points:
//!
//! * [`pre_assess_alignments`] — for `LineByLine` and `Collated`. Builds `ReadProfile`s
//!   **once** (single CIGAR+MD walk per record) then runs Tier 2.5a (aggregate subsumption)
//!   and Tier 2.5b (per-position read-space comparison) from the same data.
//!
//! * [`pre_assess_scoring_records`] — for `HashLookup`. Operates on raw BAM CIGAR bytes
//!   and MD byte slices stored in `ScoringRecord`; aggregate check only, avoiding the
//!   `RecordBuf` construction that happens later in `score_records`.

use super::read_profile::{
    build_read_profile, compare_fragment_profiles, FragmentProfile, ReadOp, ReadSpaceDecision,
};
use crate::alignment::MdCigFlags;
use crate::filter_algorithm::line_by_line::READ_CT;
use smallvec::SmallVec;
use std::cmp::Ordering;

// ---------------------------------------------------------------------------
// AlignSig — aggregate alignment quality signature (private to this module)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default)]
struct AlignSig {
    mismatches: usize,
    soft_clips: usize,
    indel_bases: usize,
}

impl AlignSig {
    fn perfect() -> Self {
        Self::default()
    }

    fn add(&mut self, other: AlignSig) {
        self.mismatches += other.mismatches;
        self.soft_clips += other.soft_clips;
        self.indel_bases += other.indel_bases;
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Derive an `AlignSig` from an already-built `FragmentProfile`.
/// No additional CIGAR/MD parsing required.
fn sig_from_fragment_profile(fp: &FragmentProfile) -> AlignSig {
    let mut sig = AlignSig::perfect();
    for mate in &fp.mates {
        for &op in &mate.ops {
            match op {
                ReadOp::Mismatch => sig.mismatches += 1,
                ReadOp::SoftClip => sig.soft_clips += 1,
                ReadOp::Insertion => sig.indel_bases += 1,
                ReadOp::Match => {}
            }
        }
        sig.indel_bases += mate.del_bases as usize;
    }
    sig
}

/// Compute an `AlignSig` from raw BAM-encoded CIGAR bytes + raw MD bytes.
/// Used by `HashLookup` where records are stored as `ScoringRecord`.
fn alignment_sig_raw(cigar_bytes: &[u8], md: &[u8]) -> AlignSig {
    let mut sig = AlignSig::perfect();
    for chunk in cigar_bytes.chunks_exact(4) {
        let enc = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        let op = enc & 0xF;
        let len = (enc >> 4) as usize;
        match op {
            4 | 5 => sig.soft_clips += len,  // SoftClip, HardClip
            1 | 2 => sig.indel_bases += len, // Insertion, Deletion
            _ => {}
        }
    }
    sig.mismatches = md_mismatches(md);
    sig
}

/// Count mismatch bases from a raw MD tag byte string.
/// Digit runs (matches) and `^`-prefixed deletion letters are skipped.
/// Never panics; returns 0 on empty or malformed input.
fn md_mismatches(md: &[u8]) -> usize {
    let mut count = 0;
    let mut i = 0;
    while i < md.len() {
        let b = md[i];
        if b.is_ascii_digit() {
            i += 1;
        } else if b == b'^' {
            i += 1;
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

/// `Some(Greater)` = `a` is better, `Some(Less)` = `b` is better,
/// `Some(Equal)` = identical metrics, `None` = incomparable (NW needed).
fn subsumes(a: &AlignSig, b: &AlignSig) -> Option<Ordering> {
    let a_dom = a.mismatches <= b.mismatches
        && a.soft_clips <= b.soft_clips
        && a.indel_bases <= b.indel_bases;
    let b_dom = b.mismatches <= a.mismatches
        && b.soft_clips <= a.soft_clips
        && b.indel_bases <= a.indel_bases;
    match (a_dom, b_dom) {
        (true, false) => Some(Ordering::Greater),
        (false, true) => Some(Ordering::Less),
        (true, true) => Some(Ordering::Equal),
        (false, false) => None,
    }
}

// ---------------------------------------------------------------------------
// PreAssessResult
// ---------------------------------------------------------------------------

/// Outcome of a pre-NW structural assessment.
#[derive(Debug)]
pub(crate) enum PreAssessResult {
    /// Structural dominance resolved the fragment; skip NW scoring.
    /// `Greater` = stream A wins, `Less` = stream B wins, `Equal` = tie.
    EarlyDecision(Ordering),
    /// No dominance found; fall through to full NW scoring.
    FullScoring,
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Unified Tier 2.5 pre-assessment for `LineByLine` and `Collated`.
///
/// Builds `ReadProfile`s from `MdCigFlags` **once** (a single `ScoreOpIter` walk
/// per primary record) and runs both tiers from the same data:
///
/// * **Tier 2.5a** — aggregate subsumption: derived from the profiles at zero
///   extra parsing cost. Fires when every alignment metric of one stream is ≤ the
///   other's on all three axes (mismatches, soft-clips, indel bases).
///
/// * **Tier 2.5b** — read-space comparison: if all per-read-position differences
///   favour the same stream and deletion counts are consistent, an early decision
///   is returned without NW DP.
///
/// Graceful fallback to `FullScoring` on segment-count mismatch, malformed MD/CIGAR,
/// insertions of differing length, or mixed-direction position evidence.
///
/// # Note on `PartialScoring`
/// When positions favouring A and positions favouring B are both present, full Tier-3
/// NW scoring is currently used. The `PartialScoring` position lists from Tier 2.5b
/// can be threaded into `Fragment::score` in a future refactor to skip identical bases.
pub(crate) fn pre_assess_alignments(
    mcfs_a: &SmallVec<[MdCigFlags<'_>; READ_CT]>,
    mcfs_b: &SmallVec<[MdCigFlags<'_>; READ_CT]>,
) -> PreAssessResult {
    if mcfs_a.len() != mcfs_b.len() || mcfs_a.is_empty() {
        return PreAssessResult::FullScoring;
    }

    // Single CIGAR+MD walk per record — both tiers consume these profiles.
    let fp_a = FragmentProfile {
        mates: mcfs_a.iter().map(build_read_profile).collect(),
    };
    let fp_b = FragmentProfile {
        mates: mcfs_b.iter().map(build_read_profile).collect(),
    };

    // Malformed MD/CIGAR: graceful fallback, no crash.
    if !fp_a.valid() || !fp_b.valid() {
        return PreAssessResult::FullScoring;
    }

    // Tier 2.5a: aggregate subsumption — free, derived from already-built profiles.
    let sig_a = sig_from_fragment_profile(&fp_a);
    let sig_b = sig_from_fragment_profile(&fp_b);
    if let Some(ord) = subsumes(&sig_a, &sig_b) {
        return PreAssessResult::EarlyDecision(ord);
    }

    // Tier 2.5b: read-space comparison — reuses same profiles, no extra I/O.
    match compare_fragment_profiles(&fp_a, &fp_b) {
        ReadSpaceDecision::EarlyDecision(ord) => PreAssessResult::EarlyDecision(ord),
        ReadSpaceDecision::PartialScoring { .. } | ReadSpaceDecision::FallThrough => {
            PreAssessResult::FullScoring
        }
    }
}

/// Tier 2.5 aggregate pre-assessment for `HashLookup`.
///
/// Operates on the raw BAM-encoded CIGAR bytes and MD byte slices stored in
/// `ScoringRecord`, short-circuiting before any `RecordBuf` construction.
/// Only the aggregate subsumption check (Tier 2.5a) is applied; building full
/// `ReadProfile`s from raw bytes requires per-base decoding that is better done
/// once inside `score_records` when full scoring is unavoidable.
pub(crate) fn pre_assess_scoring_records(
    recs_a: &[crate::filter_algorithm::hash_lookup::assemble::ScoringRecord],
    recs_b: &[crate::filter_algorithm::hash_lookup::assemble::ScoringRecord],
) -> PreAssessResult {
    let primaries_a: SmallVec<[_; 2]> = recs_a
        .iter()
        .filter(|r| r.is_primary() && !r.is_unmapped())
        .collect();
    let primaries_b: SmallVec<[_; 2]> = recs_b
        .iter()
        .filter(|r| r.is_primary() && !r.is_unmapped())
        .collect();

    if primaries_a.is_empty() || primaries_a.len() != primaries_b.len() {
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
        None => PreAssessResult::FullScoring,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn sig(mis: usize, clips: usize, indels: usize) -> AlignSig {
        AlignSig {
            mismatches: mis,
            soft_clips: clips,
            indel_bases: indels,
        }
    }

    #[test]
    fn perfect_subsumes_anything() {
        assert_eq!(
            subsumes(&sig(0, 0, 0), &sig(3, 2, 1)),
            Some(Ordering::Greater)
        );
    }

    #[test]
    fn equal_sigs_are_equal() {
        assert_eq!(
            subsumes(&sig(2, 1, 0), &sig(2, 1, 0)),
            Some(Ordering::Equal)
        );
    }

    #[test]
    fn incomparable_returns_none() {
        assert_eq!(subsumes(&sig(1, 5, 0), &sig(3, 1, 0)), None);
    }

    #[test]
    fn md_mismatches_counts_snvs_not_deletions() {
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
    fn alignment_sig_raw_clips_and_indels() {
        let cigar: Vec<u8> = {
            let mut v = Vec::new();
            v.extend_from_slice(&(5u32 << 4 | 4).to_le_bytes()); // 5S
            v.extend_from_slice(&(5u32 << 4 | 0).to_le_bytes()); // 5M
            v
        };
        let s = alignment_sig_raw(&cigar, b"5");
        assert_eq!(s.soft_clips, 5);
        assert_eq!(s.mismatches, 0);
        assert_eq!(s.indel_bases, 0);
    }
}
