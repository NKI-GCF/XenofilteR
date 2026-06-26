//! Tier 2.5 pre-assessment: unified structural + read-space alignment comparison.
//!
//! # Design
//!
//! The pre-assessment runs two sub-tiers from a SINGLE CIGAR+MD walk per record:
//!
//! **Tier 2.5a — match-count domination.**
//! `AlignSig` stores the count of exact-match bases (CIGAR M/=/X positions where
//! MD shows a digit, i.e. `BaseOp::Match` from `ScoreOpIter`) and the number of
//! pending supplementary alignments (from the `SA:Z` tag on each primary).
//! A higher match count is unambiguously better; fewer supplementaries is better.
//! A stream dominates when it is ≥ on matches AND ≤ on supplementaries.
//! This single-axis match count replaces the old 3-axis (mismatches / clips / indels)
//! approach: it is more efficient (one number instead of three comparisons) and more
//! directly correlated with the NW log-likelihood score.
//!
//! **Tier 2.5b — read-coordinate-space comparison.**
//! When Tier 2.5a cannot decide, `compare_fragment_profiles` identifies per-position
//! dominance in read space (positions where one stream matches and the other does not).
//!
//! Both tiers share the same `FragmentProfile` built by `build_read_profile`; no
//! second CIGAR+MD walk is ever required.

use super::read_profile::{
    build_read_profile, compare_fragment_profiles, FragmentProfile, ReadOp, ReadSpaceDecision,
};
use crate::alignment::MdCigFlags;
use crate::filter_algorithm::line_by_line::READ_CT;
use smallvec::SmallVec;
use std::cmp::Ordering;

// ---------------------------------------------------------------------------
// AlignSig — match-count-based quality signature (private to this module)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default)]
struct AlignSig {
    /// Exact-match bases across all PRIMARY segments (ScoreOpIter BaseOp::Match).
    /// Supplementary match bases are excluded: their gap penalty prevents
    /// safe cross-stream comparison without knowing the penalty magnitude.
    primary_match_bases: usize,
    /// Total supplementary alignments across all primary reads (SA:Z counts).
    /// Lower is better; used as the second comparison axis.
    supp_count: usize,
}

impl AlignSig {
    fn add(&mut self, other: AlignSig) {
        self.primary_match_bases += other.primary_match_bases;
        self.supp_count += other.supp_count;
    }
}

/// Derive a match-count `AlignSig` from an already-built `FragmentProfile`.
///
/// Primary records contribute their `ReadOp::Match` count and `supp_count`.
/// Supplementary records are NOT counted in `primary_match_bases` because
/// the chimeric-junction penalty (gap_open + bases × gap_extend) is applied
/// in the NW phase and may negate the match contribution; without knowing
/// the penalty magnitude we cannot compare supplementary matches across streams.
fn sig_from_fragment_profile(fp: &FragmentProfile) -> AlignSig {
    let mut sig = AlignSig::default();
    for mate in &fp.mates {
        if mate.is_supplementary {
            continue;
        }
        sig.primary_match_bases += mate.ops.iter().filter(|&&op| op == ReadOp::Match).count();
        sig.supp_count += mate.supp_count;
    }
    sig
}

/// Two-axis structural domination:
/// - more `primary_match_bases` is better (higher NW log-likelihood contribution)
/// - fewer `supp_count` is better (fewer chimeric-junction penalties)
///
/// Returns `None` when one axis favours A and the other favours B; the NW score
/// must then break the tie.
fn subsumes(a: &AlignSig, b: &AlignSig) -> Option<Ordering> {
    let a_dom = a.primary_match_bases >= b.primary_match_bases && a.supp_count <= b.supp_count;
    let b_dom = b.primary_match_bases >= a.primary_match_bases && b.supp_count <= a.supp_count;
    match (a_dom, b_dom) {
        (true, false) => Some(Ordering::Greater),
        (false, true) => Some(Ordering::Less),
        (true, true) => Some(Ordering::Equal), // identical on both axes
        (false, false) => None,                // e.g. A has more matches but more supps
    }
}

// ---------------------------------------------------------------------------
// Raw-bytes match counting for HashLookup (operates on ScoringRecord)
// ---------------------------------------------------------------------------

/// Count exact-match bases from raw BAM-encoded CIGAR bytes and a raw MD string.
///
/// Walks CIGAR M/=/X ops in parallel with the MD string. A base counts as a
/// match when the CIGAR op is M/=/X (0/7/8) AND the corresponding MD token is
/// a digit (not a mismatch letter). Gracefully returns the count accumulated so
/// far on any parse inconsistency.
fn match_count_raw(cigar_bytes: &[u8], md: &[u8]) -> usize {
    let mut matches = 0usize;
    let mut md_pos = 0usize;
    let mut md_match_remain = 0usize;

    for chunk in cigar_bytes.chunks_exact(4) {
        let enc = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        let op = enc & 0xF;
        let len = (enc >> 4) as usize;

        match op {
            // M (0), = (7), X (8): aligned bases — classify each against MD
            0 | 7 | 8 => {
                let mut cigar_remain = len;
                while cigar_remain > 0 {
                    if md_match_remain > 0 {
                        let take = cigar_remain.min(md_match_remain);
                        matches += take;
                        cigar_remain -= take;
                        md_match_remain -= take;
                    } else {
                        match md.get(md_pos) {
                            Some(&b) if b.is_ascii_digit() => {
                                let mut num = (b - b'0') as usize;
                                md_pos += 1;
                                while let Some(&d) = md.get(md_pos) {
                                    if !d.is_ascii_digit() {
                                        break;
                                    }
                                    num = num * 10 + (d - b'0') as usize;
                                    md_pos += 1;
                                }
                                md_match_remain = num;
                            }
                            Some(&b'A' | b'C' | b'G' | b'T' | b'N') => {
                                // Mismatch: advance MD and consume one CIGAR position
                                md_pos += 1;
                                cigar_remain -= 1;
                            }
                            Some(&b'^') => {
                                // Deletion block inside M run (malformed but graceful)
                                md_pos += 1;
                                while let Some(&d) = md.get(md_pos) {
                                    if d.is_ascii_digit() {
                                        break;
                                    }
                                    md_pos += 1;
                                }
                            }
                            _ => return matches, // malformed — stop here
                        }
                    }
                }
            }
            // D (2): deletion from reference — skip the ^letters MD block
            2 if md.get(md_pos) == Some(&b'^') => {
                md_pos += 1;
                while let Some(&d) = md.get(md_pos) {
                    if d.is_ascii_digit() {
                        break;
                    }
                    md_pos += 1;
                }
            }
            // I (1), N (3), S (4), H (5), P (6): no MD consumption, no matches
            _ => {}
        }
    }
    matches
}

// ---------------------------------------------------------------------------
// PreAssessResult
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub(crate) enum PreAssessResult {
    /// Structural dominance resolved the fragment without NW DP.
    /// `Greater` = stream A wins, `Less` = stream B wins, `Equal` = tie.
    EarlyDecision(Ordering),
    /// No structural dominance; fall through to full NW scoring.
    FullScoring,
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Unified Tier 2.5 pre-assessment for `LineByLine` and `Collated`.
///
/// Builds `ReadProfile`s once (single `ScoreOpIter` walk per record) and runs
/// both sub-tiers from the same data:
///
/// * **Tier 2.5a** — match-count domination via `AlignSig` (O(record) counter
///   increments derived for free from the already-built profile; no extra
///   CIGAR/MD parsing).
/// * **Tier 2.5b** — per-position read-space comparison via
///   `compare_fragment_profiles`.
///
/// Falls back to `FullScoring` on segment-count mismatch, malformed MD/CIGAR,
/// insertions, or when one stream has more matches but also more supplementaries.
pub(crate) fn pre_assess_alignments(
    mcfs_a: &SmallVec<[MdCigFlags<'_>; READ_CT]>,
    mcfs_b: &SmallVec<[MdCigFlags<'_>; READ_CT]>,
) -> PreAssessResult {
    if mcfs_a.len() != mcfs_b.len() || mcfs_a.is_empty() {
        return PreAssessResult::FullScoring;
    }

    // Single CIGAR+MD walk per record; both tiers read from these profiles.
    let fp_a = FragmentProfile {
        mates: mcfs_a.iter().map(build_read_profile).collect(),
    };
    let fp_b = FragmentProfile {
        mates: mcfs_b.iter().map(build_read_profile).collect(),
    };

    if !fp_a.valid() || !fp_b.valid() {
        return PreAssessResult::FullScoring;
    }

    // Tier 2.5a: match-count domination — derived from profiles at zero extra cost.
    let sig_a = sig_from_fragment_profile(&fp_a);
    let sig_b = sig_from_fragment_profile(&fp_b);
    if let Some(ord) = subsumes(&sig_a, &sig_b) {
        return PreAssessResult::EarlyDecision(ord);
    }

    // Tier 2.5b: per-position read-space comparison — reuses the same profiles.
    match compare_fragment_profiles(&fp_a, &fp_b) {
        ReadSpaceDecision::EarlyDecision(ord) => PreAssessResult::EarlyDecision(ord),
        ReadSpaceDecision::PartialScoring { .. } | ReadSpaceDecision::FallThrough => {
            PreAssessResult::FullScoring
        }
    }
}

/// Tier 2.5 pre-assessment for `HashLookup`.
///
/// Operates on raw BAM-encoded CIGAR bytes + MD byte slices from `ScoringRecord`,
/// avoiding `RecordBuf` construction until full NW is necessary. Uses
/// `match_count_raw` (same match-count semantics as Tier 2.5a above) and
/// respects `supp_count` from the SA:Z tag parsed in `next_scoring_record`.
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

    let mut sig_a = AlignSig::default();
    for r in &primaries_a {
        sig_a.primary_match_bases += match_count_raw(&r.cigar_bytes, &r.md);
        sig_a.supp_count += r.supp_count;
    }
    let mut sig_b = AlignSig::default();
    for r in &primaries_b {
        sig_b.primary_match_bases += match_count_raw(&r.cigar_bytes, &r.md);
        sig_b.supp_count += r.supp_count;
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

    fn sig(primary_matches: usize, supp_count: usize) -> AlignSig {
        AlignSig {
            primary_match_bases: primary_matches,
            supp_count,
        }
    }

    #[test]
    fn more_primary_matches_dominates() {
        assert_eq!(subsumes(&sig(90, 0), &sig(80, 0)), Some(Ordering::Greater));
    }

    #[test]
    fn fewer_primary_matches_loses() {
        assert_eq!(subsumes(&sig(70, 0), &sig(80, 0)), Some(Ordering::Less));
    }

    #[test]
    fn equal_matches_equal_supps_is_tie() {
        assert_eq!(subsumes(&sig(80, 0), &sig(80, 0)), Some(Ordering::Equal));
    }

    #[test]
    fn more_matches_but_more_supp_is_incomparable() {
        // A has more matches but also more supplementaries (unknown penalty) → None
        assert_eq!(subsumes(&sig(90, 2), &sig(80, 0)), None);
    }

    #[test]
    fn fewer_supp_wins_on_equal_matches() {
        // Same matches, A has no supplementaries → A wins (fewer penalty)
        assert_eq!(subsumes(&sig(80, 0), &sig(80, 1)), Some(Ordering::Greater));
    }

    #[test]
    fn more_matches_and_fewer_supp_dominates() {
        assert_eq!(subsumes(&sig(90, 0), &sig(80, 1)), Some(Ordering::Greater));
    }

    // --- match_count_raw ---

    fn cigar(ops: &[(u32, u32)]) -> Vec<u8> {
        // ops: [(op_code, length), ...]  op codes: M=0,I=1,D=2,S=4
        let mut v = Vec::new();
        for &(op, len) in ops {
            v.extend_from_slice(&((len << 4 | op).to_le_bytes()));
        }
        v
    }

    #[test]
    fn match_count_perfect_10m() {
        assert_eq!(match_count_raw(&cigar(&[(0, 10)]), b"10"), 10);
    }

    #[test]
    fn match_count_one_mismatch_in_5m() {
        // MD "2A2": 2 matches, mismatch, 2 matches → 4 matches
        assert_eq!(match_count_raw(&cigar(&[(0, 5)]), b"2A2"), 4);
    }

    #[test]
    fn match_count_soft_clip_not_counted() {
        // 3S5M, MD "5" → 5 matches (clips do not contribute)
        assert_eq!(match_count_raw(&cigar(&[(4, 3), (0, 5)]), b"5"), 5);
    }

    #[test]
    fn match_count_deletion_skipped() {
        // 5M3D5M, MD "5^ATG5" → 10 matches
        assert_eq!(
            match_count_raw(&cigar(&[(0, 5), (2, 3), (0, 5)]), b"5^ATG5"),
            10
        );
    }

    #[test]
    fn match_count_insertion_not_counted() {
        // 5M3I5M, MD "10" → 10 matches (3 insertion bases not in MD)
        assert_eq!(
            match_count_raw(&cigar(&[(0, 5), (1, 3), (0, 5)]), b"10"),
            10
        );
    }

    #[test]
    fn match_count_two_mismatches() {
        // 10M, MD "3A3C2" → 8 matches
        assert_eq!(match_count_raw(&cigar(&[(0, 10)]), b"3A3C2"), 8);
    }

    #[test]
    fn alignment_sig_raw_more_matches_wins() {
        // Perfect 10M vs 8M2-mismatch: sig_a > sig_b
        let a = sig(match_count_raw(&cigar(&[(0, 10)]), b"10"), 0);
        let b = sig(match_count_raw(&cigar(&[(0, 10)]), b"8AC"), 0);
        assert_eq!(a.primary_match_bases, 10);
        assert_eq!(b.primary_match_bases, 8);
        assert_eq!(subsumes(&a, &b), Some(Ordering::Greater));
    }
}
