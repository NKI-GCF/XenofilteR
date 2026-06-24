//! Read-coordinate-space profile for alignment comparison.
//!
//! Every physical read position `i ∈ [0, L)` is a shared index between the
//! two BAM streams: the same base was sequenced once and aligned to two
//! different references. This module builds per-read-base operation vectors
//! from the existing `MdCigFlags` abstraction and compares them to determine
//! whether a winner can be identified without (or with fewer) quality-score
//! lookups.
//!
//! # Deletions
//! Deletions consume reference bases but zero read positions. Their gap
//! penalty is quality-independent (`gap_open + len × gap_extend`), so they
//! are accumulated as counts and compared separately from per-base ops.
//!
//! # Insertions
//! Insertions of different lengths at the same read positions break the
//! per-position bijection. When either alignment contains an insertion the
//! comparison falls through to full Tier-3 scoring.
//!
//! # Orientation
//! `ScoreOpIter` (which underlies `build_read_profile`) already handles
//! reverse-strand MD orientation correctly, so no additional adjustment is
//! needed here.

use crate::alignment::{BaseOp, MdCigFlags, ScoreOpIter};
use smallvec::{Array, SmallVec};
use std::cmp::Ordering;

/// Operation at a single read position in 5′→3′ read space.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ReadOp {
    Match,
    Mismatch,
    SoftClip,
    Insertion,
}

/// Alignment summary in read coordinate space for one primary record.
#[derive(Debug)]
pub(crate) struct ReadProfile {
    /// Per-read-position operation (length == sequence length).
    pub(crate) ops: SmallVec<[ReadOp; 150]>,
    /// Number of deletion CIGAR events (each costs one `gap_open`).
    pub(crate) del_events: u32,
    /// Total deleted reference bases (total extra cost: `del_bases × gap_extend`).
    pub(crate) del_bases: u32,
    /// True when any insertion is present; triggers fall-through in comparison.
    pub(crate) has_insertions: bool,
    /// False when `ScoreOpIter` returned an error (malformed MD/CIGAR).
    pub(crate) valid: bool,
}

/// Build a `ReadProfile` by walking the existing `ScoreOpIter`.
/// Returns a profile with `valid = false` on any parse error;
/// callers must fall through to full scoring in that case.
pub(crate) fn build_read_profile(mcf: &MdCigFlags<'_>) -> ReadProfile {
    let mut ops: SmallVec<[ReadOp; 150]> = SmallVec::new();
    let mut del_events = 0u32;
    let mut del_bases = 0u32;
    let mut has_insertions = false;

    for op_result in ScoreOpIter::new(mcf) {
        match op_result {
            Err(_) => {
                return ReadProfile {
                    ops,
                    del_events,
                    del_bases,
                    has_insertions,
                    valid: false,
                };
            }
            Ok(BaseOp::Match) => ops.push(ReadOp::Match),
            Ok(BaseOp::Mis) => ops.push(ReadOp::Mismatch),
            Ok(BaseOp::Clip(n)) => {
                for _ in 0..n {
                    ops.push(ReadOp::SoftClip);
                }
            }
            Ok(BaseOp::Ins(n)) => {
                has_insertions = true;
                for _ in 0..n {
                    ops.push(ReadOp::Insertion);
                }
            }
            Ok(BaseOp::Del(n)) => {
                del_events += 1;
                del_bases += n as u32;
                // Deletions consume no read positions.
            }
            Ok(BaseOp::RefSkip(_)) => {
                // Intron skips: no read positions, no gap penalty in current model.
            }
        }
    }

    ReadProfile {
        ops,
        del_events,
        del_bases,
        has_insertions,
        valid: true,
    }
}

// ---------------------------------------------------------------------------
// Multi-segment (per-fragment) profiles
// ---------------------------------------------------------------------------

/// Per-fragment collection of per-mate profiles.
pub(crate) struct FragmentProfile {
    pub(crate) mates: SmallVec<[ReadProfile; 2]>,
}

impl FragmentProfile {
    pub(crate) fn valid(&self) -> bool {
        self.mates.iter().all(|p| p.valid)
    }
    pub(crate) fn has_insertions(&self) -> bool {
        self.mates.iter().any(|p| p.has_insertions)
    }
}

// ---------------------------------------------------------------------------
// Read-space comparison
// ---------------------------------------------------------------------------

/// Outcome of comparing two read-coordinate-space profiles.
#[derive(Debug)]
pub(crate) enum ReadSpaceDecision {
    /// All per-base positions where A ≠ B favour the same stream, and deletion
    /// counts are consistent with that direction.  No quality scores needed.
    EarlyDecision(Ordering), // Greater = A wins, Less = B wins, Equal = identical

    /// Per-base positions favouring A and B both exist; only these need
    /// quality-weighted scoring.  Deletion delta is pre-computed.
    PartialScoring {
        /// Positions (read index) where A = Match and B ≠ Match (A is better).
        a_better: SmallVec<[usize; 16]>,
        /// Positions where B = Match and A ≠ Match (B is better).
        b_better: SmallVec<[usize; 16]>,
        /// Signed counts for quality-independent deletion contribution.
        del: DelCounts,
    },

    /// Insertions present, malformed input, or mate-length mismatch.
    /// Fall through to full Tier-3 scoring.
    FallThrough,
}

/// Deletion event/base counts for both streams (quality-independent scoring).
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct DelCounts {
    pub(crate) a_events: i32,
    pub(crate) a_bases: i32,
    pub(crate) b_events: i32,
    pub(crate) b_bases: i32,
}

impl DelCounts {
    /// Signed delta: positive means A paid more deletion penalty (B better).
    pub(crate) fn delta(&self, gap_open: f64, gap_extend: f64) -> f64 {
        let a = (self.a_events as f64) * gap_open + (self.a_bases as f64) * gap_extend;
        let b = (self.b_events as f64) * gap_open + (self.b_bases as f64) * gap_extend;
        a - b // positive → A penalised more
    }

    /// Whether deletion counts alone suggest one stream is better.
    /// `None` when equal; `Some(Greater)` when A has fewer deletions.
    fn sign_only(&self) -> Option<Ordering> {
        let ord = self
            .a_events
            .cmp(&self.b_events)
            .then(self.a_bases.cmp(&self.b_bases));
        match ord {
            Ordering::Equal => None,
            Ordering::Less => Some(Ordering::Greater), // A has fewer
            Ordering::Greater => Some(Ordering::Less), // B has fewer
        }
    }
}

/// Compare single-mate profiles `a` and `b`.
///
/// Returns `FallThrough` immediately if either has insertions or if
/// the sequence lengths differ.
fn compare_mate_profiles(a: &ReadProfile, b: &ReadProfile) -> ReadSpaceDecision {
    if !a.valid || !b.valid || a.has_insertions || b.has_insertions {
        return ReadSpaceDecision::FallThrough;
    }
    if a.ops.len() != b.ops.len() {
        return ReadSpaceDecision::FallThrough;
    }

    let mut a_better: SmallVec<[usize; 16]> = SmallVec::new();
    let mut b_better: SmallVec<[usize; 16]> = SmallVec::new();

    for (i, (&op_a, &op_b)) in a.ops.iter().zip(b.ops.iter()).enumerate() {
        match (op_a, op_b) {
            // Same effective log-likelihood at this position → delta_i = 0.
            (ReadOp::Match, ReadOp::Match) => {}
            (ReadOp::Mismatch, ReadOp::Mismatch) => {}
            (ReadOp::SoftClip, ReadOp::SoftClip) => {}
            // Both score log_lik_mismatch[q_i]; same quality index → delta_i = 0.
            (ReadOp::Mismatch, ReadOp::SoftClip) => {}
            (ReadOp::SoftClip, ReadOp::Mismatch) => {}

            // A pays log_lik_match, B pays log_lik_mismatch → delta_i > 0 always.
            (ReadOp::Match, ReadOp::Mismatch) => a_better.push(i),
            (ReadOp::Match, ReadOp::SoftClip) => a_better.push(i),

            // B pays log_lik_match, A pays log_lik_mismatch → delta_i < 0 always.
            (ReadOp::Mismatch, ReadOp::Match) => b_better.push(i),
            (ReadOp::SoftClip, ReadOp::Match) => b_better.push(i),

            // Insertions filtered above; anything else is unexpected.
            _ => return ReadSpaceDecision::FallThrough,
        }
    }

    let del = DelCounts {
        a_events: a.del_events as i32,
        a_bases: a.del_bases as i32,
        b_events: b.del_events as i32,
        b_bases: b.del_bases as i32,
    };

    if a_better.is_empty() && b_better.is_empty() {
        // All per-base positions tie; winner determined by deletion counts alone.
        let ord = del.sign_only().unwrap_or(Ordering::Equal);
        return ReadSpaceDecision::EarlyDecision(ord);
    }

    if b_better.is_empty() {
        // All per-base differences favour A; check deletions are consistent.
        if del.a_events <= del.b_events && del.a_bases <= del.b_bases {
            return ReadSpaceDecision::EarlyDecision(Ordering::Greater);
        }
        // A wins on bases but B wins on deletions; need partial scoring.
    }

    if a_better.is_empty() {
        // All per-base differences favour B.
        if del.b_events <= del.a_events && del.b_bases <= del.a_bases {
            return ReadSpaceDecision::EarlyDecision(Ordering::Less);
        }
    }

    ReadSpaceDecision::PartialScoring {
        a_better,
        b_better,
        del,
    }
}

/// Compare two fragment-level profile collections.
///
/// For paired-end fragments, requires both mates to be present in the same
/// order in `a` and `b` (holds for namesorted and collated inputs).
/// Falls through if mate counts differ.
pub(crate) fn compare_fragment_profiles(
    a: &FragmentProfile,
    b: &FragmentProfile,
) -> ReadSpaceDecision {
    if a.mates.len() != b.mates.len() || a.mates.is_empty() {
        return ReadSpaceDecision::FallThrough;
    }
    if !a.valid() || !b.valid() || a.has_insertions() || b.has_insertions() {
        return ReadSpaceDecision::FallThrough;
    }

    let mut combined_a_better: SmallVec<[usize; 16]> = SmallVec::new();
    let mut combined_b_better: SmallVec<[usize; 16]> = SmallVec::new();
    let mut combined_del = DelCounts::default();
    // Read-position offset for second mate (so indices don't collide).
    let mut pos_offset = 0usize;

    for (ma, mb) in a.mates.iter().zip(b.mates.iter()) {
        match compare_mate_profiles(ma, mb) {
            ReadSpaceDecision::FallThrough => return ReadSpaceDecision::FallThrough,
            ReadSpaceDecision::EarlyDecision(Ordering::Equal) => {
                // This mate is a tie; continue to accumulate other mates.
                pos_offset += ma.ops.len();
                continue;
            }
            ReadSpaceDecision::EarlyDecision(ord) => {
                // This mate has a clear winner; check if it conflicts with prior mates.
                let a_wins = ord == Ordering::Greater;
                let no_conflict_a = combined_b_better.is_empty()
                    && (combined_del.b_events <= combined_del.a_events);
                let no_conflict_b = combined_a_better.is_empty()
                    && (combined_del.a_events <= combined_del.b_events);
                if (a_wins && !no_conflict_a) || (!a_wins && !no_conflict_b) {
                    // Conflict between mates; must score.
                    return ReadSpaceDecision::FallThrough;
                }
                // Consistent: represent as all-a_better or all-b_better sentinel
                // by adding a dummy position (the caller only checks emptiness).
                if a_wins {
                    combined_a_better.push(pos_offset);
                } else {
                    combined_b_better.push(pos_offset);
                }
                pos_offset += ma.ops.len();
            }
            ReadSpaceDecision::PartialScoring {
                a_better,
                b_better,
                del,
            } => {
                combined_a_better.extend(a_better.iter().map(|&p| p + pos_offset));
                combined_b_better.extend(b_better.iter().map(|&p| p + pos_offset));
                combined_del.a_events += del.a_events;
                combined_del.a_bases += del.a_bases;
                combined_del.b_events += del.b_events;
                combined_del.b_bases += del.b_bases;
                pos_offset += ma.ops.len();
            }
        }
    }

    // Final early-decision check after accumulating all mates.
    if combined_a_better.is_empty() && combined_b_better.is_empty() {
        let ord = combined_del.sign_only().unwrap_or(Ordering::Equal);
        return ReadSpaceDecision::EarlyDecision(ord);
    }
    if combined_b_better.is_empty()
        && combined_del.a_events <= combined_del.b_events
        && combined_del.a_bases <= combined_del.b_bases
    {
        return ReadSpaceDecision::EarlyDecision(Ordering::Greater);
    }
    if combined_a_better.is_empty()
        && combined_del.b_events <= combined_del.a_events
        && combined_del.b_bases <= combined_del.a_bases
    {
        return ReadSpaceDecision::EarlyDecision(Ordering::Less);
    }

    ReadSpaceDecision::PartialScoring {
        a_better: combined_a_better,
        b_better: combined_b_better,
        del: combined_del,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn profile(ops: &[ReadOp], del_events: u32, del_bases: u32) -> ReadProfile {
        ReadProfile {
            ops: ops.iter().copied().collect(),
            del_events,
            del_bases,
            has_insertions: false,
            valid: true,
        }
    }

    #[test]
    fn all_matching_is_equal() {
        let a = profile(&[ReadOp::Match, ReadOp::Match, ReadOp::Match], 0, 0);
        let b = profile(&[ReadOp::Match, ReadOp::Match, ReadOp::Match], 0, 0);
        assert!(matches!(
            compare_mate_profiles(&a, &b),
            ReadSpaceDecision::EarlyDecision(Ordering::Equal)
        ));
    }

    #[test]
    fn a_dominates_all_positions() {
        // A: match everywhere; B: one mismatch, one softclip.
        let a = profile(
            &[ReadOp::Match, ReadOp::Match, ReadOp::Match, ReadOp::Match],
            0,
            0,
        );
        let b = profile(
            &[
                ReadOp::Mismatch,
                ReadOp::Match,
                ReadOp::SoftClip,
                ReadOp::Match,
            ],
            0,
            0,
        );
        assert!(matches!(
            compare_mate_profiles(&a, &b),
            ReadSpaceDecision::EarlyDecision(Ordering::Greater)
        ));
    }

    #[test]
    fn b_dominates_all_positions() {
        let a = profile(&[ReadOp::SoftClip, ReadOp::Mismatch, ReadOp::Match], 0, 0);
        let b = profile(&[ReadOp::Match, ReadOp::Match, ReadOp::Match], 0, 0);
        assert!(matches!(
            compare_mate_profiles(&a, &b),
            ReadSpaceDecision::EarlyDecision(Ordering::Less)
        ));
    }

    #[test]
    fn mixed_positions_returns_partial_scoring() {
        let a = profile(&[ReadOp::Match, ReadOp::Mismatch, ReadOp::Match], 0, 0);
        let b = profile(&[ReadOp::Mismatch, ReadOp::Match, ReadOp::Match], 0, 0);
        match compare_mate_profiles(&a, &b) {
            ReadSpaceDecision::PartialScoring {
                a_better, b_better, ..
            } => {
                assert_eq!(a_better.as_slice(), &[0]);
                assert_eq!(b_better.as_slice(), &[1]);
            }
            other => panic!("expected PartialScoring, got {other:?}"),
        }
    }

    #[test]
    fn mismatch_vs_softclip_is_zero_delta() {
        // Both pay log_lik_mismatch[q_i]; delta = 0.
        let a = profile(&[ReadOp::Mismatch, ReadOp::Match], 0, 0);
        let b = profile(&[ReadOp::SoftClip, ReadOp::Match], 0, 0);
        // Position 0: Mismatch vs SoftClip → delta=0; position 1: Match vs Match → delta=0.
        assert!(matches!(
            compare_mate_profiles(&a, &b),
            ReadSpaceDecision::EarlyDecision(Ordering::Equal)
        ));
    }

    #[test]
    fn deletions_break_tie() {
        // Per-base positions are identical; A has a deletion, B does not.
        let a = profile(&[ReadOp::Match, ReadOp::Match], 1, 3);
        let b = profile(&[ReadOp::Match, ReadOp::Match], 0, 0);
        // A pays more deletion penalty → B better.
        assert!(matches!(
            compare_mate_profiles(&a, &b),
            ReadSpaceDecision::EarlyDecision(Ordering::Less)
        ));
    }

    #[test]
    fn insertion_forces_fall_through() {
        let a = ReadProfile {
            ops: smallvec::smallvec![ReadOp::Match, ReadOp::Insertion, ReadOp::Match],
            del_events: 0,
            del_bases: 0,
            has_insertions: true,
            valid: true,
        };
        let b = profile(&[ReadOp::Match, ReadOp::Match, ReadOp::Match], 0, 0);
        assert!(matches!(
            compare_mate_profiles(&a, &b),
            ReadSpaceDecision::FallThrough
        ));
    }

    #[test]
    fn a_better_on_bases_but_worse_on_deletions_is_partial() {
        // A wins all per-base positions but has more deletions → cannot early-decide.
        let a = profile(&[ReadOp::Match, ReadOp::Match], 2, 5);
        let b = profile(&[ReadOp::Mismatch, ReadOp::Mismatch], 0, 0);
        // Both b_better is empty, but del disfavours A → PartialScoring.
        assert!(matches!(
            compare_mate_profiles(&a, &b),
            ReadSpaceDecision::PartialScoring { .. }
        ));
    }
}
