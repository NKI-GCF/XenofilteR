//! Trait for fast variant-aware early assignment bypass.
//!
//! When a fragment overlaps a population-specific variant with high allele
//! frequency in one species but near-zero in the other, full Tier-3 NW
//! scoring can be bypassed. This trait defines that evaluation interface.

use crate::filter_algorithm::line_by_line::ordering::Decision;

/// Implementors can evaluate variant evidence across streams and short-circuit
/// Tier-3 scoring when the evidence is deterministic.
///
/// Intended usage: called after `cmp_perfect` returns `None` but before NW DP.
#[allow(dead_code)]
pub(crate) trait AmbiguousVcfEvaluator {
    /// Evaluate whether variant evidence at `positions` (genomic coords) is
    /// sufficient to determine a winner without per-base alignment scoring.
    ///
    /// Returns `Some(Decision)` when evidence is conclusive, `None` to fall
    /// through to full NW scoring.
    ///
    /// # Implementation note
    /// `p_variant` must exceed 0.5 for a rescue delta to be positive.
    /// Implementors should enforce this invariant before returning `Some`.
    fn evaluate_ambiguous_variants(
        &self,
        stream0_positions: &[usize],
        stream1_positions: &[usize],
    ) -> Option<Decision>;
}
