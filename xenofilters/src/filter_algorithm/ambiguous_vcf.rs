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
    /// `p_variant` must exceed 0.01 for a rescue delta to be positive.
    /// Implementors should enforce this invariant before returning `Some`.
    fn evaluate_ambiguous_variants(
        &self,
        stream0_positions: &[usize],
        stream1_positions: &[usize],
    ) -> Option<Decision>;
}


/// Conservative default: never short-circuits. Wired at the real call site
/// (see `collated::score_pair`) so the integration point is exercised and
/// tested, without asserting an unvalidated evidentiary rule for when
/// position-only overlap data is "sufficient" to skip NW scoring.
///
/// TODO(domain-owner): a real implementation needs to decide the rule --
/// e.g. "both positions fall in the SAME diagnostic VCF record, whose
/// stored p_variant exceeds some threshold, AND read bases at those
/// positions were independently confirmed (via MD) to equal ALT not just
/// REF" -- none of which is derivable from `&[usize]` alone. Extending the
/// trait to take `&[Eval]` (which carries allele identity via `Variant`)
/// rather than bare positions is the likely path forward.
pub(crate) struct NoopAmbiguousVcfEvaluator;

impl AmbiguousVcfEvaluator for NoopAmbiguousVcfEvaluator {
    fn evaluate_ambiguous_variants(&self, _: &[usize], _: &[usize]) -> Option<Decision> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn noop_evaluator_never_short_circuits() {
        let e = NoopAmbiguousVcfEvaluator;
        assert!(e.evaluate_ambiguous_variants(&[1, 2, 3], &[4, 5]).is_none());
    }
}
