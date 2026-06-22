//! Early-assignment predicate combining BED and VCF overlap checks.
//!
//! A fragment can be early-assigned when all its primary alignments are
//! perfect (single-op CIGAR, no MD mismatches) and none overlap an
//! ambiguous BED region or a diagnostic VCF position.

use crate::region::{AmbiguousRegions, DiagnosticVariants};

/// Result of the early-assignment check for one stream's set of primaries.
#[derive(Debug, PartialEq)]
pub(crate) enum EarlyCheck {
    Assignable,
    NeedsScoring,
}

/// Check whether `(ref_id, start, end, is_perfect)` tuples for all primaries
/// in a stream qualify for early assignment.
///
/// `is_perfect` must be derived from `ScoringRecord::is_perfect()` before
/// calling — the read sequence is not needed.
pub(crate) fn check_early(
    primaries: impl Iterator<Item = (usize, usize, usize, bool)>,
    bed: Option<&AmbiguousRegions>,
    vcf: Option<&DiagnosticVariants>,
) -> EarlyCheck {
    for (ref_id, start, end, is_perfect) in primaries {
        if !is_perfect {
            return EarlyCheck::NeedsScoring;
        }
        if let Some(b) = bed
            && b.overlaps(ref_id, start, end) {
                return EarlyCheck::NeedsScoring;
            }
        if let Some(v) = vcf
            && v.overlaps(ref_id, start, end) {
                return EarlyCheck::NeedsScoring;
            }
    }
    EarlyCheck::Assignable
}
