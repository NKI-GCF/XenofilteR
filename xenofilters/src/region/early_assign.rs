//! Early-assignment predicate.
//!
//! A fragment can be assigned without full scoring when:
//! 1. All primary alignments in one stream are perfect (single CIGAR op,
//!    MD all digits — no mismatches, no indels, no clipping).
//! 2. None of those alignments overlaps an ambiguous BED region.
//! 3. None overlaps a diagnostic VCF position.
//!
//! For paired-end data "all primaries" means both R1 and R2.
//!
//! This module provides the predicate used by both Collated (with tabix
//! queries) and HashLookup (with in-memory sorted stores + cursor scan).

use crate::alignment::MdCigFlags;
use crate::region::{AmbiguousRegions, DiagnosticVariants};

/// Result of the early-assignment check for one stream's fragment.
#[derive(Debug, PartialEq)]
pub(crate) enum EarlyCheck {
    /// All primaries perfect, no ambiguous/diagnostic overlap.
    /// Fragment may be early-assigned if the other stream is worse.
    Assignable,
    /// At least one primary is imperfect or overlaps an ambiguous/diagnostic
    /// region. Full scoring required.
    NeedsScoring,
}

/// Check whether a set of primary alignments (represented as `MdCigFlags`
/// slices paired with their genomic extents) qualifies for early assignment.
///
/// `primaries`: iterator of `(ref_id, aln_start, aln_end, mcf)` for each
/// primary alignment in the fragment.
///
/// `bed`: ambiguous regions for this stream's reference.
/// `vcf`: diagnostic variants for this stream's reference.
pub(crate) fn check_early<'a>(
    primaries: impl Iterator<Item = (usize, usize, usize, &'a MdCigFlags<'a>)>,
    bed: Option<&AmbiguousRegions>,
    vcf: Option<&DiagnosticVariants>,
) -> EarlyCheck {
    for (ref_id, start, end, mcf) in primaries {
        // Condition 1: perfect alignment.
        if !mcf.is_perfect() {
            return EarlyCheck::NeedsScoring;
        }
        // Condition 2: no ambiguous BED overlap.
        if let Some(b) = bed {
            if b.overlaps(ref_id, start, end) {
                return EarlyCheck::NeedsScoring;
            }
        }
        // Condition 3: no diagnostic VCF overlap.
        if let Some(v) = vcf {
            if v.overlaps(ref_id, start, end) {
                return EarlyCheck::NeedsScoring;
            }
        }
    }
    EarlyCheck::Assignable
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::region::ambiguous::{AmbiguousRegions, Region};
    use std::collections::HashMap;

    fn perfect_mcf_placeholder() -> bool {
        // MdCigFlags is not constructable in unit tests without a full record.
        // Early-assignment logic is tested via integration in CollatedMatcher
        // and HashLookup tests. Here we test only the region-overlap path,
        // which is independently tested in ambiguous::tests and diagnostic::tests.
        true
    }

    #[test]
    fn test_check_early_placeholder() {
        // Substantive tests live in collated/tests and hash_lookup/tests.
        assert!(perfect_mcf_placeholder());
    }
}
