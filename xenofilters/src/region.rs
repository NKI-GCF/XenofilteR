//! Genomic region stores for early-assignment acceleration.
//!
//! - [`AmbiguousRegions`]: BED-derived intervals forcing full scoring.
//! - [`DiagnosticVariants`]: VCF-derived diagnostic positions forcing full scoring.
//! - [`early_assign`]: predicate combining both checks.
//! - [`tabix_query`]: random-access queries for the Collated algorithm.

pub(crate) mod ambiguous;
pub(crate) mod diagnostic;
pub(crate) mod load;
pub(crate) mod scored;
pub(crate) mod tabix_load;
pub(crate) mod tabix_query;

// AmbiguousRegions is now an alias; strand-awareness comes for free.
/// Ambiguous regions force Tier-2 fast-paths off when a read overlaps.
/// Strand field respected when present in BED column 6.
pub(crate) use ambiguous::AmbiguousRegions;
pub(crate) use diagnostic::DiagnosticVariants;
pub(crate) use scored::{ScoreFn, ScoredRegion, ScoredRegions, Strand};

/// Positive-score regions add a log-likelihood bonus to the stream's NW score.
/// Strand field respected when present in BED column 6.
pub(crate) type PositiveRegions = ScoredRegions;
