//! Genomic region stores for early-assignment acceleration.
//!
//! - [`AmbiguousRegions`]: BED-derived intervals forcing full scoring.
//! - [`SegregateVariants`]: VCF-derived diagnostic positions forcing full scoring.
//! - [`early_assign`]: predicate combining both checks.
//! - [`tabix_query`]: random-access queries for the Collated algorithm.

pub(crate) mod ambiguous;
pub(crate) mod diagnostic;
pub(crate) mod scored;

// AmbiguousRegions is now an alias; strand-awareness comes for free.
/// Ambiguous regions force Tier-2 fast-paths off when a read overlaps.
/// Strand field respected when present in BED column 6.
pub(crate) use scored::{ScoreFn, ScoredRegions};

/// Positive-score regions add a log-likelihood bonus to the stream's NW score.
/// Strand field respected when present in BED column 6.
pub(crate) type PositiveRegions = ScoredRegions;
