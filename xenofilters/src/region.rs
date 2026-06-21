//! Genomic region stores for early-assignment acceleration.
//!
//! - [`AmbiguousRegions`]: BED-derived intervals forcing full scoring.
//! - [`DiagnosticVariants`]: VCF-derived diagnostic positions forcing full scoring.
//! - [`early_assign`]: predicate combining both checks.
//! - [`tabix_query`]: random-access queries for the Collated algorithm.

pub(crate) mod ambiguous;
pub(crate) mod diagnostic;
pub(crate) mod early_assign;
pub(crate) mod tabix_query;

pub(crate) use ambiguous::AmbiguousRegions;
pub(crate) use diagnostic::DiagnosticVariants;
pub(crate) use early_assign::{check_early, EarlyCheck};
