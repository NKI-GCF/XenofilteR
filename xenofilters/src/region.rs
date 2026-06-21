//! Genomic region stores for early-assignment acceleration.
//!
//! Two store types:
//! - [`AmbiguousRegions`]: BED-derived intervals where reads cannot be
//!   early-assigned even if the alignment is perfect. Forces full scoring.
//! - [`DiagnosticVariants`]: VCF-derived positions where a specific allele
//!   is diagnostic for one species. A read overlapping such a position must
//!   go through full scoring regardless of perfect-alignment status.
//!
//! Both are loaded fully into memory, sorted by position, and queried via
//! binary search (`partition_point`). For the Collated algorithm (name-ordered
//! BAM, no position guarantee) tabix-indexed files are queried instead.
//!
//! For LineByLine (namesorted) these stores are not used.

pub(crate) mod ambiguous;
pub(crate) mod diagnostic;
#[cfg(feature = "tabix")]
pub(crate) mod tabix_query;

pub(crate) use ambiguous::AmbiguousRegions;
pub(crate) use diagnostic::DiagnosticVariants;
