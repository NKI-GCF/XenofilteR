//! Reserved for scoring abstractions shared across all three backends.
//!
//! Backend-specific scoring (LineByLine, HashLookup, Collated) lives in
//! each module. Cross-backend logic that does not require `AlignmentStream`
//! access — Tier 2.5 pre-assessment, read-space comparison, WIS scheduling
//! — lives in `crate::alignment::pre_assess` and `crate::alignment::read_profile`.
