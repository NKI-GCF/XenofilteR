//! Scoring helpers intentionally not used — the new backends (collated,
//! hash_lookup) each implement scoring inline to avoid lifetime conflicts
//! arising from the `AlignmentStream` trait combining reading and writing.
//!
//! This file is retained as a placeholder. If a clean abstraction becomes
//! possible in a future refactor it can be populated here.
