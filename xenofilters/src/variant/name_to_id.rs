// ============================================================================
// src/variant/name_to_id.rs
//
// Build a reference-sequence name -> 0-based index map from a SAM header.
// Used by the indel equivalence expander and the diagnostic variant builder
// to resolve VCF CHROM strings to the integer IDs stored on variants.
// ============================================================================

use noodles::sam::Header;
use std::collections::HashMap;

/// Build a `String -> usize` map from the `@SQ` lines of a SAM/BAM header.
///
/// The index is 0-based and matches the reference-sequence ID that noodles
/// returns from `Record::reference_sequence_id()`.
pub(crate) fn header_name_to_id(header: &Header) -> HashMap<String, usize> {
    header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_string(), i))
        .collect()
}
