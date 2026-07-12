// ============================================================================
// src/variant/name_to_id.rs
//
// Build a reference-sequence name → 0-based index map from a SAM header.
// Used by the indel equivalence expander and the diagnostic variant builder
// to resolve VCF CHROM strings to the integer IDs stored on variants.
// ============================================================================

use std::collections::HashMap;
use noodles::sam::Header;

/// Build a `String → usize` map from the `@SQ` lines of a SAM/BAM header.
///
/// The index is 0-based and matches the reference-sequence ID that noodles
/// returns from `Record::reference_sequence_id()`.
///
/// # Example
/// ```
/// // Header has @SQ SN:chr1 LN:... , @SQ SN:chr2 LN:... , ...
/// let m = header_name_to_id(&header);
/// assert_eq!(m["chr1"], 0);
/// assert_eq!(m["chr2"], 1);
/// ```
pub(crate) fn header_name_to_id(header: &Header) -> HashMap<String, usize> {
    header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_string(), i))
        .collect()
}

/// Build a `name → id` map directly from an indexed FASTA `.fai` file.
///
/// Used when no BAM header is available yet (e.g., during hashlookup
/// pre-scan before the first BAM record is read).
///
/// # Format
/// Each line of the .fai:  `name\tlength\toffset\tbpl\tbpl_with_nl`
pub(crate) fn fai_name_to_id(fai_path: &std::path::Path) -> std::io::Result<HashMap<String, usize>> {
    use std::io::{BufRead, BufReader};
    let f   = std::fs::File::open(fai_path)?;
    let mut map = HashMap::new();
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line  = line?;
        let name  = line.split('\t').next().unwrap_or("").to_string();
        if !name.is_empty() {
            map.insert(name, i);
        }
    }
    Ok(map)
}
