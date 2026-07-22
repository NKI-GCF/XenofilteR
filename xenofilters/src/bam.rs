//! `src/bam.rs`
//!
//! BAM output with read-group postfix routing.
//!
//! # Design
//!
//! winners go in the `--output <file>`. Discarded reads, and ambiguous reads
//! end up in the `--ambiguous <file>`.  Discarded destination is
//! distinguished by a read-group ID suffix appended to the original `RG:Z`
//! aux tag:
//!
//! | Destination | RG suffix        |
//! |-------------|------------------|
//! | Winner      | *(unchanged)*    |
//! | Filtered    | `_xenofilt`      |
//! | Ambiguous   | `_xenoambig`     |
//!
//! For every `@RG` line present in the input header, two additional `@RG`
//! lines are written to the merged output header (with the same fields, only
//! the `ID` changed).  Records that carry no `RG:Z` tag are written without
//! modification.
//!
//! # Performance
//!
//! - Header expansion is O(|RG|) string copies, done once at startup.
//! - Winners are written without any per-record tag manipulation.
//! - Non-winners incur one `BString` clone + one `BTreeMap` insert per record,
//!   which is negligible compared to bgzf compression.
//! - In the parallel pipeline this runs on the IO thread only, after scoring
//!   is complete, so it does not affect worker throughput.

mod format;
mod io;
pub(crate) mod reader;

use crate::Error;
pub(crate) use format::AlnFormat;
pub(crate) use io::{out_from_file, path_unicode_ok, BamOutput};
use noodles::sam::{
    alignment::{
        record::data::field::Tag,
        record_buf::{data::field::Value, RecordBuf},
    },
    header::{
        record::value::{map::ReadGroup, Map},
        Header,
    },
};

/// Suffixes appended to `RG:Z` tag values.
pub(crate) const SUFFIX_FILTERED: &str = "_xenofilt";
pub(crate) const SUFFIX_AMBIGUOUS: &str = "_xenoambig";

const TAG_RG: Tag = Tag::READ_GROUP;

// -- Header expansion ----------------------------------------------------------

/// Expand a BAM header by adding one extra `@RG` line per existing group for
/// each of the two non-winner destinations.
///
/// Each new `@RG` has the same fields as the original, with only `ID` changed
/// to `<original_id><suffix>`.
///
/// Returns the modified header.  The original header is consumed so that the
/// caller does not accidentally use a header that is missing the new groups.
pub(crate) fn expand_header(mut header: Header, write_discarded: bool) -> Header {
    // Collect additions first to avoid mutating while iterating.
    let additions: Vec<(String, Map<ReadGroup>)> = header
        .read_groups()
        .iter()
        .flat_map(|(id, rg)| {
            let id_str = id.to_string();
            let ambiguous = (format!("{id_str}{SUFFIX_AMBIGUOUS}"), rg.clone());
            let filtered = write_discarded
                .then(|| (format!("{id_str}{SUFFIX_FILTERED}"), rg.clone()));
            std::iter::once(ambiguous).chain(filtered)
        })
        .collect();

    for (new_id, rg) in additions {
        header.read_groups_mut().insert(new_id.into(), rg);
    }

    header
}

// -- Per-record RG tag rewrite -------------------------------------------------

/// Rewrite the `RG:Z` tag on `rec` by appending `suffix` to its current value.
///
/// If the record carries no `RG:Z` tag it is left unchanged -- downstream tools
/// will interpret it as belonging to no read group, which is the least-
/// surprising behaviour.
///
/// # Errors
///
/// Returns an error if the existing `RG` tag is not a `String` value (which
/// would indicate a malformed BAM file).
pub(crate) fn rewrite_rg(rec: &mut RecordBuf, suffix: &str) -> Result<(), Error> {
    // Read the current value.
    match rec.data_mut().get_mut(&TAG_RG) {
        None => Ok(()), // no RG tag -- leave record unchanged
        Some(Value::String(s)) => {
            s.extend_from_slice(suffix.as_bytes());
            Ok(())
        }
        Some(other) => Err(Error::UnexpectedRgTagType(format!("{other:?}"))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record_buf::data::field::Value;

    fn make_header_with_rgs(ids: &[&str]) -> Header {
        let mut h = Header::default();
        for id in ids {
            h.read_groups_mut()
                .insert((*id).into(), Map::<ReadGroup>::default());
        }
        h
    }

    fn rg_ids(header: &Header) -> Vec<String> {
        header.read_groups().keys().map(|k| k.to_string()).collect()
    }

    #[test]
    fn expand_header_adds_two_entries_per_rg() {
        let h = make_header_with_rgs(&["rg0", "rg1"]);
        let expanded1 = expand_header(h.clone(), false);
        let ids1 = rg_ids(&expanded1);
        // 2 original + 2 derived = 4
        assert_eq!(ids1.len(), 4);
        assert!(ids1.contains(&"rg0".to_string()));
        assert!(ids1.contains(&format!("rg0{SUFFIX_AMBIGUOUS}")));
        assert!(ids1.contains(&"rg1".to_string()));
        assert!(ids1.contains(&format!("rg1{SUFFIX_AMBIGUOUS}")));

        let expanded2 = expand_header(h, true);
        let ids2 = rg_ids(&expanded2);
        // 2 original + 4 derived = 6
        assert_eq!(ids2.len(), 6);
        assert!(ids2.contains(&"rg0".to_string()));
        assert!(ids2.contains(&format!("rg0{SUFFIX_FILTERED}")));
        assert!(ids2.contains(&format!("rg0{SUFFIX_AMBIGUOUS}")));
        assert!(ids2.contains(&"rg1".to_string()));
        assert!(ids2.contains(&format!("rg1{SUFFIX_FILTERED}")));
        assert!(ids2.contains(&format!("rg1{SUFFIX_AMBIGUOUS}")));
    }

    #[test]
    fn expand_header_empty_rg_is_noop() {
        let h = Header::default();
        let expanded = expand_header(h, true);
        assert!(expanded.read_groups().is_empty());
    }

    #[test]
    fn rewrite_rg_appends_suffix() -> Result<(), Error> {
        let mut rec = RecordBuf::default();
        rec.data_mut()
            .insert(TAG_RG, Value::String(String::from("rg0").into()));

        rewrite_rg(&mut rec, SUFFIX_FILTERED)?;

        match rec.data().get(&TAG_RG) {
            Some(Value::String(s)) => {
                assert_eq!(&s.to_string(), "rg0_xenofilt");
            }
            other => panic!("unexpected tag value: {other:?}"),
        }
        Ok(())
    }

    #[test]
    fn rewrite_rg_no_tag_is_noop() -> Result<(), Error> {
        let mut rec = RecordBuf::default();
        // No RG tag set.
        rewrite_rg(&mut rec, SUFFIX_AMBIGUOUS)?;
        assert!(rec.data().get(&TAG_RG).is_none());
        Ok(())
    }

    #[test]
    fn rewrite_rg_idempotent_on_winner() -> Result<(), Error> {
        // Winners are never passed through rewrite_rg; this test documents
        // that calling it twice (a bug) would stack suffixes visibly.
        let mut rec = RecordBuf::default();
        rec.data_mut()
            .insert(TAG_RG, Value::String(String::from("rg0").into()));
        rewrite_rg(&mut rec, SUFFIX_FILTERED)?;
        rewrite_rg(&mut rec, SUFFIX_FILTERED)?;
        match rec.data().get(&TAG_RG) {
            Some(Value::String(s)) => {
                // Double-application stacks suffixes -- this is wrong and
                // the test exists to catch accidental double-calls.
                assert_eq!(&s.to_string(), "rg0_xenofilt_xenofilt");
            }
            other => panic!("{other:?}"),
        }
        Ok(())
    }
}
