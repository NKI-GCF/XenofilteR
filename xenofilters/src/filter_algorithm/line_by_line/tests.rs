use crate::filter_algorithm::line_by_line::core::AlnBuffer;
use crate::alignment::FragmentState;
use crate::tests::create_record;
use crate::config::{Config, StripReadSuffix};
use crate::LineByLine;
use anyhow::Result;
use smallvec::smallvec;
use noodles::sam::alignment::record_buf::RecordBuf;

// %s/\vmock_rec\((b".*?")/create_record(\1, "10M", &[], &[], "10", false)?/g
#[test]
fn test_qname_suffix_logic() -> Result<()> {
    let mut config = Config {
        strip_read_suffix: StripReadSuffix::Auto,
        ..Config::default()
    };

    let lbl = LineByLine::new(config.clone(), smallvec![])?;
    let best: AlnBuffer<RecordBuf> = smallvec![FragmentState::from_record(
        create_record(b"read/1", "10M", &[], &[], "10", false)?,
        0
    )?];
    assert_eq!((lbl.is_new_qname)(&best, b"read/1"), Some(false));
    assert_eq!((lbl.is_new_qname)(&best, b"other/1"), Some(true));
    assert_eq!((lbl.is_new_qname)(&best, b"read/2"), Some(false));

    // Mode: Some(true) (Always strip last 2)
    config.strip_read_suffix = StripReadSuffix::True;
    let lbl = LineByLine::new(config.clone(), smallvec![])?;
    assert_eq!((lbl.is_new_qname)(&best, b"read_suffix"), Some(true)); // "read" != "read_suff"

    // Mode: Some(false) (Exact match)
    config.strip_read_suffix = StripReadSuffix::False;
    let lbl = LineByLine::new(config, smallvec![])?;
    assert_eq!((lbl.is_new_qname)(&best, b"read/1"), Some(false));
    assert_eq!((lbl.is_new_qname)(&best, b"read/2"), Some(true));
    Ok(())
}

