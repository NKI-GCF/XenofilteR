use crate::alignment::FragmentState;
use crate::config::{Config, StripReadSuffix};
use crate::filter_algorithm::line_by_line::core::FragmentBuffer;
use crate::tests::create_record;
use crate::Error;
use crate::LineByLine;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::smallvec;

// %s/\vmock_rec\((b".*?")/create_record(\1, "10M", &[], &[], "10", false)?/g
#[test]
fn test_qname_suffix_logic() -> Result<(), Error> {
    let mut config = Config {
        strip_read_suffix: StripReadSuffix::Auto,
        ..Config::default()
    };

    let lbl = LineByLine::new(&config, smallvec![])?;
    let best: FragmentBuffer<RecordBuf> = smallvec![FragmentState::from_record(
        create_record(b"read/1", "10M", &[], &[], "10", false)?,
        0,
        false,
    )?];
    assert_eq!((lbl.is_new_qname)(&best, b"read/1"), Some(false));
    assert_eq!((lbl.is_new_qname)(&best, b"other/1"), Some(true));
    assert_eq!((lbl.is_new_qname)(&best, b"read/2"), Some(false));

    // Mode: Some(true) (Always strip last 2)
    config.strip_read_suffix = StripReadSuffix::True;
    let lbl = LineByLine::new(&config, smallvec![])?;
    assert_eq!((lbl.is_new_qname)(&best, b"read_suffix"), Some(true)); // "read" != "read_suff"

    // Mode: Some(false) (Exact match)
    config.strip_read_suffix = StripReadSuffix::False;
    let lbl = LineByLine::new(&config, smallvec![])?;
    assert_eq!((lbl.is_new_qname)(&best, b"read/1"), Some(false));
    assert_eq!((lbl.is_new_qname)(&best, b"read/2"), Some(true));
    Ok(())
}

#[test]
fn test_qname_suffix_logic_variable_mode() -> Result<(), Error> {
    let config = Config {
        strip_read_suffix: StripReadSuffix::Variable,
        ..Config::default()
    };
    let lbl = LineByLine::new(&config, smallvec![])?;

    let with_suffix: FragmentBuffer<RecordBuf> = smallvec![FragmentState::from_record(
        create_record(b"read/1", "10M", &[], &[], "10", false)?,
        0,
        false,
    )?];
    let no_suffix: FragmentBuffer<RecordBuf> = smallvec![FragmentState::from_record(
        create_record(b"read", "10M", &[], &[], "10", false)?,
        0,
        false,
    )?];

    assert_eq!((lbl.is_new_qname)(&with_suffix, b"read/2"), Some(false)); // suffix stripped
    assert_eq!((lbl.is_new_qname)(&with_suffix, b"other/2"), Some(true));
    assert_eq!((lbl.is_new_qname)(&no_suffix, b"read"), Some(false)); // exact match
    assert_eq!((lbl.is_new_qname)(&no_suffix, b"reae"), Some(true));
    Ok(())
}
