use crate::aln_stream::tests::MockStream;
use crate::aln_stream::AlignmentStream;
use crate::config::{Config, StripReadSuffix};
use crate::filter_algorithm::hash_lookup::HashLookup;
use crate::tests::create_record;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::smallvec;
use crate::Error;

fn config() -> Config {
    Config {
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        strip_read_suffix: StripReadSuffix::False,
        ..Config::default()
    }
}

fn make_lookup(
    stream0: Vec<RecordBuf>,
    stream1: Vec<RecordBuf>,
    cfg: Config,
) -> Result<HashLookup<RecordBuf>, Error> {
    let s0 = Box::new(MockStream::new(0, stream0)) as Box<dyn AlignmentStream<RecordBuf>>;
    let s1 = Box::new(MockStream::new(1, stream1)) as Box<dyn AlignmentStream<RecordBuf>>;
    HashLookup::new(cfg, smallvec![s0, s1], [None, None], [None, None])
}

#[test]
fn test_hash_perfect_vs_imperfect_stream0_wins() -> Result<(), Error> {
    let s0 = vec![create_record(b"R1", "10M", &[], &[], "10", false)?];
    let s1 = vec![create_record(b"R1", "5M5S", &[], &[], "5", false)?];
    let mut h = make_lookup(s0, s1, config())?;
    h.process()?;
    assert_eq!(h.routing_counters[1], 1); // out:0
    assert_eq!(h.routing_counters[4], 1); // discard:1
    Ok(())
}

#[test]
fn test_hash_perfect_vs_imperfect_stream1_wins() -> Result<(), Error> {
    let s0 = vec![create_record(b"R1", "5M5S", &[], &[], "5", false)?];
    let s1 = vec![create_record(b"R1", "10M", &[], &[], "10", false)?];
    let mut h = make_lookup(s0, s1, config())?;
    h.process()?;
    assert_eq!(h.routing_counters[0], 1); // discard:0
    assert_eq!(h.routing_counters[5], 1); // out:1
    Ok(())
}

#[test]
fn test_hash_tie_is_ambiguous() -> Result<(), Error> {
    let s0 = vec![create_record(b"R1", "10M", &[], &[], "10", false)?];
    let s1 = vec![create_record(b"R1", "10M", &[], &[], "10", false)?];
    let mut h = make_lookup(s0, s1, config())?;
    h.process()?;
    assert_eq!(h.routing_counters[2], 1); // ambig:0
    assert_eq!(h.routing_counters[6], 1); // ambig:1
    Ok(())
}

#[test]
fn test_hash_interleaved_streams() -> Result<(), Error> {
    // R1 arrives in stream0 first, then R2; stream1 has R2 first then R1.
    let s0 = vec![
        create_record(b"R1", "10M", &[], &[], "10", false)?,
        create_record(b"R2", "5M5S", &[], &[], "5", false)?,
    ];
    let s1 = vec![
        create_record(b"R2", "10M", &[], &[], "10", false)?,
        create_record(b"R1", "5M5S", &[], &[], "5", false)?,
    ];
    let mut h = make_lookup(s0, s1, config())?;
    h.process()?;
    assert_eq!(h.routing_counters[1], 1); // out:0 (R1)
    assert_eq!(h.routing_counters[5], 1); // out:1 (R2)
    Ok(())
}

#[test]
fn test_hash_unmapped_vs_mapped() -> Result<(), Error> {
    let s0 = vec![create_record(b"R1", "", &[b'A'; 10], &[30; 10], "", false)?];
    let s1 = vec![create_record(b"R1", "10M", &[], &[], "10", false)?];
    let mut h = make_lookup(s0, s1, config())?;
    h.process()?;
    assert_eq!(h.routing_counters[0], 1); // discard:0
    assert_eq!(h.routing_counters[5], 1); // out:1
    Ok(())
}

#[test]
fn test_hash_suffix_stripping() -> Result<(), Error> {
    let cfg = Config {
        strip_read_suffix: StripReadSuffix::True,
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Config::default()
    };
    let s0 = vec![create_record(b"R1/1", "10M", &[], &[], "10", false)?];
    let s1 = vec![create_record(b"R1/2", "5M5S", &[], &[], "5", false)?];
    let mut h = make_lookup(s0, s1, cfg)?;
    h.process()?;
    assert_eq!(h.routing_counters[1], 1); // out:0
    assert_eq!(h.routing_counters[4], 1); // discard:1
    Ok(())
}

#[test]
fn test_hash_paired_end_both_mates_grouped() -> Result<(), Error> {
    let s0 = vec![
        create_record(b"R1", "10M", &[], &[], "10", false)?,
        create_record(b"R1", "10M", &[], &[], "10", true)?,
    ];
    let s1 = vec![
        create_record(b"R1", "5M5S", &[], &[], "5", false)?,
        create_record(b"R1", "5M5S", &[], &[], "5", true)?,
    ];
    let mut h = make_lookup(s0, s1, config())?;
    h.process()?;
    assert_eq!(h.routing_counters[1], 2); // out:0 (both mates)
    assert_eq!(h.routing_counters[4], 2); // discard:1 (both mates)
    Ok(())
}

#[test]
fn test_hash_multiple_fragments() -> Result<(), Error> {
    let s0 = vec![
        create_record(b"R1", "10M", &[], &[], "10", false)?,
        create_record(b"R2", "10M", &[], &[], "10", false)?,
        create_record(b"R3", "5M5S", &[], &[], "5", false)?,
    ];
    let s1 = vec![
        create_record(b"R1", "5M5S", &[], &[], "5", false)?,
        create_record(b"R2", "5M5S", &[], &[], "5", false)?,
        create_record(b"R3", "10M", &[], &[], "10", false)?,
    ];
    let mut h = make_lookup(s0, s1, config())?;
    h.process()?;
    // R1: stream0 perfect wins; R2: stream0 perfect wins; R3: stream1 perfect wins
    assert_eq!(h.routing_counters[1], 2); // out:0 (R1, R2)
    assert_eq!(h.routing_counters[5], 1); // out:1 (R3)
    Ok(())
}
