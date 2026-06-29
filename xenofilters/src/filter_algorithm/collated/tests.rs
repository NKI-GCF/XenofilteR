use crate::aln_stream::tests::MockStream;
use crate::config::{Config, StripReadSuffix};
use crate::filter_algorithm::collated::CollatedMatcher;
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

/// Build a CollatedMatcher from two lists of records.
fn make_matcher(
    stream0: Vec<RecordBuf>,
    stream1: Vec<RecordBuf>,
    cfg: Config,
) -> Result<CollatedMatcher<RecordBuf>, Error> {
    use crate::aln_stream::AlignmentStream;
    let s0 = Box::new(MockStream::new(0, stream0)) as Box<dyn AlignmentStream<RecordBuf>>;
    let s1 = Box::new(MockStream::new(1, stream1)) as Box<dyn AlignmentStream<RecordBuf>>;
    CollatedMatcher::new(cfg, smallvec![s0, s1], [None, None], [None, None])
}

#[test]
fn test_collated_perfect_vs_imperfect_stream0_wins() -> Result<(), Error> {
    let s0 = vec![create_record(b"R1", "10M", &[], &[], "10", false)?];
    let s1 = vec![create_record(b"R1", "5M5S", &[], &[], "5", false)?];
    let mut m = make_matcher(s0, s1, config())?;
    m.process()?;
    assert_eq!(m.routing_counters[1], 1); // out:0
    assert_eq!(m.routing_counters[4], 1); // discard:1
    Ok(())
}

#[test]
fn test_collated_perfect_vs_imperfect_stream1_wins() -> Result<(), Error> {
    let s0 = vec![create_record(b"R1", "5M5S", &[], &[], "5", false)?];
    let s1 = vec![create_record(b"R1", "10M", &[], &[], "10", false)?];
    let mut m = make_matcher(s0, s1, config())?;
    m.process()?;
    assert_eq!(m.routing_counters[0], 1); // discard:0
    assert_eq!(m.routing_counters[5], 1); // out:1
    Ok(())
}

#[test]
fn test_collated_tie_is_ambiguous() -> Result<(), Error> {
    let s0 = vec![create_record(b"R1", "10M", &[], &[], "10", false)?];
    let s1 = vec![create_record(b"R1", "10M", &[], &[], "10", false)?];
    let mut m = make_matcher(s0, s1, config())?;
    m.process()?;
    assert_eq!(m.routing_counters[2], 1); // ambig:0
    assert_eq!(m.routing_counters[6], 1); // ambig:1
    Ok(())
}

#[test]
fn test_collated_streams_in_different_order() -> Result<(), Error> {
    // Stream 0: R1 then R2. Stream 1: R2 then R1.
    let s0 = vec![
        create_record(b"R1", "10M", &[], &[], "10", false)?,
        create_record(b"R2", "5M5S", &[], &[], "5", false)?,
    ];
    let s1 = vec![
        create_record(b"R2", "10M", &[], &[], "10", false)?,
        create_record(b"R1", "5M5S", &[], &[], "5", false)?,
    ];
    let mut m = make_matcher(s0, s1, config())?;
    m.process()?;
    // R1: stream0 perfect wins
    assert_eq!(m.routing_counters[1], 1); // out:0 (R1)
                                         // R2: stream1 perfect wins
    assert_eq!(m.routing_counters[5], 1); // out:1 (R2)
    Ok(())
}

#[test]
fn test_collated_paired_end_same_name_grouped() -> Result<(), Error> {
    // Both reads of a pair share a name and appear consecutively.
    let s0 = vec![
        create_record(b"R1", "10M", &[], &[], "10", false)?,
        create_record(b"R1", "10M", &[], &[], "10", true)?, // mate
    ];
    let s1 = vec![
        create_record(b"R1", "5M5S", &[], &[], "5", false)?,
        create_record(b"R1", "5M5S", &[], &[], "5", true)?,
    ];
    let mut m = make_matcher(s0, s1, config())?;
    m.process()?;
    // stream0 perfect on both reads vs imperfect stream1
    assert_eq!(m.routing_counters[1], 2); // out:0 (both mates)
    assert_eq!(m.routing_counters[4], 2); // discard:1 (both mates)
    Ok(())
}

#[test]
fn test_collated_unmapped_vs_mapped() -> Result<(), Error> {
    let s0 = vec![create_record(b"R1", "", &[b'A'; 10], &[30; 10], "", false)?];
    let s1 = vec![create_record(b"R1", "10M", &[], &[], "10", false)?];
    let mut m = make_matcher(s0, s1, config())?;
    m.process()?;
    // unmapped < mapped: stream1 wins
    assert_eq!(m.routing_counters[0], 1); // discard:0
    assert_eq!(m.routing_counters[5], 1); // out:1
    Ok(())
}

#[test]
fn test_collated_suffix_stripping() -> Result<(), Error> {
    let cfg = Config {
        strip_read_suffix: StripReadSuffix::True,
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Config::default()
    };
    // /1 and /2 should be treated as the same fragment.
    let s0 = vec![create_record(b"R1/1", "10M", &[], &[], "10", false)?];
    let s1 = vec![create_record(b"R1/2", "5M5S", &[], &[], "5", false)?];
    let mut m = make_matcher(s0, s1, cfg)?;
    m.process()?;
    assert_eq!(m.routing_counters[1], 1); // out:0
    assert_eq!(m.routing_counters[4], 1); // discard:1
    Ok(())
}

#[test]
fn test_collated_unmatched_is_emitted_as_best() -> Result<(), Error> {
    // R2 only in stream0; stream1 has nothing matching.
    let s0 = vec![
        create_record(b"R1", "10M", &[], &[], "10", false)?,
        create_record(b"R2", "10M", &[], &[], "10", false)?,
    ];
    let s1 = vec![
        create_record(b"R1", "5M5S", &[], &[], "5", false)?,
        // R2 absent
    ];
    let mut m = make_matcher(s0, s1, config())?;
    m.process()?;
    // R1: stream0 wins normally
    assert_eq!(m.routing_counters[1], 2); // out:0
    assert_eq!(m.routing_counters[0], 0); // discard:0
    assert_eq!(m.routing_counters[4], 1); // discard:1 R1
    assert_eq!(m.routing_counters[5], 0); // discard:1 R1
    Ok(())
}
