use crate::alignment::FragmentState;
use crate::aln_stream::AlignmentStream;
use crate::config::Config;
use crate::filter_algorithm::line_by_line::{core::AlnBuffer, ordering::Decision};
use crate::tests::{create_record, MockStream};
use crate::LineByLine;
use anyhow::Result;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{smallvec, SmallVec};
use std::cmp::Ordering;

pub(crate) fn setup_mock_streams() -> SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> {
    let stream1 = MockStream::new(
        0,
        vec![
            create_record(b"R0", "10M", &[], &[], "10", false).unwrap(), // perfect => out
            create_record(b"R1", "10M", &[], &[], "10", false).unwrap(), // perfect => out
            create_record(b"R2", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => discarded
            create_record(b"R3", "10M", &[], &[], "10", false).unwrap(), // perfect => out
            create_record(b"R4", "*", &[], &[], "10", false).unwrap(),   // unmapped => discarded
            create_record(b"R5", "10M", &[], &[], "10", false).unwrap(), // perfect => ambiguous
            create_record(b"R6", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => ambiguous
            create_record(b"R7", "*", &[], &[], "10", false).unwrap(),   // unmapped => ambiguous
            create_record(b"R8", "6M4S", &[], &[], "6", false).unwrap(), // less clipped => out
        ],
    );
    let stream2 = MockStream::new(
        1,
        vec![
            create_record(b"R0", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => discarded
            create_record(b"R1", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => discarded
            create_record(b"R2", "10M", &[], &[], "10", false).unwrap(), // perfect => out
            create_record(b"R3", "*", &[], &[], "10", false).unwrap(),   // unmapped => discarded
            create_record(b"R4", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => out
            create_record(b"R5", "10M", &[], &[], "10", false).unwrap(), // perfect => ambiguous
            create_record(b"R6", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => ambiguous
            create_record(b"R7", "*", &[], &[], "9", false).unwrap(),    // unmapped => ambiguous
            create_record(b"R8", "5M5S", &[], &[], "5", false).unwrap(), // less match => discarded
        ],
    );
    smallvec![
        Box::new(stream1) as Box<dyn AlignmentStream<RecordBuf>>,
        Box::new(stream2) as Box<dyn AlignmentStream<RecordBuf>>
    ]
}

fn setup_mock_streams_r4() -> SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> {
    let stream1 = MockStream::new(
        0,
        vec![
            create_record(b"R4", "*", &[], &[], "10", false).unwrap(), // unmapped => discarded
        ],
    );
    let stream2 = MockStream::new(
        1,
        vec![
            create_record(b"R4", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => out
        ],
    );
    smallvec![
        Box::new(stream1) as Box<dyn AlignmentStream<RecordBuf>>,
        Box::new(stream2) as Box<dyn AlignmentStream<RecordBuf>>
    ]
}

fn setup_mock_streams_observed_examples() -> SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> {
    // human
    // FRAGMENT1/1        100M    MD:Z:86T13
    // FRAGMENT1/2        100M    MD:Z:100
    // mouse
    // FRAGMENT1/1        23M4D37M40S     MD:Z:13A9^TAAT10C19C6
    // FRAGMENT1/2        100M    MD:Z:20T8G20T49
    let stream1 = MockStream::new(
        0,
        vec![
            create_record(b"FRAGMENT1/1", "100M", &[], &[], "86T13", false).unwrap(),
            create_record(b"FRAGMENT1/2", "100M", &[], &[], "100", true).unwrap(),
        ],
    );
    let stream2 = MockStream::new(
        1,
        vec![
            create_record(
                b"FRAGMENT1/1",
                "23M4D37M40S",
                &[],
                &[],
                "13A9^TAAT10C19C6",
                true,
            )
            .unwrap(),
            create_record(b"FRAGMENT1/2", "100M", &[], &[], "20T8G20T49", false).unwrap(),
        ],
    );
    smallvec![
        Box::new(stream1) as Box<dyn AlignmentStream<RecordBuf>>,
        Box::new(stream2) as Box<dyn AlignmentStream<RecordBuf>>
    ]
}

#[test]
fn test_branch_counters_and_skipping() -> Result<()> {
    let mut config = Config {
        discard_unmapped: true,
        skip_secondary: true,
        ..Config::default()
    };

    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config.clone(), smallvec![])?;

    let mut unmapped_fwd = create_record(b"u", "*", &[], &[], "10", false)?;
    *unmapped_fwd.flags_mut() = Flags::from_bits(0x45).unwrap(); // unmapped, paired, first in

    let mut unmapped_rev = unmapped_fwd.clone();
    *unmapped_rev.flags_mut() = Flags::from_bits(0x55).unwrap(); // reverse

    let mut secondary = create_record(b"s", "*", &[], &[], "10", false)?;
    *secondary.flags_mut() = Flags::from_bits(0x155).unwrap(); // secondary

    let mut unmapped_single = create_record(b"u2", "*", &[], &[], "10", false)?;
    *unmapped_single.flags_mut() = Flags::from_bits(0x4).unwrap(); // unmapped, single-end

    // Should return early (skipped)
    assert!(lbl.write_record(0, unmapped_fwd.clone(), None).is_ok());
    assert!(lbl.write_record(0, unmapped_rev.clone(), None).is_ok());
    assert!(lbl.write_record(0, unmapped_single, Some(false)).is_ok());
    lbl.print_counters(0);
    assert_eq!(lbl.branch_counters[16], 2); // ambiguous:0: 2
    assert_eq!(lbl.branch_counters[0], 1); // discard:0:
                                           // handle_record_is_fragment_finished should skip secondary
    let mut best: AlnBuffer<RecordBuf> = smallvec![];
    let finished = lbl
        .handle_record_is_fragment_finished(0, secondary, &mut best)
        .unwrap();
    assert!(!finished);
    assert!(best.is_empty());

    config.discard_unmapped = false;
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config, smallvec![])?;
    assert!(lbl.write_record(0, unmapped_fwd, None).is_ok());
    assert!(lbl.write_record(0, unmapped_rev, None).is_ok());
    lbl.print_counters(0);
    assert_eq!(lbl.branch_counters[16], 2); // ambiguous:0: 2

    Ok(())
}

#[test]
fn test_process_multi_stream_sync_r4() -> Result<()> {
    let config = Config {
        discard_unmapped: true,
        ..Config::default()
    };
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config, setup_mock_streams_r4())?;

    // R4 -> stream 1 (mismatch vs unmapped)

    lbl.process_sequential()?;

    assert_eq!(lbl.branch_counters[2], 0); // discard:1:
    assert_eq!(lbl.branch_counters[3], 1); // out:1: R4
    assert_eq!(lbl.branch_counters[17], 0); // ambiguous:1:
    assert_eq!(lbl.branch_counters[25], 0); // unmapped:1:
    assert_eq!(lbl.branch_counters[0], 1); // discard:0: R4
    assert_eq!(lbl.branch_counters[1], 0); // out:0:
    assert_eq!(lbl.branch_counters[16], 0); // ambiguous:0:
    assert_eq!(lbl.branch_counters[24], 0); // unmapped:0:
    Ok(())
}

#[test]
fn test_process_multi_stream_sync() -> Result<()> {
    let config = Config {
        discard_unmapped: true,
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Config::default()
    };

    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config, setup_mock_streams())?;

    // Streams now contain R1..R7
    // Expected winners:
    // R1 -> stream 0 (perfect vs mismatch)
    // R2 -> stream 1 (perfect vs mismatch)
    // R3 -> stream 0 (perfect vs unmapped)
    // R4 -> stream 1 (mismatch vs unmapped)
    // R5 -> tie perfect/perfect -> ambiguity
    // R6 -> tie mismatch/mismatch -> ambiguity
    // R7 -> unmapped/unmapped -> ambiguity
    // R8 -> stream 0 (mismatch vs more mismatches, but stream 1 discarded)

    lbl.process_sequential()?;
    // this is the order of printing, first aln 1 then aln 0
    assert_eq!(lbl.branch_counters[2], 4); // discard:1: R0, R1, R3, R8
    assert_eq!(lbl.branch_counters[3], 2); // out:1: R2, R4
    assert_eq!(lbl.branch_counters[17], 2); // ambiguous:1: R5, R6
    assert_eq!(lbl.branch_counters[25], 1); // unmapped:1: R7

    assert_eq!(lbl.branch_counters[0], 2); // discard:0: R2, R4
    assert_eq!(lbl.branch_counters[1], 4); // out:0: R0, R1, R3, R8
    assert_eq!(lbl.branch_counters[16], 2); // ambiguous:0: R5, R6
    assert_eq!(lbl.branch_counters[24], 1); // unmapped:0: R7

    Ok(())
}

#[test]
fn test_handle_ordering_logic() -> Result<()> {
    let lbl_setup: LineByLine<RecordBuf> = LineByLine::new(Config::default(), smallvec![])?;
    // Test logic in handle_ordering requires mocked AlnStream for write_record calls.
    // Direct testing of branch_counters incrementation via write_record:
    let mut lbl: LineByLine<RecordBuf> = lbl_setup;
    let rec = create_record(b"r1", "M10", &[], &[], "10", false)?;
    let rec2 = create_record(b"r2", "M10", &[], &[], "10", false)?;
    let rec3 = create_record(b"r3", "M10", &[], &[], "10", false)?;

    lbl.write_record(0, rec, Some(true))?;
    lbl.write_record(0, rec2, Some(false))?;
    lbl.write_record(0, rec3, None)?;

    assert_eq!(lbl.branch_counters[1], 1); // 1 + (0 << 1)
    assert_eq!(lbl.branch_counters[0], 1); // (0 << 1)
    assert_eq!(lbl.branch_counters[16], 1); // 16 + 0
    Ok(())
}

#[test]
fn test_fragment_finished_transitions() -> Result<()> {
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(Config::default(), smallvec![])?;
    let rec = create_record(b"R1", "M10", &[], &[], "10", false)?;
    let mut best: AlnBuffer<RecordBuf> = smallvec![FragmentState::from_record(rec.clone(), 0)?];

    // Same QName: continues fragment
    let fin = lbl.handle_record_is_fragment_finished(0, rec, &mut best)?;
    assert!(!fin);
    assert_eq!(best[0].get_records().len(), 2);

    // Different QName: finishes fragment
    // Note: this will attempt to call aln[i].un_next(), requiring a mock AlnStream.
    Ok(())
}

#[test]
fn test_handle_ordering_drain_logic() -> Result<()> {
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(Config::default(), setup_mock_streams())?;
    let mut best: AlnBuffer<RecordBuf> = smallvec![
        FragmentState::from_record(create_record(b"R1", "10M", &[], &[], "10", false)?, 0)?,
        FragmentState::from_record(create_record(b"R1", "5M5S", &[], &[], "5", false)?, 1)?,
    ];
    let mut ord = best[0].partial_cmp(&best[1]);
    assert_eq!(ord, None);

    let _ = best[0].cmp_perfect(&best[1], &mut ord)?;
    assert_eq!(ord, Some(Ordering::Greater));

    // stream 0 better than stream 1
    lbl.handle_ordering(&mut best, ord.unwrap())?;
    lbl.print_counters(1);

    assert_eq!(best.len(), 1);
    assert_eq!(best[0].get_nr(), 0);
    assert_eq!(lbl.branch_counters[1 + 1], 1); // stream 1 discarded
    eprintln!("Writing best fragment from stream {}", best[0].get_nr());
    lbl.handle_best(&mut best, None)?;
    lbl.print_counters(0);
    assert_eq!(lbl.branch_counters[1], 1); // stream 0 assigned

    Ok(())
}

#[test]
fn test_complex_fragment_grouping() -> Result<()> {
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(Config::default(), setup_mock_streams())?;
    let mut best: AlnBuffer<RecordBuf> = smallvec![];

    // paired-end style: same QNAME twice
    let r1_1 = create_record(b"READX", "10M", &[], &[], "10", false)?;
    let r1_2 = create_record(b"READX", "5M5S", &[], &[], "5", false)?;

    lbl.handle_record_is_fragment_finished(0, r1_1, &mut best)?;
    assert_eq!(best.len(), 1);
    assert_eq!(best[0].get_records().len(), 1);

    lbl.handle_record_is_fragment_finished(0, r1_2, &mut best)?;
    assert_eq!(best.len(), 1);
    assert_eq!(best[0].get_records().len(), 2);

    Ok(())
}

#[test]
fn test_line_by_line_full_flow() -> Result<()> {
    let rec1 = create_record(b"R1", "10M", &[], &[], "10", false)?;
    let rec2 = create_record(b"R2", "10M", &[], &[], "10", false)?;

    // Mocking AlnStream behavior
    // Note: You may need to wrap MockStream in AlnStream enum/trait if required by your types
    // This targets handle_record_is_fragment_finished coverage
    let config = Config::default();
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config, smallvec![])?;

    let mut best: AlnBuffer<RecordBuf> = smallvec![];

    // 1. First read
    let fin = lbl.handle_record_is_fragment_finished(0, rec1.clone(), &mut best)?;
    assert!(!fin);
    assert_eq!(best.len(), 1);

    // 2. QName change triggers finish
    let fin = lbl.handle_record_is_fragment_finished(0, rec2, &mut best)?;
    // This will attempt to call aln[0].un_next(), ensure your mock handles this.
    assert!(fin);
    Ok(())
}

#[test]
fn test_scoring_path_coverage() -> Result<()> {
    let config = Config::default();
    // Mock stream needs to exist to avoid indexing panics
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config, smallvec![])?;

    // Populate records with valid CIGAR/MD data to avoid null-pointer panics
    let mut best: AlnBuffer<RecordBuf> = smallvec![
        FragmentState::from_record(create_record(b"R1", "10M", &[], &[], "10", false)?, 0)?,
        FragmentState::from_record(create_record(b"R1", "10M", &[], &[], "10", false)?, 1)?,
    ];

    // Pass explicit Equal ordering to avoid scoring path (which needs non-empty aln)
    let result = lbl.handle_ordering(&mut best, Ordering::Equal);
    assert!(result.is_ok());
    Ok(())
}

#[test]
fn test_observed_pe_scoring1() -> Result<()> {
    let config = Config {
        discard_unmapped: true,
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Config::default()
    };

    let mut lbl: LineByLine<RecordBuf> =
        LineByLine::new(config, setup_mock_streams_observed_examples())?;

    lbl.process_sequential()?;
    assert_eq!(lbl.branch_counters[2], 2); // discard:1: both reads
    assert_eq!(lbl.branch_counters[3], 0); // out:1:
    assert_eq!(lbl.branch_counters[17], 0); // ambiguous:1:
    assert_eq!(lbl.branch_counters[25], 0); // unmapped:1:
    assert_eq!(lbl.branch_counters[0], 0); // discard:0:
    assert_eq!(lbl.branch_counters[1], 2); // out:0: both reads
    assert_eq!(lbl.branch_counters[16], 0); // ambiguous:0:
    assert_eq!(lbl.branch_counters[24], 0); // unmapped:0:
    Ok(())
}

#[cfg(test)]
impl LineByLine<RecordBuf> {
    /// Only for tests — lets us assert the converted log threshold.
    fn test_ambiguous_log_threshold(&self) -> f64 {
        self.ambiguous_log_threshold
    }
}

#[test]
fn test_ambiguous_log_threshold_conversion() -> Result<()> {
    let mut config = Config {
        ambiguous_threshold: 0,
        ..Config::default()
    };
    let aln = setup_mock_streams(); // any valid stream works for new()
    let aln_clone1 = setup_mock_streams(); // any valid stream works for new()
    let aln_clone2 = setup_mock_streams(); // any valid stream works for new()
    let aln_clone3 = setup_mock_streams(); // any valid stream works for new()

    // threshold = 0 → exactly 0.0 (or EPSILON if you changed it)
    let lbl: LineByLine<RecordBuf> = LineByLine::new(config.clone(), aln_clone1)?;
    assert_eq!(lbl.test_ambiguous_log_threshold(), 0.0);

    // standard phred values → correct natural-log ratio
    config.ambiguous_threshold = 10;
    let lbl: LineByLine<RecordBuf> = LineByLine::new(config.clone(), aln_clone2)?;
    let ln_10 = std::f64::consts::LN_10;
    assert!((lbl.test_ambiguous_log_threshold() - ln_10).abs() < 1e-9);

    config.ambiguous_threshold = 20;
    let lbl: LineByLine<RecordBuf> = LineByLine::new(config.clone(), aln_clone3)?;
    assert!((lbl.test_ambiguous_log_threshold() - ln_10 * 2.0).abs() < 1e-9);

    config.ambiguous_threshold = 3;
    let lbl: LineByLine<RecordBuf> = LineByLine::new(config, aln)?;
    assert!((lbl.test_ambiguous_log_threshold() - ln_10 * 3.0 / 10.0).abs() < 1e-9);

    Ok(())
}

#[test]
fn test_handle_ordering_quick_paths_respect_decision_tag() -> Result<()> {
    let config = Config {
        add_decision_tag: true,
        ..Config::default()
    };
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config, setup_mock_streams())?;

    let mut best: AlnBuffer<RecordBuf> = smallvec![
        FragmentState::from_record(create_record(b"R1", "10M", &[], &[], "10", false)?, 0)?,
        FragmentState::from_record(create_record(b"R1", "5M5S", &[], &[], "5", false)?, 1)?,
    ];

    let dec = lbl.handle_ordering(&mut best, Ordering::Greater)?;
    assert!(matches!(dec, Some(Decision::First)));

    let dec = lbl.handle_ordering(&mut best, Ordering::Less)?;
    assert!(matches!(dec, Some(Decision::Last)));

    let dec = lbl.handle_ordering(&mut best, Ordering::Equal)?;
    assert!(matches!(dec, Some(Decision::Ambiguous)));

    Ok(())
}

#[test]
fn test_handle_ordering_ambiguous_when_below_threshold_and_negative_delta() -> Result<()> {
    let config = Config {
        add_decision_tag: true,
        ambiguous_threshold: 30, // higher than your observed delta → forces ambiguous
        ..Config::default()
    };
    let aln = setup_mock_streams_observed_examples();
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config, aln)?;

    let mut best: AlnBuffer<RecordBuf> = smallvec![
        FragmentState::from_record(create_record(b"R1", "10M", &[], &[], "10", false)?, 0)?,
        FragmentState::from_record(create_record(b"R1", "5M5S", &[], &[], "5", false)?, 1)?,
    ];
    let mut ord = None;
    let (mcfs1, mcfs2) = best[0].cmp_perfect(&best[1], &mut ord)?;
    let (delta, s1_vd, s2_vd) = lbl.score_delta(&best, mcfs1, mcfs2)?;
    let decision = lbl.apply_delta(&mut best, delta, s1_vd, s2_vd)?;
    assert!(matches!(decision, Some(Decision::Ambiguous)));

    // Also test negative-delta case (we flip sign internally)
    // (reuse the same best – just swap first/last scores by changing order or use a second mock)
    // The code already does `delta = -delta` and then uses the positive value for phred.
    let best_flipped: AlnBuffer<RecordBuf> = smallvec![
        FragmentState::from_record(create_record(b"R1", "5M5S", &[], &[], "5", false)?, 1)?,
        FragmentState::from_record(create_record(b"R1", "10M", &[], &[], "10", false)?, 0)?,
    ];
    let (mcfs1, mcfs2) = best_flipped[0].cmp_perfect(&best_flipped[1], &mut ord)?;
    let (delta, s1_vd, s2_vd) = lbl.score_delta(&best, mcfs1, mcfs2)?;
    let decision_flipped = lbl.apply_delta(&mut best, delta, s1_vd, s2_vd)?;
    assert!(matches!(decision_flipped, Some(Decision::Ambiguous)));

    Ok(())
}

#[test]
fn test_handle_ordering_when_decision_tag_is_disabled() -> Result<()> {
    let config = Config {
        add_decision_tag: false,
        ..Config::default()
    };
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config, setup_mock_streams())?;

    let mut best: AlnBuffer<RecordBuf> = smallvec![
        FragmentState::from_record(create_record(b"R1", "10M", &[], &[], "10", false)?, 0)?,
        FragmentState::from_record(create_record(b"R1", "5M5S", &[], &[], "5", false)?, 1)?,
    ];
    let mut ord = None;
    let (mcfs1, mcfs2) = best[0].cmp_perfect(&best[1], &mut ord)?;
    let (delta, s1_vd, s2_vd) = lbl.score_delta(&best, mcfs1, mcfs2)?;
    let decision = lbl.apply_delta(&mut best, delta, s1_vd, s2_vd)?;
    assert!(decision.is_none()); // no Decision object when tag is off

    Ok(())
}
