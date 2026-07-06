use crate::{
    alignment::FragmentState,
    aln_stream::AlignmentStream,
    config::Config,
    filter_algorithm::line_by_line::core::FragmentBuffer,
    tests::{create_record, MockStream},
    Error, LineByLine,
};
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{smallvec, SmallVec};

// ---------------------------------------------------------------------------
// Builder helpers
// ---------------------------------------------------------------------------

fn lbl_from(specs: &[(&str, Vec<RecordBuf>)]) -> LineByLine<RecordBuf> {
    let aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = specs
        .iter()
        .enumerate()
        .map(|(i, (_, recs))| {
            Box::new(MockStream::new(i, recs.clone())) as Box<dyn AlignmentStream<RecordBuf>>
        })
        .collect();
    LineByLine::new(Config::default(), aln).unwrap()
}

fn lbl_chimeric(specs: &[(&str, Vec<RecordBuf>)], pairs: &[[usize; 2]]) -> LineByLine<RecordBuf> {
    let mut cfg = Config::default();
    cfg.parsed_chimeric_pairs = pairs.to_vec();
    cfg.stream_labels = specs.iter().map(|(label, _)| label.to_string()).collect();
    let aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = specs
        .iter()
        .enumerate()
        .map(|(i, (_, recs))| {
            Box::new(MockStream::new(i, recs.clone())) as Box<dyn AlignmentStream<RecordBuf>>
        })
        .collect();
    LineByLine::new(cfg, aln).unwrap()
}

fn rec(name: &[u8], cigar: &str, md: &str) -> RecordBuf {
    create_record(name, cigar, &[], &[30u8; 20], md, false).unwrap()
}

fn unmap(name: &[u8]) -> RecordBuf {
    let seq = vec![b'A'; 10];
    let q = vec![30u8; 10];
    create_record(name, "", &seq, &q, "", false).unwrap()
}

/// Set flag bits and retain name.
fn with_flags(mut r: RecordBuf, bits: u16) -> RecordBuf {
    *r.flags_mut() = Flags::from_bits(bits).unwrap();
    r
}

// Counter access helpers (layout: stream * 4 + category, 0=disc 1=out 2=ambig 3=chimeric)
fn disc(lbl: &LineByLine<RecordBuf>, stream: usize) -> u64 {
    lbl.routing_counters[stream * 4]
}
fn out(lbl: &LineByLine<RecordBuf>, stream: usize) -> u64 {
    lbl.routing_counters[stream * 4 + 1]
}
fn ambig(lbl: &LineByLine<RecordBuf>, stream: usize) -> u64 {
    lbl.routing_counters[stream * 4 + 2]
}
fn chim(lbl: &LineByLine<RecordBuf>, stream: usize) -> u64 {
    lbl.routing_counters[stream * 4 + 3]
}

// ---------------------------------------------------------------------------
// 2-stream tournament — table-driven
// ---------------------------------------------------------------------------

struct TwoStreamCase {
    label: &'static str,
    s0: Vec<RecordBuf>,
    s1: Vec<RecordBuf>,
    out0: u64,
    out1: u64,
    disc0: u64,
    disc1: u64,
    ambig0: u64,
    ambig1: u64,
}

fn run_2stream(cases: &[TwoStreamCase]) {
    for c in cases {
        let mut lbl = lbl_from(&[("a", c.s0.clone()), ("b", c.s1.clone())]);
        lbl.process_sequential().unwrap();
        assert_eq!(out(&lbl, 0), c.out0, "[{}] out[0]", c.label);
        assert_eq!(out(&lbl, 1), c.out1, "[{}] out[1]", c.label);
        assert_eq!(disc(&lbl, 0), c.disc0, "[{}] disc[0]", c.label);
        assert_eq!(disc(&lbl, 1), c.disc1, "[{}] disc[1]", c.label);
        assert_eq!(ambig(&lbl, 0), c.ambig0, "[{}] ambig[0]", c.label);
        assert_eq!(ambig(&lbl, 1), c.ambig1, "[{}] ambig[1]", c.label);
    }
}

#[test]
fn two_stream_tournament() {
    run_2stream(&[
        // ── Tier 1: unmapped ─────────────────────────────────────────────
        TwoStreamCase {
            label: "both unmapped → both ambiguous",
            s0: vec![unmap(b"R1")], s1: vec![unmap(b"R1")],
            out0:0, out1:0, disc0:0, disc1:0, ambig0:1, ambig1:1,
        },
        TwoStreamCase {
            label: "s0 unmapped s1 mapped → s1 wins",
            s0: vec![unmap(b"R1")], s1: vec![rec(b"R1", "10M", "10")],
            out0:0, out1:1, disc0:1, disc1:0, ambig0:0, ambig1:0,
        },
        // ── Tier 2: perfect-match ─────────────────────────────────────────
        TwoStreamCase {
            label: "s0 perfect s1 imperfect → s0 wins",
            s0: vec![rec(b"R1", "10M", "10")],
            s1: vec![rec(b"R1", "10M", "5A4")],
            out0:1, out1:0, disc0:0, disc1:1, ambig0:0, ambig1:0,
        },
        TwoStreamCase {
            label: "both perfect → ambiguous",
            s0: vec![rec(b"R1", "10M", "10")], s1: vec![rec(b"R1", "10M", "10")],
            out0:0, out1:0, disc0:0, disc1:0, ambig0:1, ambig1:1,
        },
        TwoStreamCase {
            label: "s1 perfect s0 softclip → s1 wins",
            s0: vec![rec(b"R1", "5S5M", "5")],
            s1: vec![rec(b"R1", "10M",  "10")],
            out0:0, out1:1, disc0:1, disc1:0, ambig0:0, ambig1:0,
        },
        // ── Tier 2.5: match-count domination ─────────────────────────────
        TwoStreamCase {
            label: "s0 more matches in CIGAR/MD → s0 wins via Tier2.5",
            // s0: 8 matches, s1: 6 matches (both imperfect so Tier2 doesn't resolve)
            s0: vec![rec(b"R1", "10M", "8AA")],
            s1: vec![rec(b"R1", "10M", "6AAAA")],
            out0:1, out1:0, disc0:0, disc1:1, ambig0:0, ambig1:0,
        },
        // ── Tier 3: NW scoring breaks tie ────────────────────────────────
        TwoStreamCase {
            label: "equal match count, NW score from quality breaks tie — not testable with MockStream quality=30 constant",
            // With all q=30 and flat MD, both identical CIGARs → ambiguous
            s0: vec![rec(b"R1", "8M2S", "8")],
            s1: vec![rec(b"R1", "8M2S", "8")],
            out0:0, out1:0, disc0:0, disc1:0, ambig0:1, ambig1:1,
        },
        // ── Multiple fragments ────────────────────────────────────────────
        TwoStreamCase {
            label: "multiple fragments independent outcomes",
            s0: vec![rec(b"R1","10M","10"), rec(b"R2","5S5M","5"),  rec(b"R3","10M","10")],
            s1: vec![rec(b"R1","5S5M","5"), rec(b"R2","10M","10"),  rec(b"R3","10M","10")],
            out0:1, out1:1, disc0:1, disc1:1, ambig0:1, ambig1:1,
        },
    ]);
}

// ---------------------------------------------------------------------------
// 3-stream tournament
// ---------------------------------------------------------------------------

#[test]
fn three_stream_tournament() {
    struct Row {
        label: &'static str,
        s: [Vec<RecordBuf>; 3],
        out: [u64; 3],
        disc: [u64; 3],
        ambig: [u64; 3],
    }
    let cases = &[
        Row {
            label: "s0 perfect others imperfect → s0 wins",
            s: [
                vec![rec(b"R1", "10M", "10")],
                vec![rec(b"R1", "10M", "5A4")],
                vec![rec(b"R1", "10M", "4A4A0")],
            ],
            out: [1, 0, 0],
            disc: [0, 1, 1],
            ambig: [0, 0, 0],
        },
        Row {
            label: "all three perfect → all ambiguous",
            s: [
                vec![rec(b"R1", "10M", "10")],
                vec![rec(b"R1", "10M", "10")],
                vec![rec(b"R1", "10M", "10")],
            ],
            out: [0, 0, 0],
            disc: [0, 0, 0],
            ambig: [1, 1, 1],
        },
        Row {
            label: "s2 has most matches → s2 wins via Tier 2.5",
            s: [
                vec![rec(b"R1", "10M", "6AAAA")], // 6 matches
                vec![rec(b"R1", "10M", "8AA")],   // 8 matches
                vec![rec(b"R1", "10M", "9A0")],   // 9 matches
            ],
            out: [0, 0, 1],
            disc: [1, 1, 0],
            ambig: [0, 0, 0],
        },
        Row {
            label: "s0 unmapped, s1/s2 mapped, s1 perfect → s1 wins",
            s: [
                vec![unmap(b"R1")],
                vec![rec(b"R1", "10M", "10")],
                vec![rec(b"R1", "10M", "8AA")],
            ],
            out: [0, 1, 0],
            disc: [1, 0, 1],
            ambig: [0, 0, 0],
        },
    ];
    for c in cases {
        let mut lbl = lbl_from(&[
            ("a", c.s[0].clone()),
            ("b", c.s[1].clone()),
            ("c", c.s[2].clone()),
        ]);
        lbl.process_sequential().unwrap();
        for i in 0..3 {
            assert_eq!(out(&lbl, i), c.out[i], "[{}] out[{i}]", c.label);
            assert_eq!(disc(&lbl, i), c.disc[i], "[{}] disc[{i}]", c.label);
            assert_eq!(ambig(&lbl, i), c.ambig[i], "[{}] ambig[{i}]", c.label);
        }
    }
}

// ---------------------------------------------------------------------------
// Chimeric detection — mate-split and read-split
// ---------------------------------------------------------------------------

/// Build a paired-end record: flags encode first/last segment.
fn pe(name: &[u8], cigar: &str, md: &str, flags_bits: u16) -> RecordBuf {
    with_flags(rec(name, cigar, md), flags_bits)
}

#[test]
fn chimeric_mate_split() {
    // read1 → stream 0 (human), read2 → stream 1 (HPV)
    // flags: 0x41 = read1+paired, 0x81 = read2+paired
    let s0 = vec![pe(b"R1", "10M", "10", 0x41)]; // read1 in human
    let s1 = vec![pe(b"R1", "10M", "10", 0x81)]; // read2 in HPV

    let mut lbl = lbl_chimeric(&[("human", s0), ("hpv", s1)], &[[0, 1]]);
    lbl.process_sequential().unwrap();

    // Both streams' records written as winners with XC:Z tag.
    assert_eq!(chim(&lbl, 0), 1, "human chimeric count");
    assert_eq!(chim(&lbl, 1), 1, "hpv chimeric count");
    assert_eq!(out(&lbl, 0), 0, "human normal-out should be 0");
    assert_eq!(out(&lbl, 1), 0, "hpv normal-out should be 0");
}

#[test]
fn chimeric_read_split_complementary_clips() {
    // read1: human maps [0,25), HPV maps [25,50) — complementary 25S25M / 25M25S
    // read2: entirely in HPV (0x81)
    // Both streams' read1 primary — segment IDs *not* disjoint (both claim read1),
    // so mate-split detection doesn't fire. Read-split detection must fire.

    use noodles::core::Position;

    let q = vec![30u8; 50];

    // Stream 0 (human): read1 25M25S, read2 unmapped
    let mut r1_human = create_record(b"R1", "25M25S", &[b'A'; 50], &q, "25", false).unwrap();
    *r1_human.flags_mut() = Flags::from_bits(0x41).unwrap();
    *r1_human.alignment_start_mut() = Some(Position::new(1).unwrap());
    *r1_human.reference_sequence_id_mut() = Some(0);

    // Stream 1 (HPV): read1 25S25M, read2 fully mapped
    let mut r1_hpv = create_record(b"R1", "25S25M", &[b'A'; 50], &q, "25", false).unwrap();
    *r1_hpv.flags_mut() = Flags::from_bits(0x41).unwrap();
    *r1_hpv.alignment_start_mut() = Some(Position::new(1).unwrap());
    *r1_hpv.reference_sequence_id_mut() = Some(0);

    let r2_hpv = pe(b"R1", "50M", "50", 0x81);

    let mut lbl = lbl_chimeric(
        &[("human", vec![r1_human]), ("hpv", vec![r1_hpv, r2_hpv])],
        &[[0, 1]],
    );
    lbl.process_sequential().unwrap();

    assert_eq!(chim(&lbl, 0), 1, "human chimeric");
    assert_eq!(chim(&lbl, 1), 2, "hpv chimeric: read1 + read2");
}

#[test]
fn chimeric_false_positive_rejected_when_supp_better() {
    // Stream A has supplementary alignment for read1 with 0 mismatches (better than
    // stream B's primary for the same region). Should NOT be called chimeric.
    use noodles::core::Position;

    let q = vec![30u8; 50];

    // Stream A primary: read1 25M25S
    let mut r1_a = create_record(b"R1", "25M25S", &[b'A'; 50], &q, "25", false).unwrap();
    *r1_a.flags_mut() = Flags::from_bits(0x41).unwrap();
    *r1_a.alignment_start_mut() = Some(Position::new(1).unwrap());
    *r1_a.reference_sequence_id_mut() = Some(0);

    // Stream A supplementary: covers the 3' portion on A's reference, 0 mismatches
    let mut r1_a_supp = create_record(b"R1", "25S25M", &[b'A'; 50], &q, "25", false).unwrap();
    *r1_a_supp.flags_mut() = Flags::from_bits(0x41 | 0x800).unwrap(); // supplementary
    *r1_a_supp.alignment_start_mut() = Some(Position::new(100).unwrap());
    *r1_a_supp.reference_sequence_id_mut() = Some(0);

    // Stream B primary for same read: 25S25M, but MD has mismatches → worse than A's supp
    let mut r1_b = create_record(b"R1", "25S25M", &[b'A'; 50], &q, "15AAAAAAAAAA", false).unwrap();
    *r1_b.flags_mut() = Flags::from_bits(0x41).unwrap();
    *r1_b.alignment_start_mut() = Some(Position::new(1).unwrap());
    *r1_b.reference_sequence_id_mut() = Some(0);

    let mut lbl = lbl_chimeric(
        &[("a", vec![r1_a, r1_a_supp]), ("b", vec![r1_b])],
        &[[0, 1]],
    );
    lbl.process_sequential().unwrap();

    assert_eq!(
        chim(&lbl, 0),
        0,
        "should NOT be chimeric: A's supplementary is better"
    );
    assert_eq!(chim(&lbl, 1), 0);
}

#[test]
fn chimeric_three_stream_pair_01_mouse_competes_normally() {
    // human=0, hpv=1, mouse=2 — chimeric pair [0,1]
    // Fragment: read1→human, read2→hpv — chimeric event on [0,1]
    // mouse records should be discarded normally
    let r1_human = pe(b"R1", "10M", "10", 0x41);
    let r2_hpv = pe(b"R1", "10M", "10", 0x81);
    let r_mouse = rec(b"R1", "10M", "5A4"); // loses normally

    let mut lbl = lbl_chimeric(
        &[
            ("human", vec![r1_human]),
            ("hpv", vec![r2_hpv]),
            ("mouse", vec![r_mouse]),
        ],
        &[[0, 1]],
    );
    lbl.process_sequential().unwrap();

    assert_eq!(chim(&lbl, 0), 1, "human chimeric");
    assert_eq!(chim(&lbl, 1), 1, "hpv chimeric");
    assert_eq!(disc(&lbl, 2), 1, "mouse discarded (not in chimeric pair)");
    assert_eq!(out(&lbl, 2), 0);
}

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
fn test_branch_counters_and_skipping() -> Result<(), Error> {
    let mut config = Config {
        discard_unmapped: true,
        skip_secondary: true,
        ..Config::default()
    };

    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config.clone(), setup_mock_streams())?;

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
    assert_eq!(lbl.routing_counters[2], 2); // ambiguous:0: 2
    assert_eq!(lbl.routing_counters[0], 1); // discard:0:
                                            // ingest_record should skip secondary
    let mut best: FragmentBuffer<RecordBuf> = smallvec![];
    let finished = lbl.ingest_record(0, secondary, &mut best).unwrap();
    assert!(!finished);
    assert!(best.is_empty());

    config.discard_unmapped = false;
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config, setup_mock_streams())?;
    assert!(lbl.write_record(0, unmapped_fwd, None).is_ok());
    assert!(lbl.write_record(0, unmapped_rev, None).is_ok());
    assert_eq!(lbl.routing_counters[2], 2); // ambiguous:0: 2

    Ok(())
}

#[test]
fn test_handle_ordering_logic() -> Result<(), Error> {
    let lbl_setup: LineByLine<RecordBuf> =
        LineByLine::new(Config::default(), setup_mock_streams())?;
    // Direct testing of routing_counters incrementation via write_record:
    let mut lbl: LineByLine<RecordBuf> = lbl_setup;
    let rec = create_record(b"r1", "M10", &[], &[], "10", false)?;
    let rec2 = create_record(b"r2", "M10", &[], &[], "10", false)?;
    let rec3 = create_record(b"r3", "M10", &[], &[], "10", false)?;

    lbl.write_record(0, rec, Some(true))?;
    lbl.write_record(0, rec2, Some(false))?;
    lbl.write_record(0, rec3, None)?;

    assert_eq!(lbl.routing_counters[1], 1); // 1 + (0 << 1)
    assert_eq!(lbl.routing_counters[0], 1); // (0 << 1)
    assert_eq!(lbl.routing_counters[2], 1); // 16 + 0
    Ok(())
}

#[test]
fn test_fragment_finished_transitions() -> Result<(), Error> {
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(Config::default(), setup_mock_streams())?;
    let rec = create_record(b"R1", "M10", &[], &[], "10", false)?;
    let mut best: FragmentBuffer<RecordBuf> =
        smallvec![FragmentState::from_record(rec.clone(), 0)?];

    // Same QName: continues fragment
    let fin = lbl.ingest_record(0, rec, &mut best)?;
    assert!(!fin);
    assert_eq!(best[0].get_records().len(), 2);

    // Different QName: finishes fragment
    // Note: this will attempt to call aln[i].un_next(), requiring a mock AlnStream.
    Ok(())
}

#[test]
fn test_complex_fragment_grouping() -> Result<(), Error> {
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(Config::default(), setup_mock_streams())?;
    let mut best: FragmentBuffer<RecordBuf> = smallvec![];

    // paired-end style: same QNAME twice
    let r1_1 = create_record(b"READX", "10M", &[], &[], "10", false)?;
    let r1_2 = create_record(b"READX", "5M5S", &[], &[], "5", false)?;

    lbl.ingest_record(0, r1_1, &mut best)?;
    assert_eq!(best.len(), 1);
    assert_eq!(best[0].get_records().len(), 1);

    lbl.ingest_record(0, r1_2, &mut best)?;
    assert_eq!(best.len(), 1);
    assert_eq!(best[0].get_records().len(), 2);

    Ok(())
}

#[test]
fn test_line_by_line_full_flow() -> Result<(), Error> {
    let rec1 = create_record(b"R1", "10M", &[], &[], "10", false)?;
    let rec2 = create_record(b"R2", "10M", &[], &[], "10", false)?;

    // Mocking AlnStream behavior
    // Note: You may need to wrap MockStream in AlnStream enum/trait if required by your types
    // This targets ingest_record coverage
    let config = Config::default();
    let mut lbl: LineByLine<RecordBuf> = LineByLine::new(config, smallvec![])?;

    let mut best: FragmentBuffer<RecordBuf> = smallvec![];

    // 1. First read
    let fin = lbl.ingest_record(0, rec1.clone(), &mut best)?;
    assert!(!fin);
    assert_eq!(best.len(), 1);

    // 2. QName change triggers finish
    let fin = lbl.ingest_record(0, rec2, &mut best)?;
    // This will attempt to call aln[0].un_next(), ensure your mock handles this.
    assert!(fin);
    Ok(())
}

#[test]
fn test_observed_pe_scoring1() -> Result<(), Error> {
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
    assert_eq!(lbl.routing_counters[4], 2); // discard:1: both reads
    assert_eq!(lbl.routing_counters[5], 0); // out:1:
    assert_eq!(lbl.routing_counters[6], 0); // ambiguous:1:
    assert_eq!(lbl.routing_counters[3], 0); // unmapped:1:
    assert_eq!(lbl.routing_counters[0], 0); // discard:0:
    assert_eq!(lbl.routing_counters[1], 2); // out:0: both reads
    assert_eq!(lbl.routing_counters[2], 0); // ambiguous:0:
    assert_eq!(lbl.routing_counters[3], 0); // unmapped:0:
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
fn test_ambiguous_log_threshold_conversion() -> Result<(), Error> {
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
