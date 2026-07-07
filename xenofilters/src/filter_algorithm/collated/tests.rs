use crate::{
    aln_stream::tests::MockStream,
    aln_stream::AlignmentStream,
    config::{Config, StripReadSuffix},
    filter_algorithm::collated::CollatedMatcher,
    tests::create_record,
};
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::smallvec;

fn cfg() -> Config {
    Config {
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        strip_read_suffix: StripReadSuffix::False,
        ..Config::default()
    }
}

fn make(s0: Vec<RecordBuf>, s1: Vec<RecordBuf>, cfg: &Config) -> CollatedMatcher<RecordBuf> {
    let a0 = Box::new(MockStream::new(0, s0)) as Box<dyn AlignmentStream<RecordBuf>>;
    let a1 = Box::new(MockStream::new(1, s1)) as Box<dyn AlignmentStream<RecordBuf>>;
    CollatedMatcher::new(cfg, smallvec![a0, a1], [None, None], [None, None]).unwrap()
}

fn r(name: &[u8], cigar: &str, md: &str) -> RecordBuf {
    create_record(name, cigar, &[], &[30u8; 20], md, false).unwrap()
}
fn u(name: &[u8]) -> RecordBuf {
    create_record(name, "", &vec![b'A'; 10], &[30u8; 10], "", false).unwrap()
}

struct Row {
    label: &'static str,
    s0: Vec<RecordBuf>,
    s1: Vec<RecordBuf>,
    rc: u64,
}

#[test]
fn collated_table() {
    let cases: &[Row] = &[
        // -- Tier 1: unmapped ----------------------------------------------
        Row {
            label: "both unmapped → ambiguous",
            s0: vec![u(b"R1")],
            s1: vec![u(b"R1")],
            // counter per byte in order: [chimeric1, ambiguous1, out1, discarded1, chimeric0,
            // ambiguous0, out0, discarded0]
            rc: 0x01000100,
        },
        Row {
            label: "s0 unmapped s1 mapped → s1 wins",
            s0: vec![u(b"R1")],
            s1: vec![r(b"R1", "10M", "10")],
            rc: 0x00100001,
        },
        // -- Tier 2: perfect -----------------------------------------------
        Row {
            label: "s0 perfect s1 imperfect → s0 wins",
            s0: vec![r(b"R1", "10M", "10")],
            s1: vec![r(b"R1", "5M5S", "5")],
            rc: 0x00010010,
        },
        Row {
            label: "both perfect → ambiguous",
            s0: vec![r(b"R1", "10M", "10")],
            s1: vec![r(b"R1", "10M", "10")],
            rc: 0x01000100,
        },
        // -- Tier 2.5: match-count domination -----------------------------
        Row {
            label: "s1 more matches (9 vs 7) → s1 wins via pre-assess",
            s0: vec![r(b"R1", "10M", "7AAA")], // 7 matches
            s1: vec![r(b"R1", "10M", "8A1")],  // 9 matches
            rc: 0x00100001,
        },
        // -- Out-of-order fragments ----------------------------------------
        Row {
            label: "streams deliver R1/R2 in opposite order — collated handles",
            s0: vec![r(b"R1", "10M", "10"), r(b"R2", "5M5S", "5")],
            s1: vec![r(b"R2", "10M", "10"), r(b"R1", "5M5S", "5")],
            rc: 0x00110011,
        },
        // -- Unmatched fragments (no counterpart in other stream) ----------
        // XXX: this case should not happen in practice, but we should handle it gracefully.
        Row {
            label: "R2 missing in s1 — s0's R2 emitted as winner",
            s0: vec![r(b"R1", "10M", "10"), r(b"R2", "10M", "10")],
            s1: vec![r(b"R1", "5M5S", "5")],
            rc: 0x00010020,
        },
        // -- Paired-end ----------------------------------------------------
        Row {
            label: "paired-end perfect both streams → ambiguous for both mates",
            s0: vec![r(b"R1", "10M", "10"), r(b"R1", "10M", "10")],
            s1: vec![r(b"R1", "10M", "10"), r(b"R1", "10M", "10")],
            rc: 0x02000200,
        },
    ];
    let mut misses = vec![];
    let routing_names = ["chimeric", "ambiguous", "out", "discarded"];
    let config = cfg();
    for c in cases {
        eprintln!("Running test case: {}", c.label);
        let mut m = make(c.s0.clone(), c.s1.clone(), &config);
        m.process(&config).unwrap();
        let rc = &m.routing_counters;
        for i in 0..8 {
            let name = routing_names[i % 4];
            let stream = i / 4;
            let expected = c.rc >> (i * 4) & 0xf;
            if rc[i] != expected {
                misses.push(format!(
                    "[{}] aln{stream}_{name}[{i}] expected {expected}, got {}",
                    c.label, rc[i]
                ));
            }
        }
    }
    if !misses.is_empty() {
        panic!("{} test cases failed:\n{}", misses.len(), misses.join("\n"));
    }
}

#[test]
fn collated_suffix_stripping() {
    let config = Config {
        strip_read_suffix: StripReadSuffix::True,
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Config::default()
    };
    let s0 = vec![create_record(b"R1/1", "10M", &[], &[30u8; 10], "10", false).unwrap()];
    let s1 = vec![create_record(b"R1/2", "5M5S", &[], &[30u8; 10], "5", false).unwrap()];
    let mut m = make(s0, s1, &config);
    m.process(&config).unwrap();
    assert_eq!(m.routing_counters[1], 1, "s0 should win");
    assert_eq!(m.routing_counters[4], 1, "s1 discarded");
}
