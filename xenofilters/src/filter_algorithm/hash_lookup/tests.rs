use crate::{
    aln_stream::tests::MockStream,
    aln_stream::AlignmentStream,
    config::{Config, StripReadSuffix},
    filter_algorithm::hash_lookup::HashLookup,
    tests::create_record,
    Error,
};
use noodles::sam::alignment::record::Flags;
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

fn make(s0: Vec<RecordBuf>, s1: Vec<RecordBuf>) -> HashLookup<RecordBuf> {
    let a0 = Box::new(MockStream::new(0, s0)) as Box<dyn AlignmentStream<RecordBuf>>;
    let a1 = Box::new(MockStream::new(1, s1)) as Box<dyn AlignmentStream<RecordBuf>>;
    HashLookup::new(cfg(), smallvec![a0, a1], [None, None], [None, None]).unwrap()
}

fn r(name: &[u8], cigar: &str, md: &str) -> RecordBuf {
    create_record(name, cigar, &[], &[30u8; 20], md, false).unwrap()
}
fn u(name: &[u8]) -> RecordBuf {
    create_record(name, "", &vec![b'A'; 10], &[30u8; 10], "", false).unwrap()
}

fn disc0(h: &HashLookup<RecordBuf>) -> u64 {
    h.routing_counters[0]
}
fn out0(h: &HashLookup<RecordBuf>) -> u64 {
    h.routing_counters[1]
}
fn ambg0(h: &HashLookup<RecordBuf>) -> u64 {
    h.routing_counters[2]
}
fn disc1(h: &HashLookup<RecordBuf>) -> u64 {
    h.routing_counters[4]
}
fn out1(h: &HashLookup<RecordBuf>) -> u64 {
    h.routing_counters[5]
}
fn ambg1(h: &HashLookup<RecordBuf>) -> u64 {
    h.routing_counters[6]
}

struct Row {
    label: &'static str,
    s0: Vec<RecordBuf>,
    s1: Vec<RecordBuf>,
    out0: u64,
    out1: u64,
    disc0: u64,
    disc1: u64,
    ambg0: u64,
    ambg1: u64,
}

#[test]
fn hash_lookup_table() {
    let cases: &[Row] = &[
        // -- EarlyKind::AllPerfect / AllUnmapped --------------------------
        Row {
            label: "s0 perfect s1 imperfect → s0 wins (AllPerfect vs Scoring)",
            s0: vec![r(b"R1", "10M", "10")],
            s1: vec![r(b"R1", "5M5S", "5")],
            out0: 1,
            out1: 0,
            disc0: 0,
            disc1: 1,
            ambg0: 0,
            ambg1: 0,
        },
        Row {
            label: "s1 perfect s0 imperfect → s1 wins",
            s0: vec![r(b"R1", "5M5S", "5")],
            s1: vec![r(b"R1", "10M", "10")],
            out0: 0,
            out1: 1,
            disc0: 1,
            disc1: 0,
            ambg0: 0,
            ambg1: 0,
        },
        Row {
            label: "both perfect → ambiguous",
            s0: vec![r(b"R1", "10M", "10")],
            s1: vec![r(b"R1", "10M", "10")],
            out0: 0,
            out1: 0,
            disc0: 0,
            disc1: 0,
            ambg0: 1,
            ambg1: 1,
        },
        Row {
            label: "both unmapped → ambiguous (AllUnmapped vs AllUnmapped)",
            s0: vec![u(b"R1")],
            s1: vec![u(b"R1")],
            out0: 0,
            out1: 0,
            disc0: 0,
            disc1: 0,
            ambg0: 1,
            ambg1: 1,
        },
        Row {
            label: "s0 unmapped s1 mapped → s1 wins",
            s0: vec![u(b"R1")],
            s1: vec![r(b"R1", "10M", "10")],
            out0: 0,
            out1: 1,
            disc0: 1,
            disc1: 0,
            ambg0: 0,
            ambg1: 0,
        },
        Row {
            label: "s0 mapped s1 unmapped → s0 wins",
            s0: vec![r(b"R1", "10M", "10")],
            s1: vec![u(b"R1")],
            out0: 1,
            out1: 0,
            disc0: 0,
            disc1: 1,
            ambg0: 0,
            ambg1: 0,
        },
        // -- Pre-assess (Tier 2.5) -----------------------------------------
        Row {
            label: "s0 more matches (9 vs 7) → s0 wins via pre-assess",
            s0: vec![r(b"R1", "10M", "9A0")],  // 9 matches
            s1: vec![r(b"R1", "10M", "6AAA")], // 7 matches
            out0: 1,
            out1: 0,
            disc0: 0,
            disc1: 1,
            ambg0: 0,
            ambg1: 0,
        },
        // -- Out-of-order streams (hash table must handle) -----------------
        Row {
            label: "streams arrive out of fragment order",
            s0: vec![r(b"R1", "10M", "10"), r(b"R2", "5S5M", "5")],
            s1: vec![r(b"R2", "10M", "10"), r(b"R1", "5S5M", "5")],
            out0: 1,
            out1: 1,
            disc0: 1,
            disc1: 1,
            ambg0: 0,
            ambg1: 0,
        },
        // -- Paired-end: both mates must be present ------------------------
        Row {
            label: "paired-end both perfect → ambiguous",
            s0: vec![
                create_record(b"R1", "10M", &[], &[30u8; 10], "10", false)
                    .map(|mut r| {
                        *r.flags_mut() = Flags::from_bits(0x41).unwrap();
                        r
                    })
                    .unwrap(),
                create_record(b"R1", "10M", &[], &[30u8; 10], "10", true)
                    .map(|mut r| {
                        *r.flags_mut() = Flags::from_bits(0x81).unwrap();
                        r
                    })
                    .unwrap(),
            ],
            s1: vec![
                create_record(b"R1", "10M", &[], &[30u8; 10], "10", false)
                    .map(|mut r| {
                        *r.flags_mut() = Flags::from_bits(0x41).unwrap();
                        r
                    })
                    .unwrap(),
                create_record(b"R1", "10M", &[], &[30u8; 10], "10", true)
                    .map(|mut r| {
                        *r.flags_mut() = Flags::from_bits(0x81).unwrap();
                        r
                    })
                    .unwrap(),
            ],
            out0: 0,
            out1: 0,
            disc0: 0,
            disc1: 0,
            ambg0: 2,
            ambg1: 2,
        },
        // -- Suffix stripping ----------------------------------------------
        Row {
            label: "suffix stripped /1 /2 matched as same fragment",
            s0: vec![create_record(b"R1/1", "10M", &[], &[30u8; 10], "10", false).unwrap()],
            s1: vec![create_record(b"R1/2", "5S5M", &[], &[30u8; 10], "5", false).unwrap()],
            // NOTE: this row needs strip_read_suffix=True; the default make() uses False.
            // Tested separately below.
            out0: 0,
            out1: 0,
            disc0: 0,
            disc1: 0,
            ambg0: 0,
            ambg1: 0,
        }, // placeholder — see below
    ];
    for c in cases {
        // Skip the suffix-stripping placeholder row.
        if c.label.starts_with("suffix") {
            continue;
        }
        let mut h = make(c.s0.clone(), c.s1.clone());
        h.process().unwrap();
        assert_eq!(out0(&h), c.out0, "[{}] out[0]", c.label);
        assert_eq!(out1(&h), c.out1, "[{}] out[1]", c.label);
        assert_eq!(disc0(&h), c.disc0, "[{}] disc[0]", c.label);
        assert_eq!(disc1(&h), c.disc1, "[{}] disc[1]", c.label);
        assert_eq!(ambg0(&h), c.ambg0, "[{}] ambig[0]", c.label);
        assert_eq!(ambg1(&h), c.ambg1, "[{}] ambig[1]", c.label);
    }
}

#[test]
fn hash_lookup_suffix_stripping() {
    let cfg = Config {
        strip_read_suffix: StripReadSuffix::True,
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Config::default()
    };
    let a0 = Box::new(MockStream::new(
        0,
        vec![create_record(b"R1/1", "10M", &[], &[30u8; 10], "10", false).unwrap()],
    )) as Box<dyn AlignmentStream<RecordBuf>>;
    let a1 = Box::new(MockStream::new(
        1,
        vec![create_record(b"R1/2", "5M5S", &[], &[30u8; 10], "5", false).unwrap()],
    )) as Box<dyn AlignmentStream<RecordBuf>>;
    let mut h = HashLookup::new(cfg, smallvec![a0, a1], [None, None], [None, None]).unwrap();
    h.process().unwrap();
    assert_eq!(out0(&h), 1, "s0 perfect should win");
    assert_eq!(disc1(&h), 1, "s1 should be discarded");
}
