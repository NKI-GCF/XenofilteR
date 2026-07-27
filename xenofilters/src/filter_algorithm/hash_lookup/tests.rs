use crate::tests::common::{r, u};
use crate::{
    aln_stream::AlignmentStream, aln_stream::tests::MockStream, config::run_config::RunConfig,
    filter_algorithm::hash_lookup::HashLookup, tests::create_record,
};
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::smallvec;

fn make(s0: Vec<RecordBuf>, s1: Vec<RecordBuf>, config: &RunConfig) -> HashLookup<RecordBuf> {
    let a0 = Box::new(MockStream::new(0, s0)) as Box<dyn AlignmentStream<RecordBuf>>;
    let a1 = Box::new(MockStream::new(1, s1)) as Box<dyn AlignmentStream<RecordBuf>>;
    HashLookup::<RecordBuf>::new(
        config,
        smallvec![a0, a1],
        [None, None],
        [None, None],
        [None, None],
    )
    .unwrap()
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
fn match_count_raw_fixture_md_strings_are_internally_consistent() {
    fn md_implied_length(md: &[u8]) -> usize {
        let mut len = 0usize;
        let mut i = 0;
        while i < md.len() {
            if md[i].is_ascii_digit() {
                let mut n = 0usize;
                while i < md.len() && md[i].is_ascii_digit() {
                    n = n * 10 + (md[i] - b'0') as usize;
                    i += 1;
                }
                len += n;
            } else {
                len += 1;
                i += 1;
            }
        }
        len
    }
    let cases: &[(&str, &[u8])] = &[("10M", b"9A0"), ("10M", b"7AAA")];
    for (cigar, md) in cases {
        let cigar_len: usize = cigar.trim_end_matches('M').parse().unwrap();
        assert_eq!(
            md_implied_length(md),
            cigar_len,
            "MD {:?} does not sum to CIGAR {cigar}'s length",
            std::str::from_utf8(md)
        );
    }
}

#[test]
fn hash_lookup_table() {
    let cases: &[Row] = &[
        // -- EarlyKind::AllPerfect / AllUnmapped --------------------------
        Row {
            label: "s0 perfect s1 imperfect -> s0 wins (AllPerfect vs Scoring)",
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
            label: "s1 perfect s0 imperfect -> s1 wins",
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
            label: "both perfect -> ambiguous",
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
            label: "both unmapped -> ambiguous (AllUnmapped vs AllUnmapped)",
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
            label: "s0 unmapped s1 mapped -> s1 wins",
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
            label: "s0 mapped s1 unmapped -> s0 wins",
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
            label: "s0 more matches (9 vs 7) -> s0 wins via pre-assess",
            s0: vec![r(b"R1", "10M", "9A0")],
            s1: vec![r(b"R1", "10M", "7AAA")],
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
            label: "paired-end both perfect -> ambiguous",
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
    ];

    crate::tests::common::run_collecting(
        cases,
        |c| c.label.to_string(),
        |c| {
            let config = RunConfig::default();
            let mut h = make(c.s0.clone(), c.s1.clone(), &config);

            if let Err(e) = h.process(&config) {
                return Err(format!("process() failed: {}", e));
            }
            if out0(&h) != c.out0 {
                return Err(format!("out[0] want {} got {}", c.out0, out0(&h)));
            }
            if out1(&h) != c.out1 {
                return Err(format!("out[1] want {} got {}", c.out1, out1(&h)));
            }
            if disc0(&h) != c.disc0 {
                return Err(format!("disc[0] want {} got {}", c.disc0, disc0(&h)));
            }
            if disc1(&h) != c.disc1 {
                return Err(format!("disc[1] want {} got {}", c.disc1, disc1(&h)));
            }
            if ambg0(&h) != c.ambg0 {
                return Err(format!("ambig[0] want {} got {}", c.ambg0, ambg0(&h)));
            }
            if ambg1(&h) != c.ambg1 {
                return Err(format!("ambig[1] want {} got {}", c.ambg1, ambg1(&h)));
            }
            Ok(())
        },
    );
}
