use super::*;
use crate::{
    alignment::{fragment::SimpleRec, MdCigFlags},
    config::args::ScoringArgs,
    filter_algorithm::line_by_line::Scratch,
    penalty::Penalty,
    tests::create_record,
    variant::{Eval, Variant},
    Error,
};
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{smallvec, SmallVec};

// ---------------------------------------------------------------------------
// Test penalty helpers
// ---------------------------------------------------------------------------

/// Setup penalties: match=0, mismatch=-1, gap_open=-2, gap_extend=-0.5.
/// Quality-independent -- use for unit-testing score logic without Q noise.
pub(crate) fn setup_penalties() -> Penalty {
    let c = ScoringArgs::default();
    let mut p = c.to_penalty();
    p.log_likelihood_match = [0.0; 93];
    p.log_likelihood_mismatch = [-1.0; 93];
    p.gap_open = -2.0;
    p.gap_extend = -0.5;
    p
}

/// Real penalties: matches near 0, mismatch = -(q/10) * scaling.
/// Use for quality-edge-case tests.
fn real_penalties() -> Penalty {
    ScoringArgs {
        mismatch_penalty: 4.0,
        gap_open: 6.0,
        gap_extend: 1.0,
        ..ScoringArgs::default()
    }
    .to_penalty()
}

fn score_one(cigar: &str, md: &str, qual: &[u8], pen: &Penalty) -> f64 {
    let rec = create_record(b"r", cigar, &[], qual, md, false).unwrap();
    let flags = rec.flags();
    let mut frag = Fragment::new(
        pen,
        smallvec![&rec],
        smallvec![MdCigFlags::try_from_record(&rec, &flags, false).unwrap()],
    )
    .unwrap();
    let mut dvnt = smallvec![smallvec![]];
    let mut scratch = Scratch::new();
    frag.score(&mut scratch, &mut dvnt).unwrap()
}

// ---------------------------------------------------------------------------
// NW scoring -- table-driven
// ---------------------------------------------------------------------------

#[test]
fn scoring_table() {
    // (label, cigar, md, score_expected_with_flat_penalties)
    // Derivations:
    //   All match: 0.0
    //   1 mismatch in 10M: -1.0
    //   2 mismatch in 10M: -2.0
    //   5S5M (5 clips): 5*(-1) + 5*0 = -5.0
    //   5M1D5M: 10*0 + gap_open + 1*gap_extend = -2.5
    //   5M2I5M: 10*0 + gap_open + 2*gap_extend = -3.0
    struct Row {
        cigar: &'static str,
        md: &'static str,
        want: f64,
    }
    impl Row {
        fn new(cigar: &'static str, md: &'static str, want: f64) -> Self {
            Row { cigar, md, want }
        }
    }

    let q30 = [30u8; 15]; // generous budget
    let flat = setup_penalties();
    let cases: &[Row] = &[
        Row::new("10M", "10", 0.0),
        Row::new("10M", "9A0", -1.0),
        Row::new("10M", "4A4A0", -2.0),
        Row::new("5S5M", "5", -5.0),
        Row::new("5M1D5M", "5^A5", -2.5),
        Row::new("5M2I3M", "10", -3.0),
        Row::new("5M1D2I2M", "5^A5", -5.5), // gap_open*2 + 1*ext + 2*ext = -2-0.5-2-1 = -5.5
    ];

    crate::tests::common::run_collecting(
        cases,
        |c| format!("cigar={} md={}", c.cigar, c.md),
        |c| {
            let got = score_one(c.cigar, c.md, &q30[..c.cigar.len().max(10)], &flat);
            if (got - c.want).abs() >= 1e-9 {
                Err(format!("want {} got {}", c.want, got))
            } else {
                Ok(())
            }
        },
    );
}

/// A single base mismatch at quality 5 should cost -0.5 with real penalties
/// (mismatch = -(q/10) * 1.0).  If stream A has mismatch at q=5 and stream B
/// at q=30, stream A pays less penalty and wins, even though both have the same
/// CIGAR/MD profile.  This is the intended quality-weighted behaviour.
#[test]
fn quality_determines_winner_at_equal_cigar_md() {
    let pen = real_penalties();
    let score_q5 = score_one("5M", "4A0", &[5u8; 5], &pen); // mismatch at q=5  -> -0.5
    let score_q30 = score_one("5M", "4A0", &[30u8; 5], &pen); // mismatch at q=30 -> -3.0
    assert!(
        score_q5 > score_q30,
        "lower-quality mismatch should score higher (less penalised): q5={score_q5} q30={score_q30}"
    );
}

// ---------------------------------------------------------------------------
// Weighted interval scheduling -- table-driven
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct V {
    pos: usize,
    ref_len: usize,
    alt_len: usize,
    delta: f64,
}

struct MockVariant {
    pos: usize,
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
}
impl Variant for MockVariant {
    fn pos(&self) -> usize {
        self.pos
    }
    fn ref_allele(&self) -> &[u8] {
        &self.ref_a
    }
    fn alt_allele(&self) -> &[u8] {
        &self.alt_a
    }
    fn p_variant(&self) -> f64 {
        0.1
    } // unused in these tests
}

fn mk_eval(v: V) -> Eval<'static> {
    let mv: &'static MockVariant = Box::leak(Box::new(MockVariant {
        pos: v.pos,
        ref_a: vec![b'A'; v.ref_len],
        alt_a: vec![b'A'; v.alt_len],
    }));
    let mut e = Eval::new();
    e.set_variant(mv);
    e.update(0.0, v.delta);
    e
}

#[test]
fn wis_table() {
    struct Row {
        label: &'static str,
        variants: &'static [V],
        want: f64,
    }

    // Interval [pos, pos + max(ref_len, alt_len))
    let cases: &[Row] = &[
        Row {
            label: "empty",
            variants: &[],
            want: 0.0,
        },
        Row {
            label: "single positive",
            variants: &[V {
                pos: 10,
                ref_len: 1,
                alt_len: 1,
                delta: 7.5,
            }],
            want: 7.5,
        },
        Row {
            label: "single zero",
            variants: &[V {
                pos: 10,
                ref_len: 1,
                alt_len: 1,
                delta: 0.0,
            }],
            want: 0.0,
        },
        Row {
            label: "single negative",
            variants: &[V {
                pos: 10,
                ref_len: 1,
                alt_len: 1,
                delta: -5.,
            }],
            want: 0.0,
        },
        // [10,11) and [20,21) -- non-overlapping -> sum
        Row {
            label: "two non-overlapping",
            variants: &[
                V {
                    pos: 10,
                    ref_len: 1,
                    alt_len: 1,
                    delta: 3.0,
                },
                V {
                    pos: 20,
                    ref_len: 1,
                    alt_len: 1,
                    delta: 2.0,
                },
            ],
            want: 5.0,
        },
        // [10,20) and [15,25) -- overlapping; pick the larger
        Row {
            label: "two overlapping pick best",
            variants: &[
                V {
                    pos: 10,
                    ref_len: 10,
                    alt_len: 10,
                    delta: 5.0,
                },
                V {
                    pos: 15,
                    ref_len: 10,
                    alt_len: 10,
                    delta: 8.0,
                },
            ],
            want: 8.0,
        },
        // [10,20) and [15,17) and [40,41) -- inner + disjoint
        // Best chain: inner(20.0) + disjoint(3.0) = 23.0 vs outer(5.0) + disjoint(3.0)
        Row {
            label: "nested variants best chain",
            variants: &[
                V {
                    pos: 10,
                    ref_len: 20,
                    alt_len: 20,
                    delta: 5.0,
                },
                V {
                    pos: 15,
                    ref_len: 2,
                    alt_len: 2,
                    delta: 20.0,
                },
                V {
                    pos: 40,
                    ref_len: 1,
                    alt_len: 1,
                    delta: 3.0,
                },
            ],
            want: 23.0,
        },
        // Touching but non-overlapping: [10,20) [20,30) -> both count
        Row {
            label: "touching boundaries non-overlapping",
            variants: &[
                V {
                    pos: 10,
                    ref_len: 10,
                    alt_len: 10,
                    delta: 5.0,
                },
                V {
                    pos: 20,
                    ref_len: 10,
                    alt_len: 10,
                    delta: 7.0,
                },
            ],
            want: 12.0,
        },
        // Insertion: end = pos + alt_len when alt > ref
        Row {
            label: "insertion span respected",
            variants: &[
                V {
                    pos: 10,
                    ref_len: 1,
                    alt_len: 20,
                    delta: 10.0,
                }, // [10, 30)
                V {
                    pos: 25,
                    ref_len: 1,
                    alt_len: 1,
                    delta: 50.0,
                }, // overlaps insertion span
            ],
            want: 50.0,
        },
    ];
    for c in cases {
        let mut dvnt: FragEvalVec<'_> =
            smallvec![c.variants.iter().copied().map(mk_eval).collect()];
        let mut dp = smallvec![];
        let got = wis_max_rescue_delta(&mut dvnt, &mut dp);
        assert!(
            (got - c.want).abs() < 1e-9,
            "[{}] want {} got {}",
            c.label,
            c.want,
            got
        );
    }
}

// ---------------------------------------------------------------------------
// Variant rescue: p_variant threshold
// ---------------------------------------------------------------------------

#[test]
fn variant_rescue_p_variant_table() {
    // delta = alt_score - incurred = (p*lm + (1-p)*lmm) - ((1-p)*lm + p*lmm)
    //       = (2p-1)*(lm-lmm)
    // With flat penalties (lm=0, lmm=-1): delta = (2p-1)*1 = 2p-1
    // p < 0.5 -> negative -> no rescue; p=0.5 -> 0; p > 0.5 -> positive.
    struct Row {
        label: &'static str,
        p: f64,
        want_positive: bool,
    }
    let cases: &[Row] = &[
        Row {
            label: "p=0.1 no rescue",
            p: 0.1,
            want_positive: false,
        },
        Row {
            label: "p=0.49 no rescue",
            p: 0.49,
            want_positive: false,
        },
        Row {
            label: "p=0.5  zero",
            p: 0.5,
            want_positive: false,
        }, // exactly zero
        Row {
            label: "p=0.51 rescue",
            p: 0.51,
            want_positive: true,
        },
        Row {
            label: "p=1.0  rescue",
            p: 1.0,
            want_positive: true,
        },
    ];

    // Read: "AAGAA" (base at position 2 matches alt 'G')
    let pen = setup_penalties();
    for c in cases {
        let rec = create_record(b"r", "5M", b"AAGAA", &[30u8; 5], "2G2", false).unwrap();
        let flags = rec.flags();
        let _mv: &'static _ = Box::leak(Box::new(MockVariant {
            pos: 2,
            ref_a: vec![b'A'],
            alt_a: vec![b'G'],
        }));
        // Override p_variant via a wrapper; easier to leak a struct with p embedded.
        struct PVariant {
            inner: MockVariant,
            p: f64,
        }
        impl Variant for PVariant {
            fn pos(&self) -> usize {
                self.inner.pos()
            }
            fn ref_allele(&self) -> &[u8] {
                self.inner.ref_allele()
            }
            fn alt_allele(&self) -> &[u8] {
                self.inner.alt_allele()
            }
            fn p_variant(&self) -> f64 {
                self.p
            }
        }
        let pv: &'static _ = Box::leak(Box::new(PVariant {
            inner: MockVariant {
                pos: 2,
                ref_a: vec![b'A'],
                alt_a: vec![b'G'],
            },
            p: c.p,
        }));
        let mut ev = Eval::new();
        ev.set_variant(pv as &dyn Variant);

        let mut dvnt: FragEvalVec<'_> = smallvec![smallvec![ev]];
        let mut scratch = Scratch::new();
        let mut frag = Fragment::new(
            &pen,
            smallvec![&rec],
            smallvec![MdCigFlags::try_from_record(&rec, &flags, false).unwrap()],
        )
        .unwrap();
        let _ = frag.score(&mut scratch, &mut dvnt).unwrap();
        let delta = scratch.last_variant_delta;
        if c.want_positive {
            assert!(
                delta > 0.0,
                "[{}] expected positive delta, got {delta}",
                c.label
            );
        } else {
            assert!(
                delta <= 0.0,
                "[{}] expected non-positive delta, got {delta}",
                c.label
            );
        }
    }
}

struct TestVariant {
    pos: usize,
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
    p_variant: f64,
}
impl TestVariant {
    fn new(pos: usize, len: usize, p_variant: f64) -> Self {
        TestVariant {
            pos,
            ref_a: vec![b'A'; len],
            alt_a: vec![b'G'; len],
            p_variant,
        }
    }
}
impl Variant for TestVariant {
    fn pos(&self) -> usize {
        self.pos
    }
    fn ref_allele(&self) -> &[u8] {
        &self.ref_a
    }
    fn alt_allele(&self) -> &[u8] {
        &self.alt_a
    }
    fn p_variant(&self) -> f64 {
        self.p_variant
    }
}

fn make_eval<'a>(v: &'a dyn Variant) -> Eval<'a> {
    let mut e = Eval::new();
    e.set_variant(v);
    e
}
/* conversion to Record not possible
#[test]
fn test_simple_rec_for_bam_record() -> Result<(), Error> {
    use noodles::bam;
    use noodles::sam::Header;

    let header = Header::default();
    let mut buf = RecordBuf::default();
    *buf.reference_sequence_id_mut() = Some(42);
    *buf.quality_scores_mut() =
        noodles::sam::alignment::record_buf::QualityScores::from(vec![30, 31, 32]);

    // Convert to a BAM record
    let bam_record = bam::Record::try_from_alignment_record(&header, &buf)?;

    assert_eq!(bam_record.ref_seq_id().unwrap()?, 42);

    assert_eq!(bam_record.quality_at(0), Some(30));
    assert_eq!(bam_record.quality_at(2), Some(32));
    assert_eq!(bam_record.quality_at(3), None);

    let converted = bam_record.as_record_buf(&header)?;
    assert_eq!(converted.reference_sequence_id(), Some(42));

    Ok(())
}*/

#[test]
fn test_fragment_requires_revcmp() -> Result<(), Error> {
    let rec_fwd = create_record(b"read1", "5M", b"AAAAA", &[30; 5], "5", false)?;
    let mut rec_rev = create_record(b"read1", "5M", b"AAAAA", &[30; 5], "10", false)?;
    *rec_rev.flags_mut() = Flags::from_bits(0x10).unwrap(); // Reverse complemented

    let flags_fwd = rec_fwd.flags();
    let flags_rev = rec_rev.flags();

    let p = setup_penalties();
    let mut fragment = Fragment::new(
        &p,
        smallvec![&rec_fwd, &rec_rev],
        smallvec![
            MdCigFlags::try_from_record(&rec_fwd, &flags_fwd, false)?,
            MdCigFlags::try_from_record(&rec_rev, &flags_rev, false)?
        ],
    )?;

    // fragment is initialized with seg_i = 0 (the forward segment)

    assert!(!fragment.requires_revcmp(0)); // Same segment -> false

    assert!(fragment.requires_revcmp(1)); // Other segment has different ori -> true

    // Move to the revcmp segment
    fragment.seg_i = 1;
    assert!(fragment.requires_revcmp(0)); // Other segment has different ori -> true
    assert!(!fragment.requires_revcmp(1)); // Same segment -> false

    Ok(())
}

#[test]
fn test_score_variants_in_window_boundaries() -> Result<(), Error> {
    // Start the read at position 1. Spans reference [1, 11).
    let rec = create_record(b"read1", "10M", &[b'A'; 10], &[30; 10], "1", false)?;
    let flags = rec.flags();
    let p = setup_penalties();
    let fragment = Fragment::new(
        &p,
        smallvec![&rec],
        smallvec![MdCigFlags::try_from_record(&rec, &flags, false)?],
    )?;

    // We create three variants around the right edge of a window ending at 6:
    // v1: fully inside the window (start=3, ref_end=4)
    // v2: overlapping the right edge (start=5, ref_end=7)
    // v3: starts exactly on the right edge (start=6, ref_end=7) -> should NOT be processed
    let v1 = TestVariant::new(3, 1, 0.1);
    let v2 = TestVariant::new(5, 2, 0.1);
    let v3 = TestVariant::new(6, 1, 0.1);

    let mut dvnt: FragEvalVec =
        smallvec![smallvec![make_eval(&v1), make_eval(&v2), make_eval(&v3)]];
    let mut finished: FragEvalVec = smallvec![smallvec![]];
    let mut scratch = Scratch::new();

    // Create a WindowCtx that covers [1, 6)
    let ctx = WindowCtx::new(0, 0, 1, 6, 0.0);

    fragment.evaluate_variants_in_window(&mut scratch, &mut dvnt, &mut finished, ctx)?;

    // v1 is fully processed (ref_end 4 <= 6). It gets moved to `finished`.
    assert_eq!(finished[0].len(), 1);
    assert_eq!(finished[0][0].start(), 3);

    // dvnt should now hold v2 and v3.
    // v2 was partially processed (start 5 < 6, but ref_end 7 > 6), so it remains in `dvnt`.
    // v3 was NOT processed because its start (6) is NOT < ctx.ref_end (6).
    assert_eq!(dvnt[0].len(), 2);
    assert_eq!(dvnt[0][0].start(), 5);
    assert_eq!(dvnt[0][1].start(), 6);

    Ok(())
}

#[test]
fn maximize_delta_exact_zero_mutant() {
    let mut dvnt: FragEvalVec = smallvec![smallvec![mk_eval2(10, 1, 1, 0.0)]];
    let mut dp = SmallVec::new();

    // If it were >=, the variant with delta 0.0 would be included and processed.
    assert_eq!(wis_max_rescue_delta(&mut dvnt, &mut dp), 0.0);
    assert!(
        dp.is_empty(),
        "dp should be empty because variants with <= 0 delta should be filtered out"
    );
}

#[test]
fn test_simple_rec_for_sam_record_buf() {
    let mut buf = RecordBuf::default();
    *buf.reference_sequence_id_mut() = Some(15);

    assert_eq!(buf.ref_seq_id().unwrap().unwrap(), 15);
}

#[test]
fn test_stitched_fragment_creation() -> Result<(), Error> {
    let record1 = create_record(b"read1", "5M3S", &[b'A'; 8], &[30; 8], "5", false)?;
    let record2 = create_record(b"read1", "4M4S", &[b'A'; 8], &[30; 8], "4", false)?;
    let flags1 = record1.flags();
    let flags2 = record2.flags();

    let records: SmallVec<[&RecordBuf; READ_CT]> = smallvec![&record1, &record2];
    let p = setup_penalties();
    let mut md_cig_flags: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();
    md_cig_flags.push(MdCigFlags::try_from_record(&record1, &flags1, false)?);
    md_cig_flags.push(MdCigFlags::try_from_record(&record2, &flags2, false)?);
    let _stitched = Fragment::new(&p, records, md_cig_flags)?;
    Ok(())
}
#[test]
fn test_complement() {
    assert_eq!(complement(b'A'), b'T');
    assert_eq!(complement(b'C'), b'G');
    assert_eq!(complement(b'G'), b'C');
    assert_eq!(complement(b'T'), b'A');
    assert_eq!(complement(b'N'), b'N');
}

#[test]
fn test_q() -> Result<(), Error> {
    let record = create_record(
        b"read1",
        "5M",
        &[b'A'; 5],
        &[30, 31, 32, 33, 34],
        "5",
        false,
    )?;
    let flags = record.flags();
    let seg = smallvec![&record];
    let md_cig_flags = smallvec![MdCigFlags::try_from_record(&record, &flags, false)?];
    let p = setup_penalties();
    let fragment = Fragment::new(&p, seg, md_cig_flags)?;
    assert_eq!(fragment.q(0, 0)?, 30);
    assert_eq!(fragment.q(0, 4)?, 34);
    assert!(fragment.q(0, 5).is_err());
    Ok(())
}
fn mk_eval2(pos: usize, ref_len: usize, alt_len: usize, delta: f64) -> Eval<'static> {
    let v: &'static TestVariant = Box::leak(Box::new(TestVariant {
        pos,
        ref_a: vec![b'A'; ref_len],
        alt_a: vec![b'A'; alt_len],
        p_variant: 0.1,
    }));
    let mut eval = Eval::new();
    eval.set_variant(v);
    eval.update(0.0, delta); // delta() = alt_score - incurred = delta - 0 = delta
    eval
}

#[test]
fn test_test_variant_p_variant_reflects_constructed_field() {
    let v = TestVariant {
        pos: 0,
        ref_a: b"A".to_vec(),
        alt_a: b"G".to_vec(),
        p_variant: 1.0,
    };
    assert_eq!(v.p_variant(), 1.0);
    let v2 = TestVariant {
        pos: 0,
        ref_a: b"A".to_vec(),
        alt_a: b"G".to_vec(),
        p_variant: 0.25,
    };
    assert_eq!(v2.p_variant(), 0.25);
}
