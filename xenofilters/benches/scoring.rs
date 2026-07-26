//! Microbenchmarks for the per-fragment scoring hot path.
//! Run: cargo bench --bench scoring

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use smallvec::smallvec;
use xenofilters::{
    alignment::{fragment::Fragment, pre_assess::pre_assess_alignments, MdCigFlags},
    filter_algorithm::line_by_line::core::Scratch,
    penalty::Penalty,
    tests::create_record,
};

fn bench_penalties() -> Penalty {
    xenofilters::config::args::ScoringArgs {
        mismatch_penalty: 4.0,
        gap_open: 6.0,
        gap_extend: 1.0,
        chimeric_junction_bases: 20,
        ..Default::default()
    }
    .to_penalty()
}

// ---------------------------------------------------------------------------
// ScoreOpIter walk — Tier 2.5a inner loop
// ---------------------------------------------------------------------------

fn bench_score_op_iter(c: &mut Criterion) {
    // Representative cases: perfect, 1 mismatch, 10% mismatch, softclip, deletion.
    let cases: &[(&str, &str)] = &[
        ("perfect_150M", "150"),
        ("1mm_150M", "74A75"),
        ("10pct_mm_150M", "0A14A14A14A14A14A14A14A14A14A14"),
        ("softclip_10S140M", "140"),
        ("deletion_70M5D75M", "70^ACGTA75"),
    ];
    let pen = bench_penalties();
    let mut g = c.benchmark_group("ScoreOpIter");
    for &(label, md) in cases {
        let cigar = if label.starts_with("softclip") {
            "10S140M"
        } else if label.starts_with("deletion") {
            "70M5D75M"
        } else {
            "150M"
        };
        let qual = vec![30u8; 150];
        let rec = create_record(b"r", cigar, &[], &qual, md, false).unwrap();
        let flags = rec.flags();
        g.bench_with_input(BenchmarkId::new("score_one", label), label, |b, _| {
            b.iter(|| {
                let mcf = MdCigFlags::try_from_record(&rec, &flags, false).unwrap();
                let mut frag =
                    Fragment::new(&pen, smallvec![black_box(&rec)], smallvec![mcf]).unwrap();
                let mut dvnt = smallvec![smallvec![]];
                let mut scratch = Scratch::new();
                black_box(frag.score(&mut scratch, &mut dvnt).unwrap())
            })
        });
    }
    g.finish();
}

// ---------------------------------------------------------------------------
// pre_assess_alignments — Tier 2.5 dispatch
// ---------------------------------------------------------------------------

fn bench_pre_assess(c: &mut Criterion) {
    // Pairs: (a_cigar, a_md, b_cigar, b_md, expected_outcome)
    let cases: &[(&str, &str, &str, &str)] = &[
        ("150M", "150", "150M", "148AA"),   // a dominates → EarlyDecision
        ("150M", "148AA", "150M", "148AA"), // equal → Equal
        ("150M", "74A75", "150M", "74C75"), // incomparable → FullScoring
    ];
    let qual = vec![30u8; 150];
    let mut g = c.benchmark_group("pre_assess_alignments");
    let pen = bench_penalties();
    for (i, &(ac, am, bc, bm)) in cases.iter().enumerate() {
        let ra = create_record(b"r", ac, &[], &qual, am, false).unwrap();
        let rb = create_record(b"r", bc, &[], &qual, bm, false).unwrap();
        let fa = ra.flags();
        let fb = rb.flags();
        g.bench_with_input(BenchmarkId::from_parameter(i), &i, |b, _| {
            b.iter(|| {
                let ma = smallvec![Some(MdCigFlags::try_from_record(&ra, &fa, false).unwrap())];
                let mb = smallvec![Some(MdCigFlags::try_from_record(&rb, &fb, false).unwrap())];
                black_box(pre_assess_alignments(&ma, &mb, &pen, 0.5))
            })
        });
    }
    g.finish();
}

// ---------------------------------------------------------------------------
// wis_max_rescue_delta — WIS scheduling
// ---------------------------------------------------------------------------

fn bench_wis(c: &mut Criterion) {
    use xenofilters::{
        alignment::fragment::wis_max_rescue_delta,
        variant::{eval::Eval, Variant},
    };

    struct V {
        pos: usize,
        len: usize,
        delta: f64,
    }
    impl Variant for V {
        fn pos(&self) -> usize {
            self.pos
        }
        fn ref_allele(&self) -> &[u8] {
            b"A"
        }
        fn alt_allele(&self) -> &[u8] {
            b"G"
        }
        fn p_variant(&self) -> f64 {
            0.9
        }
    }

    // Build N non-overlapping variants (worst-case WIS is O(N log N)).
    let make_dvnt = |n: usize| {
        let variants: Vec<&'static V> = (0..n)
            .map(|i| {
                let v: &'static V = Box::leak(Box::new(V {
                    pos: i * 10,
                    len: 5,
                    delta: 2.0,
                }));
                v
            })
            .collect();
        let evals: smallvec::SmallVec<[Eval<'_>; 4]> = variants
            .iter()
            .map(|v| {
                let mut e = Eval::new();
                e.set_variant(*v as &dyn Variant);
                e.update(0.0, 2.0);
                e
            })
            .collect();
        smallvec![evals]
    };

    let mut g = c.benchmark_group("wis_max_rescue_delta");
    for n in [1, 10, 50, 200] {
        g.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            b.iter(|| {
                let mut dvnt = make_dvnt(n);
                let mut dp = smallvec![];
                black_box(wis_max_rescue_delta(&mut dvnt, &mut dp, 0.0))
            })
        });
    }
    g.finish();
}

criterion_group!(benches, bench_score_op_iter, bench_pre_assess, bench_wis);
criterion_main!(benches);
