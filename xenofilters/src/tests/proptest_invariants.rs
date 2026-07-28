//! Property-based invariants for the tournament/scoring pipeline.
//!
//! Scope: these tests operate at the `RecordBuf` / `MockStream` /
//! `LineByLine` level, the same abstraction level as the existing unit
//! tests (`src/tests/exploratory.rs`, `ordering/tests.rs`). They do not
//! exercise real BAM/CRAM I/O.
//!
//! Each `proptest!` block targets exactly one invariant so failures shrink
//! to the smallest counterexample for that invariant specifically, rather
//! than a shared mega-test.

use crate::{
    alignment::fragment::wis_max_rescue_delta,
    alignment::FragmentState,
    config::{args::ScoringArgs, run_config::RunConfig},
    filter_algorithm::{
        collated::CollatedMatcher,
        hash_lookup::HashLookup,
        line_by_line::{
            core::FragmentBuffer, detect_chimeric_event, ChimericDecision, ChimericThresholds,
            LineByLine, MAX_STREAMS,
        },
    },
    tests::{create_record, MockStream},
    variant::{Eval, Variant},
    Error,
};
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::RecordBuf;
use proptest::prelude::*;
use smallvec::{smallvec, SmallVec};

// ===========================================================================
// Custom strategies
// ===========================================================================

/// A CIGAR/MD pair that is internally consistent for `MdCigFlags`/`ScoreOpIter`:
/// the CIGAR is a single `NM` op, and the MD string's implied length
/// (match-run digits + one letter per mismatch) equals `N`.
///
/// Deliberately avoids leading/duplicate zero-length digit runs (see module
/// doc: `next_md_base` mis-scores a leading `"0"` as a Match). Indels are out
/// of scope here; other tests already cover indel CIGAR/MD combinations with
/// fixed patterns.
#[derive(Debug, Clone)]
struct MdCigFixture {
    cigar: String,
    md: String,
    read_len: usize,
    /// number of bases classified as Match by this MD string (informational,
    /// used to derive "more matches" vs "fewer matches" fixtures)
    match_bases: usize,
}

fn md_cig_strategy(max_len: usize) -> impl Strategy<Value = MdCigFixture> {
    (2..=max_len.max(2)).prop_flat_map(|len| {
        // Force position 0 to be a match, avoiding the leading-zero quirk.
        prop::collection::vec(any::<bool>(), len - 1).prop_map(move |rest| {
            let mismatches: Vec<bool> = std::iter::once(false).chain(rest).collect();

            let mut md = String::new();
            let mut run = 0usize;
            let mut prev_mismatch = false;
            let mut match_bases = 0usize;
            for &is_mm in &mismatches {
                if is_mm {
                    if !prev_mismatch {
                        md.push_str(&run.to_string());
                        run = 0;
                    }
                    md.push('A');
                    prev_mismatch = true;
                } else {
                    run += 1;
                    match_bases += 1;
                    prev_mismatch = false;
                }
            }
            md.push_str(&run.to_string());

            MdCigFixture {
                cigar: format!("{len}M"),
                md,
                read_len: len,
                match_bases,
            }
        })
    })
}

impl MdCigFixture {
    fn record(&self, qname: &[u8], qual: u8, is_rev: bool) -> Result<RecordBuf, Error> {
        let q = vec![qual; self.read_len];
        create_record(qname, &self.cigar, &[], &q, &self.md, is_rev)
    }
}

/// `Arbitrary`-style wrapper demonstrating how to keep generated flags
/// within valid domain constraints: a primary/secondary/supplementary
/// record can never simultaneously be secondary AND supplementary AND
/// unmapped in a way that would make downstream code panic, so we encode
/// the *legal* combinations directly in the strategy rather than filtering
/// post-hoc with `prop_filter` (which wastes shrinking effort).
#[derive(Debug, Clone, Copy)]
enum RecRole {
    PrimaryMapped { reverse: bool },
    PrimaryUnmapped,
    Secondary,
}

fn rec_role_strategy() -> impl Strategy<Value = RecRole> {
    prop_oneof![
        3 => any::<bool>().prop_map(|reverse| RecRole::PrimaryMapped { reverse }),
        1 => Just(RecRole::PrimaryUnmapped),
        1 => Just(RecRole::Secondary),
    ]
}

fn n_stream_aln(
    per_stream: Vec<Vec<RecordBuf>>,
) -> SmallVec<[Box<dyn crate::aln_stream::AlignmentStream<RecordBuf>>; 2]> {
    per_stream
        .into_iter()
        .enumerate()
        .map(|(i, recs)| {
            Box::new(MockStream::new(i, recs))
                as Box<dyn crate::aln_stream::AlignmentStream<RecordBuf>>
        })
        .collect()
}

fn default_run_config(add_decision_tag: bool, quiet: bool, no_program_line: bool) -> RunConfig {
    RunConfig {
        io: crate::config::args::IoArgs {
            add_decision_tag,
            quiet,
            no_program_line,
            ..Default::default()
        },
        scoring: ScoringArgs {
            gap_open: 6.0,
            gap_extend: 1.0,
            mismatch_penalty: 4.0,
            ..Default::default()
        },
        ..Default::default()
    }
}

fn frag_op_strategy() -> impl Strategy<Value = FragOp> {
    prop_oneof![Just(FragOp::AddTiedPerfect), Just(FragOp::AddCleanWin)]
}

struct PropVariant {
    pos: usize,
    p: f64,
}

impl Variant for PropVariant {
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
        self.p
    }
}

#[derive(Debug, Clone)]
enum FragOp {
    /// Add a new independent (tied-perfect) fragment: contributes 1 to
    /// ambiguous on both streams.
    AddTiedPerfect,
    /// Add a fragment where stream 0 wins cleanly.
    AddCleanWin,
}

/// Build a decisive winner/loser MD pair: `n_mismatch` mismatches, each
/// separated by at least one match base, with no zero-length digit run
/// anywhere in the MD string.
///
/// NOTE: `next_md_base` (src/alignment/ops.rs) advances `read_pos` and
/// returns `BaseOp::Match` unconditionally whenever it parses a digit
/// run -- including a run of `"0"` -- so a `"0"` token anywhere (leading,
/// trailing, or between two adjacent mismatch letters, e.g. "5A0C10")
/// injects a phantom match base with no corresponding read base skipped.
/// This generator avoids emitting any zero-length digit run so the test
/// exercises cross-algorithm agreement, not that parser quirk.
fn decisive_loser_md(len: usize, n_mismatch: usize) -> String {
    debug_assert!(len >= 2 * n_mismatch + 1);
    let mut md = "1A".repeat(n_mismatch);
    md.push_str(&(len - 2 * n_mismatch).to_string());
    md
}

/// A pair of fixtures with a match-count delta comfortably above the
/// default Phred-10 ambiguous threshold, so all three backends' Tier-2.5
/// / Tier-3 floors resolve the same winner (see module doc above p9).
#[derive(Debug, Clone)]
struct DecisivePair {
    len: usize,
    n_mismatch: usize,
    winner_on_stream0: bool,
}

fn decisive_pair_strategy() -> impl Strategy<Value = DecisivePair> {
    (3usize..=10, any::<bool>()).prop_flat_map(|(n_mismatch, winner_on_stream0)| {
        // len >= 2*n_mismatch + 1 keeps every MD digit run >= 1 (see decisive_loser_md).
        (2 * n_mismatch + 1..=2 * n_mismatch + 40).prop_map(move |len| DecisivePair {
            len,
            n_mismatch,
            winner_on_stream0,
        })
    })
}

impl DecisivePair {
    fn records(&self) -> Result<(RecordBuf, RecordBuf), Error> {
        let q = vec![30u8; self.len];
        let winner = create_record(
            b"R1",
            &format!("{}M", self.len),
            &[],
            &q,
            &self.len.to_string(),
            false,
        )?;
        let loser = create_record(
            b"R1",
            &format!("{}M", self.len),
            &[],
            &q,
            &decisive_loser_md(self.len, self.n_mismatch),
            false,
        )?;
        if self.winner_on_stream0 {
            Ok((winner, loser))
        } else {
            Ok((loser, winner))
        }
    }
}

fn two_stream_run_config() -> RunConfig {
    default_run_config(false, true, false)
}

struct StreamCountCase {
    label: &'static str,
    n: usize,
}

// ===========================================================================
// Table-driven edge cases for cross-algorithm / cross-parallelism agreement
// ===========================================================================
//
// Complements p9/p10: named rows for boundary conditions that random
// generation covers only by luck (exact ties, exact threshold cutoffs,
// minimal lengths). Failures report the row label, not a shrunk seed.

struct AgreementCase {
    label: &'static str,
    len: usize,
    n_mismatch: usize,
    winner_on_stream0: bool,
    /// If `Some(false)`, this case is a KNOWN boundary where Tier-2.5's
    /// floor may legitimately diverge from Tier-3 across backends (see
    /// `min_delta_for_early_decision`); we assert only that LineByLine's
    /// own sequential/parallel outputs agree, not full 3-backend equality.
    expect_full_backend_agreement: bool,
}

fn agreement_cases() -> Vec<AgreementCase> {
    vec![
        AgreementCase {
            label: "exact_tie_zero_mismatch_delta",
            len: 20,
            n_mismatch: 0,
            winner_on_stream0: true,
            expect_full_backend_agreement: true, // both perfect -> all backends ambiguous
        },
        AgreementCase {
            label: "minimal_length_single_mismatch_boundary",
            // len=2*1+1=3 is the smallest legal len for n_mismatch=1 per
            // decisive_loser_md's debug_assert.
            len: 3,
            n_mismatch: 1,
            winner_on_stream0: true,
            expect_full_backend_agreement: false, // delta too small vs Phred-10 floor
        },
        AgreementCase {
            label: "decisive_winner_on_stream0_large_delta",
            len: 50,
            n_mismatch: 8,
            winner_on_stream0: true,
            expect_full_backend_agreement: true,
        },
        AgreementCase {
            label: "decisive_winner_on_stream1_large_delta",
            len: 50,
            n_mismatch: 8,
            winner_on_stream0: false,
            expect_full_backend_agreement: true,
        },
        AgreementCase {
            label: "delta_just_above_phred10_threshold",
            // min_delta_for_early_decision at Phred-10 default; pick a
            // mismatch count comfortably (not marginally) above it so this
            // row is stable rather than flaky against penalty tuning --
            // see p9's module doc on why near-tie cases are excluded from
            // the "must fully agree" bucket.
            len: 30,
            n_mismatch: 4,
            winner_on_stream0: true,
            expect_full_backend_agreement: true,
        },
        AgreementCase {
            label: "long_read_many_mismatches",
            len: 150,
            n_mismatch: 20,
            winner_on_stream0: false,
            expect_full_backend_agreement: true,
        },
        AgreementCase {
            label: "short_read_max_mismatch_density",
            // len must satisfy len >= 2*n_mismatch+1; this is the tightest
            // packing (every other base a mismatch) at a small scale.
            len: 9,
            n_mismatch: 4,
            winner_on_stream0: true,
            expect_full_backend_agreement: false, // delta likely below floor at this scale
        },
    ]
}

// ===========================================================================
// P1 — Determinism: identical input -> identical routing counters
// ===========================================================================

proptest! {
    #[test]
    fn p1_deterministic_routing_for_identical_input(
        f0 in md_cig_strategy(60),
        f1 in md_cig_strategy(60),
        tag in any::<bool>(),
    ) {
        let cfg = default_run_config(tag, true, false);
        let rec0 = f0.record(b"R1", 30, false)?;
        let rec1 = f1.record(b"R1", 30, false)?;

        let run_once = |cfg: &RunConfig| -> Result<SmallVec<[u64; 8]>, Error> {
            let aln = n_stream_aln(vec![vec![rec0.clone()], vec![rec1.clone()]]);
            let mut lbl: LineByLine<RecordBuf> = LineByLine::new(cfg, aln, vec![], vec![])?;
            lbl.process_sequential(cfg)?;
            Ok(lbl.routing_counters)
        };

        let a = run_once(&cfg)?;
        let b = run_once(&cfg)?;
        prop_assert_eq!(a, b, "identical input must yield identical routing counters");
    }
}

// ===========================================================================
// P2 — Up to MAX_STREAMS (32) streams: no read is silently dropped
// ===========================================================================

proptest! {
    #![proptest_config(ProptestConfig { cases: 32, .. ProptestConfig::default() })]
    #[test]
    fn p2_all_streams_conserve_the_fragment(
        n in 2usize..=MAX_STREAMS,
        fixtures in prop::collection::vec(md_cig_strategy(40), 2..=MAX_STREAMS),
    ) {
        let n = n.min(fixtures.len());
        prop_assume!(n >= 2);

        let per_stream: Vec<Vec<RecordBuf>> = fixtures
            .iter()
            .take(n)
            .map(|f| vec![f.record(b"R1", 30, false).unwrap()])
            .collect();

        let cfg = default_run_config(false, true, false);
        let aln = n_stream_aln(per_stream);
        let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])?;
        lbl.process_sequential(&cfg)?;

        // Structural invariant, independent of *which* tier wins:
        // every stream received exactly one record, so summed over all
        // buckets (discard/out/ambiguous/chimeric) the total must equal n.
        let total: u64 = lbl.routing_counters.iter().sum();
        prop_assert_eq!(total, n as u64,
            "n={} streams each contributed one record; total routed must equal n, got {}", n, total);
    }
}

#[test]
fn table_stream_count_boundaries_conserve_fragment() {
    let cases = vec![
        StreamCountCase {
            label: "minimum_two_streams",
            n: 2,
        },
        StreamCountCase {
            label: "three_streams",
            n: 3,
        },
        StreamCountCase {
            label: "max_streams_32",
            n: MAX_STREAMS,
        },
    ];

    crate::tests::common::run_collecting(
        &cases,
        |c| c.label.to_string(),
        |c| {
            let cfg = default_run_config(false, true, false);
            let per_stream: Vec<Vec<RecordBuf>> = (0..c.n)
                .map(|i| {
                    // Vary MD slightly per stream so it's not a degenerate
                    // all-identical-perfect fixture at every n.
                    let md = if i == 0 {
                        "20".to_string()
                    } else {
                        format!("{}A{}", i.min(19), 19 - i.min(19))
                    };
                    create_record(b"R1", "20M", &[], &[30u8; 20], &md, false)
                        .map(|r| vec![r])
                        .map_err(|e| format!("create_record failed for stream {i}: {e:?}"))
                })
                .collect::<Result<Vec<Vec<_>>, String>>()?;

            let aln = n_stream_aln(per_stream);
            let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])
                .map_err(|e| format!("LineByLine::new failed: {e:?}"))?;
            lbl.process_sequential(&cfg)
                .map_err(|e| format!("process_sequential failed: {e:?}"))?;

            let total: u64 = lbl.routing_counters.iter().sum();
            if total != c.n as u64 {
                Err(format!("expected total={}, got {}", c.n, total))
            } else {
                Ok(())
            }
        },
    );
}

// ===========================================================================
// P3 — Confidence-tag / decision-outcome consistency
// ===========================================================================
//
// `apply_decision_tag` is the single place that turns a `Decision` into an
// XF/XR aux tag. The invariant: a tag is written iff the decision is
// PhredConfidence/VariantRescued, and the tag identity (XF vs XR) matches
// the decision variant exactly. Ambiguous/First/Last decisions must never
// write a tag (those states are conveyed by output stream, not aux tag).

proptest! {
    #[test]
    fn p3_decision_tag_matches_decision_variant(v in 0u8..=255u8, rescued in any::<bool>()) {
        use crate::filter_algorithm::line_by_line::ordering::{Decision, apply_decision_tag};
        use noodles::sam::alignment::record::data::field::Tag;

        let decision = if rescued {
            Decision::VariantRescued(v)
        } else {
            Decision::PhredConfidence(v)
        };

        let mut rec = RecordBuf::default();
        apply_decision_tag(&mut rec, Some(&decision));

        let xf = rec.data().get(&Tag::new(b'X', b'F'));
        let xr = rec.data().get(&Tag::new(b'X', b'R'));

        if rescued {
            prop_assert!(xr.is_some() && xf.is_none(), "VariantRescued must set XR only");
        } else {
            prop_assert!(xf.is_some() && xr.is_none(), "PhredConfidence must set XF only");
        }

        // Ambiguous/None must never write either tag.
        let mut rec2 = RecordBuf::default();
        apply_decision_tag(&mut rec2, Some(&Decision::Ambiguous));
        prop_assert!(rec2.data().get(&Tag::new(b'X', b'F')).is_none());
        prop_assert!(rec2.data().get(&Tag::new(b'X', b'R')).is_none());
    }
}

// ===========================================================================
// P4 — Sequential vs. parallel scoring equivalence
// ===========================================================================
//
// Output *order* is documented as nondeterministic under score_threads>1,
// but the aggregate routing decision per fragment must not depend on
// whether scoring happened on the IO thread or a worker thread.

proptest! {
    #![proptest_config(ProptestConfig { cases: 24, .. ProptestConfig::default() })]
    #[test]
    fn p4_sequential_and_parallel_agree(
        frags in prop::collection::vec((md_cig_strategy(50), md_cig_strategy(50)), 1..=6),
        threads in 2usize..=4,
    ) {
        let mut s0 = Vec::new();
        let mut s1 = Vec::new();
        for (i, (fa, fb)) in frags.iter().enumerate() {
            let name = format!("R{i}");
            s0.push(fa.record(name.as_bytes(), 30, false)?);
            s1.push(fb.record(name.as_bytes(), 30, false)?);
        }

        let cfg = default_run_config(true, true, false);

        let seq_counters = {
            let aln = n_stream_aln(vec![s0.clone(), s1.clone()]);
            let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])?;
            lbl.process_sequential(&cfg)?;
            lbl.routing_counters
        };
        let par_counters = {
            let aln = n_stream_aln(vec![s0, s1]);
            let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])?;
            lbl.process_parallel(&cfg, threads)?;
            lbl.routing_counters
        };

        prop_assert_eq!(seq_counters, par_counters,
            "sequential and parallel scoring must agree on aggregate routing");
    }
}

// ===========================================================================
// P5 — Ambiguous stability under non-informative noise
// ===========================================================================
//
// Two identical perfect records tie -> ambiguous for both streams. Adding
// extra, *unambiguously resolved* fragments (perfect vs clearly imperfect,
// so they can never land in the ambiguous bucket themselves) must not
// change the ambiguous count contributed by the original tied pair.

proptest! {
    #[test]
    fn p5_ambiguous_pair_unaffected_by_unrelated_noise(
        noise_count in 0usize..8,
    ) {
        let cfg = default_run_config(false, true, false);

        let tied0 = create_record(b"TIED", "10M", &[], &[30u8; 10], "10", false)?;
        let tied1 = create_record(b"TIED", "10M", &[], &[30u8; 10], "10", false)?;

        let mut s0 = vec![tied0];
        let mut s1 = vec![tied1];
        for i in 0..noise_count {
            let name = format!("NOISE{i}");
            // Stream 0 perfect, stream 1 imperfect -> resolves cleanly, never ambiguous.
            s0.push(create_record(name.as_bytes(), "10M", &[], &[30u8; 10], "10", false)?);
            s1.push(create_record(name.as_bytes(), "10M", &[], &[30u8; 10], "5A4", false)?);
        }
        // Streams must stay name-sorted per the contract; both lists share
        // the same qname order (TIED first, then NOISE0..N in lockstep).

        let aln = n_stream_aln(vec![s0, s1]);
        let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])?;
        lbl.process_sequential(&cfg)?;

        // ambiguous bucket is index [nr*4 + 2]
        let ambig0 = lbl.routing_counters[2];
        let ambig1 = lbl.routing_counters[6];
        prop_assert_eq!(ambig0, 1, "exactly the TIED fragment should be ambiguous on stream 0");
        prop_assert_eq!(ambig1, 1, "exactly the TIED fragment should be ambiguous on stream 1");
    }
}

// ===========================================================================
// P6 — Chimeric-pair contract: both streams emit, XC:Z present
// ===========================================================================

proptest! {
    #![proptest_config(ProptestConfig { cases: 20, .. ProptestConfig::default() })]
    #[test]
    fn p6_chimeric_pair_emits_both_streams(
        len in 40usize..=200,
        split_frac in 0.20f64..0.80, // keeps both arms well above min_mapped_frac(0.10)
    ) {
        let s = ((len as f64) * split_frac).round() as usize;
        prop_assume!(s >= 4 && len - s >= 4);

        let q = vec![30u8; len];
        // Stream A: matches read[0..s), soft-clips the rest.
        let rec_a = create_record(b"R1", &format!("{s}M{}S", len - s), &vec![b'A'; len], &q, &s.to_string(), false)?;
        // Stream B: soft-clips [0..s), matches the rest.
        let rec_b = create_record(b"R1", &format!("{s}S{}M", len - s), &vec![b'A'; len], &q, &(len - s).to_string(), false)?;

        let mut best: FragmentBuffer<RecordBuf> = smallvec![
            FragmentState::from_record(rec_a, 0, false)?,
            FragmentState::from_record(rec_b, 1, false)?,
        ];
        let decision = detect_chimeric_event(&best, &[[0, 1]], &ChimericThresholds::default());

        prop_assert!(
            matches!(decision, ChimericDecision::Chimeric { stream_a: 0, stream_b: 1, .. }),
            "complementary read-split with 0 overlap / full union must be detected as chimeric \
             (len={}, split={}), got {:?}", len, s, decision
        );

        // Full-pipeline contract: both outputs actually receive the record,
        // tagged XC:Z, and the chimeric counter increments on both sides.
        let cfg = default_run_config(false, true, false);
        let rec_a2 = create_record(b"R1", &format!("{s}M{}S", len - s), &vec![b'A'; len], &q, &s.to_string(), false)?;
        let rec_b2 = create_record(b"R1", &format!("{s}S{}M", len - s), &vec![b'A'; len], &q, &(len - s).to_string(), false)?;
        let aln = n_stream_aln(vec![vec![rec_a2], vec![rec_b2]]);
        let mut lbl: LineByLine<RecordBuf> = LineByLine::new(
            &cfg, aln, vec!["a".into(), "b".into()], vec![[0, 1]],
        )?;
        lbl.process_sequential(&cfg)?;
        prop_assert_eq!(lbl.routing_counters[3], 1, "stream 0 chimeric counter");
        prop_assert_eq!(lbl.routing_counters[7], 1, "stream 1 chimeric counter");

        drop(best.drain(..));
    }
}

// ===========================================================================
// P7 — Variant rescue never contributes below the evidence gate
// ===========================================================================
//
// wis_max_rescue_delta must only ever be moved by variants with
// delta() > 0.0 AND p_variant() >= min_p_variant. This is the pure-function
// core of `--min-rescue-p-variant`; test it directly rather than through
// the full NW pipeline so failures shrink to a single (delta, p) pair.

proptest! {
    #[test]
    fn p7_rescue_respects_min_p_variant_gate(
        delta in -5.0f64..5.0,
        p_variant in 0.0f64..1.0,
        min_p_variant in 0.0f64..1.0,
    ) {
        let v: &'static PropVariant = Box::leak(Box::new(PropVariant { pos: 0, p: p_variant }));
        let mut ev = Eval::new();
        ev.set_variant(v);
        ev.update(0.0, delta); // delta() == alt_score - incurred == delta

        let mut dvnt = smallvec![smallvec![ev]];
        let mut dp = smallvec![];
        let got = wis_max_rescue_delta(&mut dvnt, &mut dp, min_p_variant);

        let should_contribute = delta > 0.0 && p_variant >= min_p_variant;
        if should_contribute {
            prop_assert!((got - delta).abs() < 1e-9,
                "eligible variant (delta={delta}, p={p_variant}, gate={min_p_variant}) must contribute exactly its delta, got {got}",
                delta=delta, p_variant=p_variant, min_p_variant=min_p_variant, got=got);
        } else {
            prop_assert_eq!(got, 0.0,
                "ineligible variant (delta={delta}, p={p_variant}, gate={min_p_variant}) must contribute nothing, got {got}",
                delta=delta, p_variant=p_variant, min_p_variant=min_p_variant, got=got);
        }
    }
}

// ===========================================================================
// P8 — Cosmetic flags don't change classification
// ===========================================================================
//
// add_decision_tag / quiet / no_program_line are output-formatting concerns
// only. Toggling any of them must not change discard/out/ambiguous/chimeric
// counts for identical input.

proptest! {
    #[test]
    fn p8_cosmetic_flags_do_not_affect_routing(
        f0 in md_cig_strategy(50),
        f1 in md_cig_strategy(50),
        tag in any::<bool>(),
        quiet in any::<bool>(),
        no_pg in any::<bool>(),
    ) {
        let rec0 = f0.record(b"R1", 30, false)?;
        let rec1 = f1.record(b"R1", 30, false)?;

        let baseline_cfg = default_run_config(false, true, false);
        let variant_cfg = default_run_config(tag, quiet, no_pg);

        let run = |cfg: &RunConfig, r0: RecordBuf, r1: RecordBuf| -> Result<SmallVec<[u64; 8]>, Error> {
            let aln = n_stream_aln(vec![vec![r0], vec![r1]]);
            let mut lbl: LineByLine<RecordBuf> = LineByLine::new(cfg, aln, vec![], vec![])?;
            lbl.process_sequential(cfg)?;
            Ok(lbl.routing_counters)
        };

        let a = run(&baseline_cfg, rec0.clone(), rec1.clone())?;
        let b = run(&variant_cfg, rec0, rec1)?;
        prop_assert_eq!(a, b,
            "toggling add_decision_tag/quiet/no_program_line must not change routing counters");
    }
}

proptest! {
    #[test]
    fn stateful_ambiguous_and_win_counts_track_ops(
        ops in prop::collection::vec(frag_op_strategy(), 0..12),
    ) {
        let cfg = default_run_config(false, true, false);
        let mut s0 = Vec::new();
        let mut s1 = Vec::new();
        let mut expect_ambig = 0u64;
        let mut expect_win0 = 0u64;

        for (i, op) in ops.iter().enumerate() {
            let name = format!("F{i}");
            match op {
                FragOp::AddTiedPerfect => {
                    s0.push(create_record(name.as_bytes(), "10M", &[], &[30u8; 10], "10", false)?);
                    s1.push(create_record(name.as_bytes(), "10M", &[], &[30u8; 10], "10", false)?);
                    expect_ambig += 1;
                }
                FragOp::AddCleanWin => {
                    s0.push(create_record(name.as_bytes(), "10M", &[], &[30u8; 10], "10", false)?);
                    s1.push(create_record(name.as_bytes(), "10M", &[], &[30u8; 10], "5A4", false)?);
                    expect_win0 += 1;
                }
            }

            // Re-run the whole pipeline on the prefix built so far and check
            // the model still matches — this is what makes it "stateful":
            // the invariant must hold after every prefix, not just at the end.
            let aln = n_stream_aln(vec![s0.clone(), s1.clone()]);
            let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])?;
            lbl.process_sequential(&cfg)?;
            prop_assert_eq!(lbl.routing_counters[2], expect_ambig, "ambiguous count after {} ops", i + 1);
            prop_assert_eq!(lbl.routing_counters[1], expect_win0, "stream-0 win count after {} ops", i + 1);
        }
    }
}

// ===========================================================================
// P9 — Namesorted-compatible input: all three algorithms agree
// ===========================================================================
//
// LineByLine (namesorted), HashLookup, and CollatedMatcher all claim to
// implement the same tournament semantics on 2-stream input, differing
// only in I/O strategy (streaming vs. two-pass seek vs. hash-buffered).
// For a decisively-scored fragment (match-count delta far above the
// default ambiguous threshold, so no backend's Tier-2.5 floor disagrees
// with what Tier-3 NW would resolve), the aggregate routing counters
// (stride-4 per stream: discard/out/ambiguous/chimeric) must be identical
// across all three.
proptest! {
    #![proptest_config(ProptestConfig { cases: 40, .. ProptestConfig::default() })]
    #[test]
    fn p9_all_three_algorithms_agree_on_namesorted_input(pair in decisive_pair_strategy()) {
        let (rec0, rec1) = pair.records()?;
        let cfg = two_stream_run_config();

        let lbl_counters = {
            let aln = n_stream_aln(vec![vec![rec0.clone()], vec![rec1.clone()]]);
            let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])?;
            lbl.process_sequential(&cfg)?;
            lbl.routing_counters
        };

        let hl_counters = {
            let aln = crate::tests::common::two_stream_aln(vec![rec0.clone()], vec![rec1.clone()]);
            let mut hl = HashLookup::<RecordBuf>::new(&cfg, aln, [None, None], [None, None])?;
            hl.process(&cfg)?;
            hl.routing_counters
        };

        let coll_counters = {
            let aln = crate::tests::common::two_stream_aln(vec![rec0], vec![rec1]);
            let mut coll = CollatedMatcher::<RecordBuf>::new(&cfg, aln, [None, None], [None, None])?;
            coll.process(&cfg)?;
            coll.routing_counters
        };

        prop_assert_eq!(
            &lbl_counters, &hl_counters,
            "LineByLine vs HashLookup diverged for len={} n_mismatch={} winner_on_stream0={}",
            pair.len, pair.n_mismatch, pair.winner_on_stream0
        );
        prop_assert_eq!(
            &lbl_counters, &coll_counters,
            "LineByLine vs CollatedMatcher diverged for len={} n_mismatch={} winner_on_stream0={}",
            pair.len, pair.n_mismatch, pair.winner_on_stream0
        );
    }
}

// ===========================================================================
// P10 — Cross-algorithm agreement extended to sequential vs parallel LineByLine
// ===========================================================================
//
// Composes with p4: the same decisive fixture must agree across all of
// {sequential, parallel, HashLookup, Collated} simultaneously -- a
// regression here pinpoints whether a divergence is backend-specific
// (only HashLookup or only Collated disagrees with LineByLine) or a
// parallelism-specific regression (only process_parallel disagrees),
// which the pairwise p4/p9 tests can each miss if run in isolation but
// this one catches jointly.
proptest! {
    #![proptest_config(ProptestConfig { cases: 24, .. ProptestConfig::default() })]
    #[test]
    fn p10_sequential_parallel_and_all_backends_agree(
        pair in decisive_pair_strategy(),
        threads in 2usize..=4,
    ) {
        let cfg = two_stream_run_config();

        let (r0, r1) = pair.records()?;
        let seq = {
            let aln = n_stream_aln(vec![vec![r0.clone()], vec![r1.clone()]]);
            let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])?;
            lbl.process_sequential(&cfg)?;
            lbl.routing_counters
        };
        let par = {
            let aln = n_stream_aln(vec![vec![r0.clone()], vec![r1.clone()]]);
            let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])?;
            lbl.process_parallel(&cfg, threads)?;
            lbl.routing_counters
        };
        let hl = {
            let aln = crate::tests::common::two_stream_aln(vec![r0.clone()], vec![r1.clone()]);
            let mut hl = HashLookup::<RecordBuf>::new(&cfg, aln, [None, None], [None, None])?;
            hl.process(&cfg)?;
            hl.routing_counters
        };
        let coll = {
            let aln = crate::tests::common::two_stream_aln(vec![r0], vec![r1]);
            let mut coll = CollatedMatcher::<RecordBuf>::new(&cfg, aln, [None, None], [None, None])?;
            coll.process(&cfg)?;
            coll.routing_counters
        };

        prop_assert_eq!(&seq, &par, "sequential vs parallel diverged");
        prop_assert_eq!(&seq, &hl, "sequential vs HashLookup diverged");
        prop_assert_eq!(&seq, &coll, "sequential vs CollatedMatcher diverged");
    }
}

#[test]
fn table_cross_algorithm_agreement() {
    let cases = agreement_cases();
    crate::tests::common::run_collecting(
        &cases,
        |c| c.label.to_string(),
        |c| {
            let pair = DecisivePair {
                len: c.len,
                n_mismatch: c.n_mismatch,
                winner_on_stream0: c.winner_on_stream0,
            };
            let cfg = two_stream_run_config();

            let (r0, r1) = pair
                .records()
                .map_err(|e| format!("record construction failed: {e:?}"))?;

            let lbl_counters = {
                let aln = n_stream_aln(vec![vec![r0.clone()], vec![r1.clone()]]);
                let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])
                    .map_err(|e| format!("LineByLine::new failed: {e:?}"))?;
                lbl.process_sequential(&cfg)
                    .map_err(|e| format!("process_sequential failed: {e:?}"))?;
                lbl.routing_counters
            };
            let par_counters = {
                let aln = n_stream_aln(vec![vec![r0.clone()], vec![r1.clone()]]);
                let mut lbl: LineByLine<RecordBuf> = LineByLine::new(&cfg, aln, vec![], vec![])
                    .map_err(|e| format!("LineByLine::new (parallel) failed: {e:?}"))?;
                lbl.process_parallel(&cfg, 3)
                    .map_err(|e| format!("process_parallel failed: {e:?}"))?;
                lbl.routing_counters
            };
            if lbl_counters != par_counters {
                return Err(format!(
                    "sequential {:?} != parallel {:?}",
                    lbl_counters, par_counters
                ));
            }

            if !c.expect_full_backend_agreement {
                // Known boundary case: only sequential/parallel LineByLine
                // agreement is asserted above; cross-backend Tier-2.5 floor
                // behavior is allowed to differ here. Nothing further to check.
                return Ok(());
            }

            let hl_counters = {
                let aln = crate::tests::common::two_stream_aln(vec![r0.clone()], vec![r1.clone()]);
                let mut hl = HashLookup::<RecordBuf>::new(&cfg, aln, [None, None], [None, None])
                    .map_err(|e| format!("HashLookup::new failed: {e:?}"))?;
                hl.process(&cfg)
                    .map_err(|e| format!("HashLookup::process failed: {e:?}"))?;
                hl.routing_counters
            };
            if lbl_counters != hl_counters {
                return Err(format!(
                    "LineByLine {:?} != HashLookup {:?}",
                    lbl_counters, hl_counters
                ));
            }

            let coll_counters = {
                let aln = crate::tests::common::two_stream_aln(vec![r0], vec![r1]);
                let mut coll =
                    CollatedMatcher::<RecordBuf>::new(&cfg, aln, [None, None], [None, None])
                        .map_err(|e| format!("CollatedMatcher::new failed: {e:?}"))?;
                coll.process(&cfg)
                    .map_err(|e| format!("CollatedMatcher::process failed: {e:?}"))?;
                coll.routing_counters
            };
            if lbl_counters != coll_counters {
                return Err(format!(
                    "LineByLine {:?} != CollatedMatcher {:?}",
                    lbl_counters, coll_counters
                ));
            }

            Ok(())
        },
    );
}
