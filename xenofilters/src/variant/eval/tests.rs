// Replacement for src/variant/eval/tests.rs
//
// Table-driven tests that reproduce the original per-function assertions.
// Each case includes a human-readable test name (original test function name)
// which is printed if the case fails. All failures are collected and reported
// together at the end so one run shows all misses.

use super::*;
use crate::variant::Variant;

// A small mock Variant used in the position/allele tests.
#[derive(Clone)]
struct MockVariant {
    pos: usize,
    ref_allele: Vec<u8>,
    alt_allele: Vec<u8>,
    p_variant: f64,
}
impl MockVariant {
    fn boxed(pos: usize, ref_a: &[u8], alt_a: &[u8], p: f64) -> &'static Self {
        Box::leak(Box::new(MockVariant {
            pos,
            ref_allele: ref_a.to_vec(),
            alt_allele: alt_a.to_vec(),
            p_variant: p,
        }))
    }
}
impl Variant for MockVariant {
    fn pos(&self) -> usize {
        self.pos
    }
    fn ref_allele(&self) -> &[u8] {
        &self.ref_allele
    }
    fn alt_allele(&self) -> &[u8] {
        &self.alt_allele
    }
    fn p_variant(&self) -> f64 {
        self.p_variant
    }
}

// ---------------------------
// Group A: accumulation tests
// ---------------------------
#[test]
fn table_eval_accumulation_collect_misses() {
    struct AccCase {
        name: &'static str,
        updates: &'static [(f64, f64)],
        want_incurred: f64,
        want_alt: f64,
        want_delta: f64,
    }

    let cases = &[
        // corresponds to `test_update_accumulates_scores`
        AccCase {
            name: "test_update_accumulates_scores",
            updates: &[(1.0, 2.0), (3.0, 4.0)],
            want_incurred: 4.0,
            want_alt: 6.0,
            want_delta: 2.0,
        },
        // corresponds to `delta_is_alt_minus_incurred`
        AccCase {
            name: "delta_is_alt_minus_incurred",
            updates: &[(2.5, 10.0), (1.0, 4.0)],
            // incurred = 2.5 + 1.0 = 3.5 ; alt = 10.0 + 4.0 = 14.0 ; delta = 10.5
            want_incurred: 3.5,
            want_alt: 14.0,
            want_delta: 10.5,
        },
        // corresponds to `eval_update_accumulates_values`
        AccCase {
            name: "eval_update_accumulates_values",
            updates: &[(2.0, 10.0), (3.0, 7.0)],
            // incurred = 5.0 ; alt = 17.0 ; delta = 12.0
            want_incurred: 5.0,
            want_alt: 17.0,
            want_delta: 12.0,
        },
        // baseline: new Eval zero initialized
        AccCase {
            name: "eval_new_is_zero_initialized",
            updates: &[],
            want_incurred: 0.0,
            want_alt: 0.0,
            want_delta: 0.0,
        },
    ];

    crate::tests::common::run_collecting(
        cases,
        |c| c.name.to_string(),
        |c| {
            let mut ev = Eval::new();
            for &(inc, alt) in c.updates {
                ev.update(inc, alt);
            }
            if (ev.incurred - c.want_incurred).abs() > 1e-12 {
                return Err(format!(
                    "incurred mismatch\n  expected incurred = {}\n  got      incurred = {}",
                    c.want_incurred, ev.incurred
                ));
            }
            if (ev.alt_score - c.want_alt).abs() > 1e-12 {
                return Err(format!(
                    "alt_score mismatch\n  expected alt_score = {}\n  got      alt_score = {}",
                    c.want_alt, ev.alt_score
                ));
            }
            if (ev.delta() - c.want_delta).abs() > 1e-12 {
                return Err(format!(
                    "delta mismatch\n  expected delta = {}\n  got      delta = {}",
                    c.want_delta,
                    ev.delta()
                ));
            }
            Ok(())
        },
    );
}

// -------------------------------------------------------
// Group B: variant position/allele endpoints and boundaries
// -------------------------------------------------------
#[test]
fn table_variant_position_and_endpoints_collect_misses() {
    struct PosCase {
        name: &'static str,
        pos: usize, // 1-based per Variant trait docs in repo; tests used values accordingly
        ref_a: &'static [u8],
        alt_a: &'static [u8],
        want_start: usize,
        want_ref_end: usize,
        want_alt_end: usize,
        want_end: usize,
    }

    let cases = &[
        // start_matches_variant_position
        PosCase {
            name: "start_matches_variant_position",
            pos: 123,
            ref_a: b"A",
            alt_a: b"G",
            want_start: 123,
            want_ref_end: 124,
            want_alt_end: 124,
            want_end: 124,
        },
        // ref_end_uses_reference_length
        PosCase {
            name: "ref_end_uses_reference_length",
            pos: 100,
            ref_a: b"ACTG",
            alt_a: b"A",
            want_start: 100,
            want_ref_end: 104,
            want_alt_end: 101,
            want_end: 104,
        },
        // alt_end_uses_alternate_length
        PosCase {
            name: "alt_end_uses_alternate_length",
            pos: 100,
            ref_a: b"A",
            alt_a: b"ACTG",
            want_start: 100,
            want_ref_end: 101,
            want_alt_end: 104,
            want_end: 104,
        },
        // insertion: end prefers longest allele span
        PosCase {
            name: "end_prefers_longest_allele_span_for_insertion",
            pos: 50,
            ref_a: b"A",
            alt_a: b"ACTGA",
            want_start: 50,
            want_ref_end: 51,
            want_alt_end: 55,
            want_end: 55,
        },
        // deletion: end prefers longest allele span
        PosCase {
            name: "end_prefers_longest_allele_span_for_deletion",
            pos: 50,
            ref_a: b"ACTGA",
            alt_a: b"A",
            want_start: 50,
            want_ref_end: 55,
            want_alt_end: 51,
            want_end: 55,
        },
    ];

    crate::tests::common::run_collecting(
        cases,
        |c| c.name.to_string(),
        |c| {
            let mv = MockVariant::boxed(c.pos, c.ref_a, c.alt_a, 0.5);
            let mut ev = Eval::new();
            ev.set_variant(mv as &dyn Variant);

            if ev.start() != c.want_start {
                return Err(format!(
                    "start mismatch\n  expected start = {}\n  got      start = {}",
                    c.want_start,
                    ev.start()
                ));
            }
            if ev.ref_end() != c.want_ref_end {
                return Err(format!(
                    "ref_end mismatch\n  expected ref_end = {}\n  got      ref_end = {}",
                    c.want_ref_end,
                    ev.ref_end()
                ));
            }
            if ev.alt_end() != c.want_alt_end {
                return Err(format!(
                    "alt_end mismatch\n  expected alt_end = {}\n  got      alt_end = {}",
                    c.want_alt_end,
                    ev.alt_end()
                ));
            }
            if ev.end() != c.want_end {
                return Err(format!(
                    "end() mismatch\n  expected end = {}\n  got      end = {}",
                    c.want_end,
                    ev.end()
                ));
            }
            Ok(())
        },
    );
}
