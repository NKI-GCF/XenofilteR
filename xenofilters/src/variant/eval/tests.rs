use super::*;
use crate::variant::Variant;

struct MockVariant {
    pos: usize,
    ref_allele: &'static [u8],
    alt_allele: &'static [u8],
    p_variant: f64,
}

impl Variant for MockVariant {
    fn pos(&self) -> usize {
        self.pos
    }

    fn ref_allele(&self) -> &[u8] {
        self.ref_allele
    }

    fn alt_allele(&self) -> &[u8] {
        self.alt_allele
    }

    fn p_variant(&self) -> f64 {
        self.p_variant
    }
}

#[test]
fn delta_is_alt_minus_incurred() {
    let mut ev = Eval::new();

    ev.update(2.5, 10.0);
    ev.update(1.0, 4.0);

    assert!((ev.delta() - 10.5).abs() < 1e-10);
}

#[test]
fn start_matches_variant_position() {
    let v = MockVariant {
        pos: 123,
        ref_allele: b"A",
        alt_allele: b"G",
        p_variant: 0.5,
    };

    let mut ev = Eval::new();
    ev.set_variant(&v);

    assert_eq!(ev.start(), 123);
}

#[test]
fn ref_end_uses_reference_length() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"ACTG",
        alt_allele: b"A",
        p_variant: 0.5,
    };

    let mut ev = Eval::new();
    ev.set_variant(&v);

    assert_eq!(ev.ref_end(), 104);
}

#[test]
fn alt_end_uses_alternate_length() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"A",
        alt_allele: b"ACTG",
        p_variant: 0.5,
    };

    let mut ev = Eval::new();
    ev.set_variant(&v);

    assert_eq!(ev.alt_end(), 104);
}

#[test]
fn end_prefers_longest_allele_span_for_insertion() {
    let v = MockVariant {
        pos: 50,
        ref_allele: b"A",
        alt_allele: b"ACTGA",
        p_variant: 0.5,
    };

    let mut ev = Eval::new();
    ev.set_variant(&v);

    assert_eq!(ev.ref_end(), 51);
    assert_eq!(ev.alt_end(), 55);
    assert_eq!(ev.end(), 55);
}

#[test]
fn end_prefers_longest_allele_span_for_deletion() {
    let v = MockVariant {
        pos: 50,
        ref_allele: b"ACTGA",
        alt_allele: b"A",
        p_variant: 0.5,
    };

    let mut ev = Eval::new();
    ev.set_variant(&v);

    assert_eq!(ev.ref_end(), 55);
    assert_eq!(ev.alt_end(), 51);
    assert_eq!(ev.end(), 55);
}

#[test]
fn update_accumulates_scores() {
    let mut ev = Eval::new();

    ev.update(1.0, 2.0);
    ev.update(3.0, 4.0);

    assert!((ev.incurred - 4.0).abs() < 1e-10);
    assert!((ev.alt_score - 6.0).abs() < 1e-10);
    assert!((ev.delta() - 2.0).abs() < 1e-10);
}
#[test]
fn eval_new_is_zero_initialized() {
    let ev = Eval::new();

    assert_eq!(ev.delta(), 0.0);
}

#[test]
fn eval_update_accumulates_values() {
    let mut ev = Eval::new();

    ev.update(2.0, 10.0);
    ev.update(3.0, 7.0);

    assert!((ev.delta() - 12.0).abs() < 1e-12);
}

#[test]
fn eval_delta_is_alt_minus_incurred() {
    let mut ev = Eval::new();

    ev.update(4.5, 11.0);

    assert!((ev.delta() - 6.5).abs() < 1e-12);
}

#[test]
fn eval_reports_variant_position() {
    let v = MockVariant {
        pos: 1234,
        ref_allele: b"A",
        alt_allele: b"G",
        p_variant: 0.5,
    };

    let mut ev = Eval::new();
    ev.set_variant(&v);

    assert_eq!(ev.start(), 1234);
}

#[test]
fn eval_ref_end_for_snp() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"A",
        alt_allele: b"G",
        p_variant: 0.5,
    };

    let mut ev = Eval::new();
    ev.set_variant(&v);

    assert_eq!(ev.ref_end(), 101);
}

#[test]
fn eval_alt_end_for_snp() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"A",
        alt_allele: b"G",
        p_variant: 0.5,
    };

    let mut ev = Eval::new();
    ev.set_variant(&v);

    assert_eq!(ev.alt_end(), 101);
}

#[test]
fn eval_end_tracks_longer_alt_for_insertion() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"A",
        alt_allele: b"ATGC",
        p_variant: 0.5,
    };

    let mut ev = Eval::new();
    ev.set_variant(&v);

    assert_eq!(ev.ref_end(), 101);
    assert_eq!(ev.alt_end(), 104);
    assert_eq!(ev.end(), 104);
}

#[test]
fn eval_end_tracks_longer_ref_for_deletion() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"ATGC",
        alt_allele: b"A",
        p_variant: 0.5,
    };

    let mut ev = Eval::new();
    ev.set_variant(&v);

    assert_eq!(ev.ref_end(), 104);
    assert_eq!(ev.alt_end(), 101);
    assert_eq!(ev.end(), 104);
}

#[test]
fn variant_overlap_detects_internal_overlap() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"ATG",
        alt_allele: b"A",
        p_variant: 0.5,
    };

    assert!(v.overlaps(101, 102));
}

#[test]
fn variant_overlap_detects_left_boundary_overlap() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"ATG",
        alt_allele: b"A",
        p_variant: 0.5,
    };

    assert!(v.overlaps(99, 101));
}

#[test]
fn variant_overlap_detects_right_boundary_overlap() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"ATG",
        alt_allele: b"A",
        p_variant: 0.5,
    };

    assert!(v.overlaps(102, 105));
}

#[test]
fn variant_overlap_excludes_touching_left_edge() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"ATG",
        alt_allele: b"A",
        p_variant: 0.5,
    };

    assert!(!v.overlaps(90, 100));
}

#[test]
fn variant_overlap_excludes_touching_right_edge() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"ATG",
        alt_allele: b"A",
        p_variant: 0.5,
    };

    assert!(!v.overlaps(103, 110));
}

#[test]
fn variant_overlap_excludes_disjoint_regions() {
    let v = MockVariant {
        pos: 100,
        ref_allele: b"ATG",
        alt_allele: b"A",
        p_variant: 0.5,
    };

    assert!(!v.overlaps(1, 50));
    assert!(!v.overlaps(200, 300));
}

struct FakeVariant {
    pos: usize,
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
}
impl Variant for FakeVariant {
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
    }
}

#[test]
fn test_start_ref_end_alt_end_and_end() {
    let v = FakeVariant {
        pos: 100,
        ref_a: vec![b'A'; 3],
        alt_a: vec![b'A'; 5],
    };
    let mut eval = Eval::new();
    eval.set_variant(&v);
    assert_eq!(eval.start(), 100);
    assert_eq!(eval.ref_end(), 103);
    assert_eq!(eval.alt_end(), 105);
    assert_eq!(eval.end(), 105);
}

#[test]
fn test_update_accumulates_and_delta_is_alt_minus_incurred() {
    let v = FakeVariant {
        pos: 0,
        ref_a: vec![b'A'],
        alt_a: vec![b'A'],
    };
    let mut eval = Eval::new();
    eval.set_variant(&v);
    eval.update(1.0, 2.0);
    eval.update(0.5, 0.5);
    assert!((eval.delta() - 1.0).abs() < 1e-12); // (2.0+0.5) - (1.0+0.5)
}

#[test]
#[should_panic(expected = "VariantEval should always have a variant reference")]
fn test_vnt_panics_when_unset() {
    let _ = Eval::new().vnt();
}
