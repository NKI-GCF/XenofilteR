use super::*;
use crate::variant::Variant;

struct MockVariant {
    pos: usize,
    ref_allele: &'static [u8],
    alt_allele: &'static [u8],
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
        0.5
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
