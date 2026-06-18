use super::*;
use crate::variant::Variant;

struct FakeVariant {
    pos: usize,
    ref_a: Vec<u8>,
}
impl Variant for FakeVariant {
    fn pos(&self) -> usize {
        self.pos
    }
    fn ref_allele(&self) -> &[u8] {
        &self.ref_a
    }
    fn alt_allele(&self) -> &[u8] {
        b"A"
    }
    fn p_variant(&self) -> f64 {
        0.5
    }
}
fn fv(pos: usize, ref_len: usize) -> FakeVariant {
    FakeVariant {
        pos,
        ref_a: vec![b'A'; ref_len],
    }
}

fn make_store(variants: Vec<(usize, FakeVariant)>) -> Store<FakeVariant> {
    let mut per_chr: HashMap<usize, Vec<FakeVariant>> = HashMap::new();
    let mut max_variant_len = 1;
    for (chr, v) in variants {
        max_variant_len = max_variant_len.max(v.ref_allele().len());
        per_chr.entry(chr).or_default().push(v);
    }
    for vs in per_chr.values_mut() {
        vs.sort_by_key(|v| v.pos());
    }
    Store {
        per_chr,
        max_variant_len,
    }
}

#[test]
fn test_overlap_basic_hit() {
    let store = make_store(vec![(1, fv(100, 1))]); // [100,101)
    assert_eq!(store.overlapping(1, 100, 101).len(), 1);
}

#[test]
fn test_overlap_wrong_chromosome_no_hit() {
    let store = make_store(vec![(1, fv(100, 1))]);
    assert!(store.overlapping(2, 100, 101).is_empty());
}

#[test]
fn test_overlap_read_starts_exactly_at_variant_end_is_excluded() {
    let store = make_store(vec![(1, fv(100, 1))]); // ends at 101 (exclusive)
    assert!(store.overlapping(1, 101, 110).is_empty());
}

#[test]
fn test_overlap_read_ends_exactly_at_variant_start_is_excluded() {
    let store = make_store(vec![(1, fv(100, 1))]);
    assert!(store.overlapping(1, 90, 100).is_empty());
}

#[test]
fn test_overlap_partial_overlap_on_both_edges() {
    let store = make_store(vec![(1, fv(100, 5))]); // [100,105)
    assert_eq!(store.overlapping(1, 95, 101).len(), 1); // overlaps start
    assert_eq!(store.overlapping(1, 104, 110).len(), 1); // overlaps end
    assert_eq!(store.overlapping(1, 102, 103).len(), 1); // fully inside
}

#[test]
fn test_overlap_lower_bound_uses_max_variant_len_correctly() {
    let store = make_store(vec![(1, fv(50, 10))]); // [50,60), max_variant_len=10
    assert_eq!(store.overlapping(1, 55, 70).len(), 1);
    assert!(store.overlapping(1, 61, 70).is_empty());
}

#[test]
fn test_overlap_multi_variant_binary_search_selects_correct_one() {
    let store = make_store(vec![(1, fv(10, 1)), (1, fv(200, 1)), (1, fv(500, 1))]);
    let hits = store.overlapping(1, 195, 205);
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].pos(), 200);
}
