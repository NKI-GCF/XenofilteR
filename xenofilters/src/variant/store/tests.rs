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
    let mut store = Store::new();
    for (chr, v) in variants {
        store.insert(chr, v);
    }
    store.dedup();
    store
}

#[test]
fn test_overlap_basic_hit() {
    let store = make_store(vec![(1, fv(100, 1))]); // [100,101)
    assert_eq!(store.overlapping(1, 100, 101).len(), 1);
}

#[test]
fn test_overlap_exclusions() {
    let store = make_store(vec![(1, fv(100, 1))]); // start 100, ends at 101

    let cases = [
        (2, 100, 101), // wrong chromosome
        (1, 101, 110), // starts exactly at variant end
        (1, 90, 100),  // ends exactly at variant start
    ];

    for (q_ref, q_start, q_end) in cases {
        assert!(store.overlapping(q_ref, q_start, q_end).is_empty());
    }
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

#[test]
fn insert_preserves_stable_order_for_equal_positions() {
    let mut store = Store::<FakeVariant>::new();
    store.insert(1, fv(100, 1));
    store.insert(1, fv(100, 2)); // same pos, different ref_len -- must land after
    store.insert(1, fv(100, 3));
    let bucket =  &store.per_chr_data[&1];
    let lens: Vec<usize> = bucket.iter().map(|v| v.ref_allele().len()).collect();
    assert_eq!(
        lens,
        vec![1, 2, 3],
        "insertion order at equal pos must be preserved"
    );
}

#[test]
fn dedup_only_merges_byte_identical_adjacent_entries() {
    let mut store = Store::<FakeVariant>::new();
    store.insert(1, fv(100, 1));
    store.insert(1, fv(100, 1)); // exact duplicate
    store.insert(1, fv(100, 2)); // different ref_len -- must survive
    store.dedup();
    assert_eq!(store.per_chr_data[&1].len(), 2);
}
