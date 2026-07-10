// src/variant/store_insert.rs
//
// Corrected Store<V>::insert signatures.
//
// The earlier indel_equiv_impls.rs called `v.ref_id()` which requires
// Variant::ref_id().  This file supersedes that approach with an explicit
// ref_id parameter so that the Variant trait change is isolated to adding
// one method, not restructuring call sites.
//
// Merge instructions:
//   1. Replace the `insert` and `insert_expanded` bodies in indel_equiv_impls.rs
//      with the versions below (or add this file and remove those bodies).
//   2. Ensure Store<V> has a `per_ref: HashMap<usize, Vec<V>, RandomState>` field.
//   3. Ensure Variant trait gains `fn ref_id(&self) -> usize`.

use ahash::RandomState;
use std::collections::HashMap;

use crate::variant::{store::Store, Variant};

// ---------------------------------------------------------------------------
// Store<V> insert methods
// ---------------------------------------------------------------------------

impl<V: Variant + Clone> Store<V> {
    /// Insert one variant into the per-chromosome sorted bucket.
    ///
    /// `ref_id` is the 0-based chromosome index from the reference header.
    /// Uses `partition_point` (binary search) to find the insertion index;
    /// O(log n + n) per call — acceptable at startup, never called on the
    /// hot path.
    pub(crate) fn insert(&mut self, ref_id: usize, v: V) {
        let bucket = self.per_ref.entry(ref_id).or_insert_with(Vec::new);
        let pos = v.pos();
        // Stable sort: existing entries with the same pos keep their order;
        // new entry inserted after them.
        let idx = bucket.partition_point(|e| e.pos() <= pos);
        bucket.insert(idx, v);
    }

    /// Insert a batch of pre-expanded variants.  All share the same ref_id.
    /// Delegates to `insert` per element; no separate finalization required
    /// because each insert maintains the sorted invariant.
    pub(crate) fn insert_expanded(&mut self, ref_id: usize, variants: Vec<V>) {
        // Reserve capacity for the whole batch upfront.
        let bucket = self.per_ref.entry(ref_id).or_insert_with(Vec::new);
        bucket.reserve(variants.len());
        drop(bucket); // release mutable borrow before calling insert in loop

        for v in variants {
            self.insert(ref_id, v);
        }
    }

    /// Remove byte-exact duplicates introduced when two VCF records in a
    /// repeat region expand to the same (pos, ref_a, alt_a) tuple.
    ///
    /// Called once after all inserts are complete.  The sorted invariant
    /// ensures adjacent equal elements are the candidates; `dedup_by` is
    /// O(n) after sorting.
    pub(crate) fn dedup(&mut self) {
        for bucket in self.per_ref.values_mut() {
            bucket.dedup_by(|a, b| {
                // b comes first in sorted order; dedup_by retains b.
                a.pos()         == b.pos()
                    && a.ref_allele() == b.ref_allele()
                    && a.alt_allele() == b.alt_allele()
            });
        }
    }
}
