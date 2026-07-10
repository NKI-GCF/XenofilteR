// src/variant/indel_equiv_impls.rs
//
// Concrete implementations of WithAlleles for Sample and Population,
// and the Store<V>::insert + insert_expanded methods.
//
// Merge notes:
//   - Add `#[derive(Clone)]` to Sample and Population if not present.
//   - Add `pub(crate) fn insert(&mut self, v: V)` to Store<V> if not present.
//   - Re-export `WithAlleles` from `crate::variant`.

use crate::variant::{
    indel_equiv::WithAlleles,
    population::Population,
    sample::Sample,
    store::Store,
    Variant,
};

// ---------------------------------------------------------------------------
// WithAlleles for Sample
// ---------------------------------------------------------------------------
//
// Sample stores: pos, ref_a, alt_a, genotype_quality, alt_haplotype,
// phase_set, gamete_mode (from phasing work).  All fields except the
// three allele fields are cloned verbatim.

impl WithAlleles for Sample {
    fn with_alleles(&self, pos_0based: usize, ref_a: &[u8], alt_a: &[u8]) -> Self {
        Self {
            pos:               pos_0based,
            ref_a:             ref_a.to_vec(),
            alt_a:             alt_a.to_vec(),
            // -- Scoring parameters propagated unchanged ------------------
            genotype_quality:  self.genotype_quality,
            alt_haplotype:     self.alt_haplotype,
            phase_set:         self.phase_set,
            gamete_mode:       self.gamete_mode,
        }
    }
}

// ---------------------------------------------------------------------------
// WithAlleles for Population
// ---------------------------------------------------------------------------

impl WithAlleles for Population {
    fn with_alleles(&self, pos_0based: usize, ref_a: &[u8], alt_a: &[u8]) -> Self {
        Self {
            pos:               pos_0based,
            ref_a:             ref_a.to_vec(),
            alt_a:             alt_a.to_vec(),
            allele_frequency:  self.allele_frequency,
        }
    }
}

// ---------------------------------------------------------------------------
// Store<V> - insert and insert_expanded
// ---------------------------------------------------------------------------
//
// The existing Store<V> is built in bulk from a VCF file.  We add two
// incremental methods used by the expander at startup.  The hot-path
// `overlapping_multi` is unchanged.
//
// Internal invariant: per_ref[ref_id] is sorted by Variant::pos() and
// secondarily by ref_allele length (longer first, for correct interval
// overlap).  `finalize()` must be called after all inserts are complete.

impl<V: Variant + Ord> Store<V> {
    /// Insert one variant into the store, maintaining sort order.
    ///
    /// Uses binary search + insert (O(log n + n) per call).  Acceptable
    /// at startup; not called on the hot path.
    pub(crate) fn insert(&mut self, v: V) {
        let ref_id = v.ref_id();   // chromosome index stored on the variant
        let bucket = self.per_ref.entry(ref_id).or_insert_with(Vec::new);
        // Find insertion point by position + allele for stable deterministic order.
        let idx = bucket.partition_point(|existing| {
            existing.pos() < v.pos()
                || (existing.pos() == v.pos() && existing.ref_allele() <= v.ref_allele())
        });
        bucket.insert(idx, v);
    }

    /// Insert a batch of expanded equivalent variants.
    ///
    /// Delegates to `insert` per element; the caller should call
    /// `finalize()` after all records have been processed.
    pub(crate) fn insert_expanded(&mut self, variants: Vec<V>) {
        for v in variants {
            self.insert(v);
        }
    }

    /// Deduplicate exact-position duplicates introduced by expansion.
    ///
    /// Two equivalents for two *different* VCF indels that happen to
    /// map to the same (pos, ref, alt) are kept (both contribute to
    /// rescue scoring).  Exact byte-for-byte duplicates from the same
    /// VCF record are removed.
    pub(crate) fn dedup(&mut self) {
        for bucket in self.per_ref.values_mut() {
            bucket.dedup_by(|a, b| {
                a.pos()         == b.pos()
                    && a.ref_allele() == b.ref_allele()
                    && a.alt_allele() == b.alt_allele()
            });
        }
    }
}
