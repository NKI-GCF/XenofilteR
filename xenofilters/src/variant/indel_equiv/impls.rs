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
            ref_id:            self.ref_id,
            pos:               pos_0based,
            ref_a:             ref_a.to_vec(),
            alt_a:             alt_a.to_vec(),
            // -- Scoring parameters propagated unchanged ------------------
            genotype_quality:  self.genotype_quality,
            is_called:         self.is_called,
        }
    }
}

// ---------------------------------------------------------------------------
// WithAlleles for Population
// ---------------------------------------------------------------------------

impl WithAlleles for Population {
    fn with_alleles(&self, pos_0based: usize, ref_a: &[u8], alt_a: &[u8]) -> Self {
        Self {
            ref_id:            self.ref_id,
            pos:               pos_0based,
            ref_a:             ref_a.to_vec(),
            alt_a:             alt_a.to_vec(),
            allele_frequency:  self.allele_frequency,
        }
    }
}
