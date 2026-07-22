// src/variant/store.rs

use crate::Error;
use crate::variant::{Eval, Variant};
use noodles::bcf::io::reader::Builder;
use noodles::vcf::Header;
use smallvec::SmallVec;
use std::path::Path;
use noodles::vcf::variant::record_buf::RecordBuf;
use crate::region::interval_store::{Interval, IntervalStore};

pub(crate) const VNT_CT: usize = 4;

pub(crate) type EvalVec<'s> = SmallVec<[Eval<'s>; VNT_CT]>;

/// Any object that can answer overlap queries for variant scoring.
///
/// `Send + Sync` are required so that stores can be placed behind `Arc` and
/// shared across scoring worker threads.
pub(crate) trait StoreTrait: Send + Sync {
    fn overlapping_multi<'s>(&'s self, rid: usize, start: usize, end: usize) -> EvalVec<'s>;
}

impl<V: Variant> Interval for V {
    fn start(&self) -> usize { self.pos() }
    fn end(&self) -> usize { Variant::end(self) }
}

impl<V: Variant> StoreTrait for Store<V> {
    fn overlapping_multi<'s>(&'s self, id: usize, start: usize, end: usize) -> EvalVec<'s> {
        let mut hits = SmallVec::new();
        for v in self.overlapping(id, start, end) {
            let mut eval = Eval::new();
            eval.set_variant(v as &dyn Variant);
            hits.push(eval);
        }
        hits
    }
}

/// `Vec<V>` per chrom, sorted by `pos`.
#[derive(Debug)]
pub(crate) struct Store<V: Variant> {
    inner: IntervalStore<V>,
}

impl<V: Variant> Store<V> {
    pub(crate) fn new() -> Self { Self { inner: IntervalStore::new() } }

    /// Insert one variant into the per-chromosome sorted bucket.
    ///
    /// `ref_id` is the 0-based chromosome index from the reference header.
    /// Uses `partition_point` (binary search) to find the insertion index;
    /// O(log n + n) per call -- acceptable at startup, never called on the
    /// hot path.
    pub(crate) fn insert(&mut self, ref_id: usize, v: V) { self.inner.insert(ref_id, v); }

    /// Insert a batch of pre-expanded variants.  All share the same ref_id.
    /// Delegates to `insert` per element; no separate finalization required
    /// because each insert maintains the sorted invariant.
    pub(crate) fn insert_expanded(&mut self, ref_id: usize, variants: Vec<V>) {
        for v in variants { self.insert(ref_id, v); }
    }
    /// Remove byte-exact duplicates introduced when two VCF records in a
    /// repeat region expand to the same (pos, ref_a, alt_a) tuple.
    ///
    /// Called once after all inserts are complete.  The sorted invariant
    /// ensures adjacent equal elements are the candidates; `dedup_by` is
    /// O(n) after sorting.
    pub(crate) fn dedup(&mut self) {
        self.inner.sort();
        self.inner.dedup_by(|a, b| a.pos() == b.pos() && a.ref_allele() == b.ref_allele() && a.alt_allele() == b.alt_allele());
        // dedup_by needs direct access; add IntervalStore::dedup_per_ref if desired,
        // or keep dedup here operating on inner.per_ref via a crate-visible accessor.
    }
    pub(crate) fn overlapping(&self, id: usize, start: usize, end: usize) -> SmallVec<[&V; 0]> {
        self.inner.overlapping(id, start, end).collect()
    }
    pub(crate) fn new_from_path(
        f: &Path,
        parser: impl Fn(&mut RecordBuf, &Header) -> Result<Vec<V>, Error>,
    ) -> Result<Store<V>, Error> {
        let mut bcf_reader =
            Builder::default()
                .build_from_path(f)
                .map_err(|e| Error::FailedToOpenVcfBcf {
                    path: f.to_path_buf(),
                    source: e,
                })?;

        let mut inner = IntervalStore::new();
        let mut is_sorted = true;
        let header = bcf_reader.read_header()?;

        for record_result in bcf_reader.records() {
            let record = record_result?;
            let mut record_buf = RecordBuf::try_from_variant_record(&header, &record)?;

            let id = record.reference_sequence_id()?;
            let variants = parser(&mut record_buf, &header)?;

            let mut last_pos = None;

            for v in variants {
                let pos = v.pos();
                if is_sorted {
                    if let Some(last) = last_pos
                        && pos < last
                    {
                        tracing::warn!(
                            path = %f.display(),
                            "Variants are not sorted by position; sorting now."
                        );
                        is_sorted = false;
                    }
                    last_pos = Some(pos);
                }
                inner.insert(id, v);
            }
        }
        if !is_sorted {
            inner.sort();
        }

        Ok(Store { inner })
    }

}

#[cfg(test)]
mod tests;
