// src/variant/store.rs

use crate::Error;
use crate::variant::{Eval, Variant};
use noodles::bcf::io::reader::Builder;
use noodles::vcf::Header;
use smallvec::SmallVec;
use std::collections::HashMap;
use std::path::Path;
use noodles::vcf::variant::record_buf::RecordBuf;

pub(crate) const VNT_CT: usize = 4;

pub(crate) type EvalVec<'s> = SmallVec<[Eval<'s>; VNT_CT]>;

/// Any object that can answer overlap queries for variant scoring.
///
/// `Send + Sync` are required so that stores can be placed behind `Arc` and
/// shared across scoring worker threads.
pub(crate) trait StoreTrait: Send + Sync {
    fn overlapping_multi<'s>(&'s self, rid: usize, start: usize, end: usize) -> EvalVec<'s>;
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
    pub(crate) per_chr: HashMap<usize, Vec<V>>,
    /// Maximum reference span of any variant in this store.
    pub(crate) max_variant_len: usize,
}

impl<V: Variant> Store<V> {
    pub(crate) fn new() -> Store<V> {
        Store {
            per_chr: HashMap::new(),
            max_variant_len: 1,
        }
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

        let mut per_chr = HashMap::new();
        let mut max_variant_len: usize = 1;
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
                let span = v.end() - pos;
                if span > max_variant_len {
                    max_variant_len = span;
                }
                per_chr.entry(id).or_insert_with(Vec::new).push(v);
            }
        }
        if !is_sorted {
            for chr in per_chr.values_mut() {
                chr.sort_by_key(|v| v.pos());
            }
        }

        Ok(Store {
            per_chr,
            max_variant_len,
        })
    }

    pub(crate) fn overlapping(
        &self,
        id: usize,
        read_start: usize,
        read_end: usize,
    ) -> SmallVec<[&V; 0]> {
        let Some(chr_vars) = self.per_chr.get(&id) else {
            return SmallVec::new();
        };

        let lower = read_start.saturating_sub(self.max_variant_len);
        let upper = read_end;

        let start_idx = chr_vars.partition_point(|v| v.pos() < lower);
        let end_idx = chr_vars.partition_point(|v| v.pos() < upper);

        let mut hits = SmallVec::new();
        for v in &chr_vars[start_idx..end_idx] {
            if v.overlaps(read_start, read_end) {
                hits.push(v);
            }
        }
        hits
    }
    /// Insert one variant into the per-chromosome sorted bucket.
    ///
    /// `ref_id` is the 0-based chromosome index from the reference header.
    /// Uses `partition_point` (binary search) to find the insertion index;
    /// O(log n + n) per call — acceptable at startup, never called on the
    /// hot path.
    pub(crate) fn insert(&mut self, ref_id: usize, v: V) {
        let bucket = self.per_chr.entry(ref_id).or_default();
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
        let bucket = self.per_chr.entry(ref_id).or_default();
        bucket.reserve(variants.len());

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
        for bucket in self.per_chr.values_mut() {
            bucket.dedup_by(|a, b| {
                // b comes first in sorted order; dedup_by retains b.
                a.pos() == b.pos()
                    && a.ref_allele() == b.ref_allele()
                    && a.alt_allele() == b.alt_allele()
            });
        }
    }
}

#[cfg(test)]
mod tests;
