// src/variant/store.rs

use crate::Error;
use crate::variant::{Eval, Variant};
use noodles::bcf::io::reader::Builder;
use noodles::vcf::Header;
use noodles::vcf::variant::record_buf::RecordBuf;
use rust_lapper::{Interval, Lapper};
use smallvec::SmallVec;
use std::collections::HashMap;
use std::path::Path;

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

/// A generic helper to iterate over records, map their sequence names to IDs,
/// extract intervals, and build a `Lapper` interval tree per chromosome.
pub(crate) fn load_lappers<I, R, FName, FParse, T, Iter>(
    records: I,
    name_to_id: &HashMap<String, usize>,
    mut get_name: FName,
    mut parse_record: FParse,
) -> Result<HashMap<usize, Lapper<usize, T>>, Error>
where
    I: IntoIterator<Item = Result<R, Error>>,
    FName: FnMut(&R) -> String,
    FParse: FnMut(&R, usize) -> Result<Iter, Error>,
    Iter: IntoIterator<Item = Interval<usize, T>>,
    T: Send + Sync + Clone + Eq,
{
    let mut raw_intervals: HashMap<usize, Vec<Interval<usize, T>>> = HashMap::new();

    for record_res in records {
        let record = record_res?;
        let name = get_name(&record);

        if let Some(&ref_id) = name_to_id.get(&name) {
            // Because of IntoIterator, this handles both Some(iv) and vec![iv1, iv2]
            for interval in parse_record(&record, ref_id)? {
                raw_intervals.entry(ref_id).or_default().push(interval);
            }
        }
    }
    Ok(raw_intervals
        .into_iter()
        .map(|(rid, ivs)| (rid, Lapper::new(ivs)))
        .collect())
}

/// `Vec<V>` per chrom, sorted by `pos`.
#[derive(Debug)]
pub(crate) struct Store<V: Variant> {
    per_chr_lapper: HashMap<usize, Lapper<usize, usize>>,
    per_chr_data: HashMap<usize, Vec<V>>,
}

impl<V: Variant> Store<V> {
    pub(crate) fn new() -> Self {
        Self {
            per_chr_lapper: HashMap::new(),
            per_chr_data: HashMap::new(),
        }
    }

    /// Insert one variant into the per-chromosome sorted bucket.
    ///
    /// `ref_id` is the 0-based chromosome index from the reference header.
    /// Uses `partition_point` (binary search) to find the insertion index;
    /// O(log n + n) per call -- acceptable at startup, never called on the
    /// hot path.
    pub(crate) fn insert(&mut self, ref_id: usize, v: V) {
        self.per_chr_data.entry(ref_id).or_default().push(v);
    }

    /// Insert a batch of pre-expanded variants.  All share the same ref_id.
    /// Delegates to `insert` per element; no separate finalization required
    /// because each insert maintains the sorted invariant.
    pub(crate) fn insert_expanded(&mut self, ref_id: usize, variants: Vec<V>) {
        self.per_chr_data
            .entry(ref_id)
            .or_default()
            .extend(variants);
    }

    /// Remove byte-exact duplicates introduced when two VCF records in a
    /// repeat region expand to the same (pos, ref_a, alt_a) tuple.
    ///
    /// Called once after all inserts are complete.  The sorted invariant
    /// ensures adjacent equal elements are the candidates; `dedup_by` is
    /// O(n) after sorting.
    pub(crate) fn dedup(&mut self) {
        for (&ref_id, data) in &mut self.per_chr_data {
            data.sort_by_key(|v| v.pos());
            data.dedup_by(|a, b| {
                a.pos() == b.pos()
                    && a.ref_allele() == b.ref_allele()
                    && a.alt_allele() == b.alt_allele()
            });

            let intervals = data
                .iter()
                .enumerate()
                .map(|(idx, v)| Interval {
                    start: v.pos(),
                    stop: Variant::end(v),
                    val: idx,
                })
                .collect();

            self.per_chr_lapper.insert(ref_id, Lapper::new(intervals));
        }
    }

    pub(crate) fn overlapping(&self, ref_id: usize, start: usize, end: usize) -> SmallVec<[&V; 0]> {
        let mut hits = SmallVec::new();
        if let (Some(lapper), Some(data)) = (
            self.per_chr_lapper.get(&ref_id),
            self.per_chr_data.get(&ref_id),
        ) {
            for interval in lapper.find(start, end) {
                hits.push(&data[interval.val]);
            }
        }
        hits
    }

    pub(crate) fn new_from_path(
        f: &Path,
        parser: impl Fn(&mut RecordBuf, &Header, f64) -> Result<Vec<V>, Error>,
        min_af: f64,
    ) -> Result<Store<V>, Error> {
        let mut bcf_reader =
            Builder::default()
                .build_from_path(f)
                .map_err(|e| Error::FailedToOpenVcfBcf {
                    path: f.to_path_buf(),
                    source: e,
                })?;

        let mut store = Store::new();
        let header = bcf_reader.read_header()?;

        for record_result in bcf_reader.records() {
            let record = record_result?;
            let mut record_buf = RecordBuf::try_from_variant_record(&header, &record)?;

            let id = record.reference_sequence_id()?;
            let variants = parser(&mut record_buf, &header, min_af)?;

            store.insert_expanded(id, variants);
        }

        store.dedup();
        Ok(store)
    }
}

#[cfg(test)]
mod tests;
