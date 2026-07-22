//! Generic sorted-interval container + shared BED/VCF loading loop.
//! Backs AmbiguousRegions, ScoredRegions, SegregateVariants, and Store<V>.

use crate::Error;
use std::collections::HashMap;

/// Anything with a half-open reference span `[start, end)`.
pub(crate) trait Interval {
    fn start(&self) -> usize;
    fn end(&self) -> usize;
}

/// Per-reference-id, sorted-by-start collection of intervals with
/// binary-search overlap queries. Coordinates are 0-based half-open.
#[derive(Debug)]
pub(crate) struct IntervalStore<T: Interval> {
    per_ref: Vec<Vec<T>>,
    max_span: usize,
}

impl<T: Interval> Default for IntervalStore<T> {
    fn default() -> Self {
        Self {
            per_ref: Vec::new(),
            max_span: 1,
        }
    }
}

impl<T: Interval> IntervalStore<T> {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn insert(&mut self, ref_id: usize, item: T) {
        let span = item.end().saturating_sub(item.start()).max(1);
        self.max_span = self.max_span.max(span);
        if self.per_ref.len() <= ref_id {
            self.per_ref.resize_with(ref_id + 1, Vec::new);
        }
        self.per_ref[ref_id].push(item);
    }

    /// Call once after all inserts. O(n log n) per ref bucket.
    pub(crate) fn sort(&mut self) {
        for v in &mut self.per_ref {
            v.sort_unstable_by_key(|i| i.start());
        }
    }

    /// Items overlapping `[start, end)`. `max_span` bounds how far back
    /// the binary search needs to look (an interval starting before
    /// `start - max_span` cannot possibly reach into `start`).
    pub(crate) fn overlapping(
        &self,
        ref_id: usize,
        start: usize,
        end: usize,
    ) -> impl Iterator<Item = &T> {
        let items: &[T] = self.per_ref.get(ref_id).map(Vec::as_slice).unwrap_or(&[]);
        let lower = start.saturating_sub(self.max_span);
        let lo = items.partition_point(|i| i.start() < lower);
        items[lo..]
            .iter()
            .take_while(move |i| i.start() < end)
            .filter(move |i| i.end() > start)
    }

    pub(crate) fn overlaps(&self, ref_id: usize, start: usize, end: usize) -> bool {
        self.overlapping(ref_id, start, end).next().is_some()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.per_ref.iter().all(Vec::is_empty)
    }
    pub(crate) fn dedup_by(&mut self, mut same: impl FnMut(&T, &T) -> bool) {
        for v in &mut self.per_ref {
            v.dedup_by(|a, b| same(a, b));
        }
    }
    #[cfg(test)]
    pub(crate) fn per_ref(&self) -> &Vec<Vec<T>> {
        &self.per_ref
    }
}

/// Shared "open records, resolve chrom→ref_id, skip on miss, insert, sort" loop.
///
/// `build` returns anything iterable (`Option<T>` for 0-or-1 items,
/// `Vec<T>`/`SmallVec<T>` for indel-equivalence-style expansion into
/// multiple items per record).
pub(crate) fn load_into_store<T, Rec, F, C, I>(
    records: impl Iterator<Item = Result<Rec, Error>>,
    name_to_id: &HashMap<String, usize>,
    mut chrom_of: C,
    mut build: F,
) -> Result<IntervalStore<T>, Error>
where
    T: Interval,
    I: IntoIterator<Item = T>,
    C: FnMut(&Rec) -> String,
    F: FnMut(&Rec, usize) -> Result<I, Error>,
{
    let mut store = IntervalStore::new();
    for rec in records {
        let rec = rec?;
        let chrom = chrom_of(&rec);
        let Some(&ref_id) = name_to_id.get(&chrom) else {
            continue;
        };
        for item in build(&rec, ref_id)? {
            store.insert(ref_id, item);
        }
    }
    store.sort();
    Ok(store)
}
