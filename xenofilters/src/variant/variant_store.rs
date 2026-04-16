use smallvec::SmallVec;
use crate::variant::Variant;

pub(crate) trait VariantStoreTrait {
    fn variants_overlapping_multi<'s>(
        &'s self,
        ranges: &[(usize, i64, i64)],
    ) -> SmallVec<[&'s dyn Variant; 8]>;
}

impl<V: Variant> VariantStoreTrait for VariantStore<V> {
    fn variants_overlapping_multi<'s>(
        &'s self,
        ranges: &[(usize, i64, i64)],
    ) -> SmallVec<[&'s dyn Variant; 8]> {
        let mut hits = SmallVec::new();
        for &(rid, start, end) in ranges {
            for v in self.variants_overlapping(rid, start, end) {
                // Coerce &V → &dyn Variant
                if !hits.iter().any(|h :&&dyn Variant| h.pos() == v.pos()) {  // cheap dedup by pos
                    hits.push(v as &dyn Variant);
                }
            }
        }
        hits
    }
}

/// `Vec<V>` per chrom, sorted by `pos`.
#[derive(Debug)]
pub(crate) struct VariantStore<V: Variant> {
    pub(crate) per_chr: Vec<Vec<V>>,
    /// Maximum reference span of any variant
    pub(crate) max_variant_len: i64,
}

impl<V: Variant> VariantStore<V> {
    /// Returns a small borrow-only slice of variants that **overlap** any part of the read.
    pub(crate) fn variants_overlapping(
        &self,
        rid: usize,          // htslib chromosome id
        read_start: i64,     // 0-based, inclusive
        read_end: i64,       // exclusive
    ) -> SmallVec<[&V; 8]> {
        let Some(chr_vars) = self.per_chr.get(rid) else {
            return SmallVec::new();
        };

        // Look back far enough for any variant that could touch the read
        let lower = read_start.saturating_sub(self.max_variant_len);
        let upper = read_end;

        // partition_point is stable since Rust 1.52 (edition 2021+)
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

    /// Convenience: union of several ranges (mates, supplementaries, etc.).
    /// Still returns a single SmallVec — almost always 0–5 variants.
    pub(crate) fn variants_overlapping_multi(
        &self,
        ranges: &[(usize, i64, i64)],   // (rid, start, end) tuples
    ) -> SmallVec<[&V; 8]> {
        let mut hits = SmallVec::new();
        for &(rid, start, end) in ranges {
            for v in self.variants_overlapping(rid, start, end) {
                if !hits.iter().any(|h: &&V| h.pos() == v.pos()) {  // cheap dedup by pos
                    hits.push(v);
                }
            }
        }
        hits
    }
}
