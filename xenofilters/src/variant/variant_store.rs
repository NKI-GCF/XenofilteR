use smallvec::SmallVec;
use crate::variant::{Variant, VariantEval};

pub(crate) trait VariantStoreTrait {
    fn variants_overlapping_multi<'s>(
        &'s self,
        ranges: &[(usize, i64, i64)],
    ) -> SmallVec<[VariantEval<'s>; 0]>;
}

impl<V: Variant> VariantStoreTrait for VariantStore<V> {
    fn variants_overlapping_multi<'s>(
        &'s self,
        ranges: &[(usize, i64, i64)],
    ) -> SmallVec<[VariantEval<'s>; 0]> {
        let mut hits = SmallVec::new();
        for &(rid, start, end) in ranges {
            for  v in self.variants_overlapping(rid, start, end) {
                let mut eval = VariantEval::new();
                eval.set_variant(v as &dyn Variant);
                hits.push(eval);
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
    pub(crate) fn variants_overlapping(
        &self,
        rid: usize,
        read_start: i64,
        read_end: i64,
    ) -> SmallVec<[&V; 0]> {
        let Some(chr_vars) = self.per_chr.get(rid) else {
            return SmallVec::new();
        };

        let lower = read_start.saturating_sub(self.max_variant_len);
        let upper = read_end;

        let start_idx = chr_vars.partition_point(|v| v.pos() < lower);
        let end_idx   = chr_vars.partition_point(|v| v.pos() < upper);

        let mut hits = SmallVec::new();
        for v in &chr_vars[start_idx..end_idx] {
            if v.overlaps(read_start, read_end) {
                hits.push(v);
            }
        }
        hits
    }
}
