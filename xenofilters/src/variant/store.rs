use crate::variant::{Variant, Eval};
use anyhow::{Result, anyhow};
use noodles::bcf::{record::Record, io::reader::Builder};
use std::path::Path;
use smallvec::SmallVec;
use std::collections::HashMap;
use noodles::vcf::Header;

pub(crate) trait StoreTrait {
    fn overlapping_multi<'s>(&'s self, rid: usize, start: usize, end: usize) -> SmallVec<[Eval<'s>; 0]>;
}

impl<V: Variant> StoreTrait for Store<V> {
    fn overlapping_multi<'s>(&'s self, id: usize, start: usize, end: usize) -> SmallVec<[Eval<'s>; 0]> {
        let mut hits = SmallVec::new();
        for  v in self.overlapping(id, start, end) {
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
    /// Maximum reference span of any variant
    pub(crate) max_variant_len: usize,
}

impl<V: Variant> Store<V> {
    pub(crate) fn new(
        f: &Path,
        parser: impl Fn(&mut Record, &Header) -> Result<Vec<V>>,
    ) -> Result<Store<V>> {
        let mut bcf_reader = Builder::default().build_from_path(f)
            .map_err(|e| anyhow!("Failed to open VCF/BCF {}: {}", f.display(), e))?;

        let mut per_chr = HashMap::new();
        let mut max_variant_len: usize = 1; // at least 1
        let mut is_sorted = true;
        let header = bcf_reader.read_header()?;

        for record_result in bcf_reader.records() {
            let mut record = record_result?;
            let id = record.reference_sequence_id()?;
            let variants = parser(&mut record, &header)?;

            let mut last_pos = None;

            for v in variants {
                let pos = v.pos();
                if is_sorted {
                    if let Some(last) = last_pos
                        && pos < last {
                            eprintln!("Variants in {} are not sorted by position.", f.display());
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
            // Sort each chromosome once, for binary search.
            for chr in per_chr.values_mut() {
                chr.sort_by_key(|v| v.pos());
            }
        }

        Ok(Store {
            per_chr,
            max_variant_len,
        })
    }
    pub(crate) fn overlapping(&self, id: usize, read_start: usize, read_end: usize) -> SmallVec<[&V; 0]> {
        let Some(chr_vars) = self.per_chr.get(&id) else {
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
