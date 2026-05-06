mod population_variant;
mod sample_variant;
mod variant_store;
mod variant_eval;

pub(super) use population_variant::{PopulationVariant, parse_population_record};
pub(super) use sample_variant::{SampleVariant, parse_sample_record};
pub(super) use variant_store::{VariantStore, VariantStoreTrait};
use anyhow::{Result, anyhow};
use rust_htslib::bcf::{self, Read};
use std::path::Path;
use smallvec::SmallVec;
pub(crate) use variant_eval::VariantEval;

pub(crate) type DeltaPerRec<'a> = SmallVec<[VariantEval<'a>; 0]>;
pub(crate) type VntPerRec<'a> = SmallVec<[DeltaPerRec<'a>; 2]>;

/// Trait for any object that can be scored against an alignment.
pub trait Variant: Sync + Send {
    /// The 1-based reference position of the variant.
    fn pos(&self) -> i64;

    /// The reference allele
    fn ref_allele(&self) -> &[u8];

    /// The alternate allele
    fn alt_allele(&self) -> &[u8];

    /// End position (1-based, exclusive) — allows easy overlap check
    fn end(&self) -> i64 {
        self.pos() + self.ref_allele().len() as i64
    }

    /// Does this variant overlap [read_start, read_end) ?
    fn overlaps(&self, read_start: i64, read_end: i64) -> bool {
        let v_start = self.pos();
        let v_end   = self.end();
        v_start < read_end && v_end > read_start
    }
    fn p_variant(&self) -> f64;
}

/// Generic VCF reader that accepts a parser function.
pub fn vcf_reader<V: Variant>(
    f: &Path,
    parser: impl Fn(&mut bcf::Record) -> Result<Vec<V>>,
) -> Result<VariantStore<V>> {
    let mut bcf_reader = bcf::Reader::from_path(f)
        .map_err(|e| anyhow!("Failed to open VCF/BCF {}: {}", f.display(), e))?;

    let mut per_chr: Vec<Vec<V>> = Vec::new();
    let mut max_variant_len: i64 = 1; // at least 1
    let mut is_sorted = true;

    for record_result in bcf_reader.records() {
        let mut record = record_result?;
        let rid = record.rid().ok_or_else(|| anyhow!("Record missing RID"))? as usize;
        let variants = parser(&mut record)?;

        while per_chr.len() <= rid {
            per_chr.push(Vec::new());
        }
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
            per_chr[rid].push(v);
        }
    }
    if !is_sorted {
        // Sort each chromosome once, for binary search.
        for chr in &mut per_chr {
            chr.sort_by_key(|v| v.pos());
        }
    }

    Ok(VariantStore {
        per_chr,
        max_variant_len,
    })
}
