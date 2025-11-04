use anyhow::{Result, anyhow};
use rust_htslib::bcf::{Read, Reader, record::Record};
use std::collections::HashMap;
use std::path::Path;

use crate::MAX_Q;

/// Trait for any object that can be scored against an alignment.
pub trait Variant {
    /// The 1-based reference position of the variant.
    fn pos(&self) -> i64;

    /// The reference allele
    fn ref_allele(&self) -> &[u8];

    /// The alternate allele
    fn alt_allele(&self) -> &[u8];

    /// Checks if a read chunk matches this variant's ALT allele.
    /// Default implementation handles simple string equality.
    fn matches_alt(&self, read_bases: &[u8]) -> bool {
        self.alt_allele() == read_bases
    }

    /// Provides an adjusted score for a read chunk that **matches the ALT allele**.
    /// This is called when the read *disagrees* with the reference but *agrees* with the variant.
    fn score_alt_match(&self, quals: &[u8], log_likelihood_mismatch: &[f64; MAX_Q + 2]) -> f64;

    /// Provides an adjusted score for a read chunk that **matches the REF allele**.
    /// This is called when a variant is present, but the read *agrees* with the reference.
    fn score_ref_match(&self, quals: &[u8], log_likelihood_mismatch: &[f64; MAX_Q + 2]) -> f64;
}

/// Generic VCF reader that accepts a parser function.
pub fn vcf_reader<V: Variant>(
    f: &Path,
    // The parser function defines which use case (Population vs. Sample)
    parser: impl Fn(&mut Record) -> Result<Vec<V>>,
) -> Result<Vec<HashMap<i64, Vec<V>>>> {
    let mut chr_variants: Vec<HashMap<i64, Vec<V>>> = Vec::new();
    let mut bcf_reader = Reader::from_path(f)
        .map_err(|e| anyhow!("Failed to open BCF file {}: {}", f.display(), e))?;

    for record_result in bcf_reader.records() {
        let mut record = record_result?;
        let rid = record.rid().ok_or_else(|| anyhow!("Record missing RID"))? as usize;
        let pos = record.pos();

        // Delegate parsing logic to the provided function
        let variants = parser(&mut record)?;

        while chr_variants.len() <= rid {
            chr_variants.push(HashMap::new());
        }
        chr_variants[rid].insert(pos, variants);
    }
    Ok(chr_variants)
}
