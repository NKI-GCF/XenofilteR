use crate::variant::Variant;
use anyhow::{Result, anyhow};
use rust_htslib::bcf::record::Record;

/// FIXME a variant could have multiple ALT alleles, but for simplicity we only consider one here.
/// We can extend this later if needed.
pub struct PopulationVariant {
    pos: i64,
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
    /// Allele frequency, e.g., 0.01 (1%)
    allele_frequency: f64,
}

impl Variant for PopulationVariant {
    fn pos(&self) -> i64 {
        self.pos
    }
    fn ref_allele(&self) -> &[u8] {
        &self.ref_a
    }
    fn alt_allele(&self) -> &[u8] {
        &self.alt_a
    }
    fn p_variant(&self) -> f64 {
        self.allele_frequency
    }
}

/// Example parser for Population VCF (checks INFO tag "AF")
pub fn parse_population_record(record: &mut Record) -> Result<Vec<PopulationVariant>> {
    // 1. Get AF from INFO
    let af_values = record
        .info(b"AF")
        .float()?
        .ok_or_else(|| anyhow!("Missing AF tag"))?;

    let alleles = record.alleles();
    let ref_a = alleles[0].to_vec();
    let mut variants = Vec::new();

    for (i, alt_a) in alleles[1..].iter().enumerate() {
        variants.push(PopulationVariant {
            pos: record.pos(),
            ref_a: ref_a.clone(),
            alt_a: alt_a.to_vec(),
            allele_frequency: af_values[i] as f64,
        });
    }
    Ok(variants)
}
