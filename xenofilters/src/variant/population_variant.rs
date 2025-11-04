use crate::vcf_format::Variant;
use crate::{LOG_LIKELIHOOD_MATCH, MAX_Q};
use anyhow::{Result, anyhow};
use rust_htslib::bcf::record::Record;

pub struct PopulationVariant {
    pos: i64,
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
    /// Allele frequency, e.g., 0.01 (1%)
    allele_frequency: f32,
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

    fn score_alt_match(&self, quals: &[u8], log_likelihood_mismatch: &[f64; MAX_Q + 2]) -> f64 {
        let p_variant = self.allele_frequency as f64;

        // This read matches the ALT.
        // Score = P(Variant is truth) * (Score for Match) + P(Ref is truth) * (Score for Mismatch)
        let len = quals.len();
        if len == 0 {
            return 0.0;
        }
        let len = len as f64;
        let mut score_match = 0.0;
        let mut score_mismatch = 0.0;
        for q in quals {
            score_match += LOG_LIKELIHOOD_MATCH[*q as usize];
            score_mismatch += log_likelihood_mismatch[*q as usize];
        }
        p_variant * (score_match / len) + (1.0 - p_variant) * (score_mismatch / len)
    }

    fn score_ref_match(&self, quals: &[u8], log_likelihood_mismatch: &[f64; MAX_Q + 2]) -> f64 {
        let p_variant = self.allele_frequency as f64;

        // This read matches the REF.
        // Score = P(Ref is truth) * (Score for Match) + P(Variant is truth) * (Score for Mismatch)

        let len = quals.len();
        if len == 0 {
            return 0.0;
        }
        let len = len as f64;

        let mut score_match = 0.0;
        let mut score_mismatch = 0.0;
        for q in quals {
            score_match += LOG_LIKELIHOOD_MATCH[*q as usize];
            score_mismatch += log_likelihood_mismatch[*q as usize];
        }
        (1.0 - p_variant) * (score_match / len) + p_variant * (score_mismatch / len)
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
            allele_frequency: af_values[i],
        });
    }
    Ok(variants)
}
