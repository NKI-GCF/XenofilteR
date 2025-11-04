use crate::vcf_format::Variant;
use crate::{LOG_LIKELIHOOD_MATCH, MAX_Q};
use anyhow::Result;
use rust_htslib::bcf::record::Record;

#[allow(dead_code)]
pub struct SampleVariant {
    pos: i64,
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
    /// Genotype Quality (GQ) from the FORMAT field
    genotype_quality: f32,
    /// True if GT is 0/1 or 1/1
    is_called: bool,
}

impl Variant for SampleVariant {
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
        let len = quals.len();
        if len == 0 {
            return 0.0;
        }
        let len = len as f64;
        // P(Genotype is correct)
        let p_gt_correct = 1.0 - 10f64.powf(-(self.genotype_quality as f64) / 10.0);

        // P(Variant is truth) = P(GT is correct) if ALT is called,
        // OR 1-P(GT is correct) if REF is called.
        let p_variant = if self.is_called {
            p_gt_correct
        } else {
            1.0 - p_gt_correct
        };

        let mut score_match = 0.0;
        let mut score_mismatch = 0.0;
        for q in quals {
            score_match += LOG_LIKELIHOOD_MATCH[*q as usize];
            score_mismatch += log_likelihood_mismatch[*q as usize];
        }
        p_variant * (score_match / len) + (1.0 - p_variant) * (score_mismatch / len)
    }

    fn score_ref_match(&self, quals: &[u8], log_likelihood_mismatch: &[f64; MAX_Q + 2]) -> f64 {
        let len = quals.len();
        if len == 0 {
            return 0.0;
        }
        let len = len as f64;
        let p_gt_correct = 1.0 - 10f64.powf(-(self.genotype_quality as f64) / 10.0);
        let p_variant = if self.is_called {
            p_gt_correct
        } else {
            1.0 - p_gt_correct
        };

        let mut score_match = 0.0;
        let mut score_mismatch = 0.0;
        for q in quals {
            score_match += LOG_LIKELIHOOD_MATCH[*q as usize];
            score_mismatch += log_likelihood_mismatch[*q as usize];
        }
        (1.0 - p_variant) * (score_match / len) + p_variant * (score_mismatch / len)
    }
}

/// Example parser for Sample-Specific VCF (checks FORMAT tags "GT" and "GQ")
pub fn parse_sample_record(record: &mut Record) -> Result<Vec<SampleVariant>> {
    // Genotype representation as a vector of GenotypeAllele.
    // 1. Get GT and GQ from FORMAT
    let gq = record
        .format(b"GQ")
        .integer()?
        .first()
        .and_then(|v| v.first().map(|i| *i as f32))
        .unwrap_or(0.0);
    let gt = record
        .format(b"GT")
        .integer()?
        .first()
        .and_then(|v| v.first());

    let alleles = record.alleles();
    let ref_a = alleles[0].to_vec();
    let mut variants = Vec::new();

    for (i, alt_a) in alleles[1..].iter().enumerate() {
        let alt_index = (i + 1) as i32;
        // Check if this allele (alt_index) is present in the GT
        let is_called = gt.is_some_and(|&g| g == alt_index);

        variants.push(SampleVariant {
            pos: record.pos(),
            ref_a: ref_a.clone(),
            alt_a: alt_a.to_vec(),
            genotype_quality: gq,
            is_called,
        });
    }
    Ok(variants)
}
