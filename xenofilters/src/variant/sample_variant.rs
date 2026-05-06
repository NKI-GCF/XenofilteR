use crate::variant::Variant;
use anyhow::Result;
use rust_htslib::bcf::record::Record;

// FIXME, a variant could have multiple ALT alleles, and the GT could be 0/2, so we should ideally
// have one SampleVariant per ALT allele, and check if each ALT allele is present in the GT. For
// simplicity, this example assumes only one ALT allele and GT is either 0/1 or 1/1 for the ALT.
pub struct SampleVariant {
    pos: i64,
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
    /// Genotype Quality (GQ) from the FORMAT field
    genotype_quality: f64,
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
    fn p_variant(&self) -> f64 {
        let p_gt_correct = 1.0 - 10f64.powf(-self.genotype_quality / 10.0);
        if self.is_called {
            p_gt_correct
        } else {
            1.0 - p_gt_correct
        }
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
            genotype_quality: gq as f64,
            is_called,
        });
    }
    Ok(variants)
}
