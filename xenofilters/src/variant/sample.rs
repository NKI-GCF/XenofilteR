use crate::Error;
use crate::variant::Variant;
use noodles::bcf::record::Record;
use noodles::vcf::Header;
use noodles::vcf::variant::record::samples::{Sample as NoodlesSample, series::Value};

// FIXME, a variant could have multiple ALT alleles, and the GT could be 0/2, so we should ideally
// have one Sample per ALT allele, and check if each ALT allele is present in the GT. For
// simplicity, this example assumes only one ALT allele and GT is either 0/1 or 1/1 for the ALT.
pub(crate) struct Sample {
    pos: usize,
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
    /// Genotype Quality (GQ) from the FORMAT field
    genotype_quality: f64,
    /// True if GT is 0/1 or 1/1
    is_called: bool,
}

impl Variant for Sample {
    fn pos(&self) -> usize {
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
pub(crate) fn parse_sample_record(
    record: &mut Record,
    header: &Header,
) -> Result<Vec<Sample>, Error> {
    // Genotype representation as a vector of GenotypeAllele.
    // 1. Get GT and GQ from FORMAT
    let samples = record.samples()?;
    let mut it = samples.iter();
    let sample = it.next().ok_or(Error::NoSampleData)?;
    if it.next().is_some() {
        return Err(Error::MultipleSamplesNotSupported);
    }
    let gq = match sample.get(header, "GQ").transpose()? {
        Some(Some(Value::Integer(gq))) => gq,
        _ => return Err(Error::MissingOrInvalidGqTag),
    };

    let gt = match sample.get(header, "GT").transpose()? {
        Some(Some(Value::Integer(gt))) => gt,
        _ => return Err(Error::MissingOrInvalidGtTag),
    };

    let alleles = record.alternate_bases();
    let alleles = alleles.as_ref().split(|&b| b == b',');
    let ref_a = record.reference_bases().as_ref().to_vec();
    let mut variants = Vec::new();

    for (i, alt_a) in alleles.enumerate() {
        let alt_index = (i + 1) as i32;
        // Check if this allele (alt_index) is present in the GT
        let is_called = gt == alt_index;

        variants.push(Sample {
            pos: record
                .variant_start()
                .transpose()?
                .map(|p| p.get())
                .unwrap_or(0),
            ref_a: ref_a.clone(),
            alt_a: alt_a.to_vec(),
            genotype_quality: gq as f64,
            is_called,
        });
    }
    Ok(variants)
}
