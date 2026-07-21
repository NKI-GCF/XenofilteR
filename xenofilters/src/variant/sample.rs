use crate::variant::Variant;
use crate::Error;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record_buf::samples::sample::Value;
use noodles::vcf::variant::record_buf::RecordBuf;
use noodles::vcf::Header;

// FIXME, a variant could have multiple ALT alleles, and the GT could be 0/2, so we should ideally
// have one Sample per ALT allele, and check if each ALT allele is present in the GT. For
// simplicity, this example assumes only one ALT allele and GT is either 0/1 or 1/1 for the ALT.
#[derive(Debug, Clone)]
pub(crate) struct Sample {
    pub(crate) ref_id: usize,
    pub(crate) pos: usize,
    pub(crate) ref_a: Vec<u8>,
    pub(crate) alt_a: Vec<u8>,
    /// Genotype Quality (GQ) from the FORMAT field
    pub(crate) genotype_quality: f64,
    /// True if GT is 0/1 or 1/1
    pub(crate) is_called: bool,
}

impl Variant for Sample {
    fn ref_id(&self) -> usize {
        self.ref_id
    }
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
    record: &mut RecordBuf,
    header: &Header,
) -> Result<Vec<Sample>, Error> {
    let core = crate::variant::parse_core::parse_variant_core(record, header)?;

    let samples = record.samples();
    if samples.get_index(2).is_some() {
        return Err(Error::MultipleSamplesNotSupported);
    }
    let sample = samples.get_index(1).ok_or(Error::NoSampleData)?;
    let gq = match sample.get("GQ").flatten() {
        Some(Value::Integer(gq)) => *gq,
        _ => return Err(Error::MissingOrInvalidGqTag),
    };
    let gt = match sample.get("GT").flatten() {
        Some(Value::Integer(gt)) => *gt,
        _ => return Err(Error::MissingOrInvalidGtTag),
    };

    core.alts
        .iter()
        .enumerate()
        .map(|(i, alt_a)| {
            let is_called = gt == (i + 1) as i32;
            Ok(Sample {
                ref_id: core.ref_id,
                pos: core.pos,
                ref_a: core.ref_a.clone(),
                alt_a: alt_a.to_vec(),
                genotype_quality: gq as f64,
                is_called,
            })
        })
        .collect()
}
