use crate::variant::Variant;
use noodles::bcf::record::Record;
use noodles::vcf::{Header, variant::record::info::field::Value};
use crate::{Error, Result};

pub(crate) struct Population {
    pos: usize,
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
    /// Allele frequency, e.g., 0.01 (1%)
    allele_frequency: f64,
}

impl Variant for Population {
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
        self.allele_frequency
    }
}

/// Example parser for Population VCF (checks INFO tag "AF")
pub(crate) fn parse_population_record(record: &mut Record, header: &Header) -> Result<Vec<Population>, Error> {

    let pos = record.variant_start().transpose()?.map(|p| p.get()).unwrap_or(0);
    // FIXME a variant could have multiple ALT alleles, but for simplicity we only consider one here.
    // We can extend this later if needed.
    let alleles = record.alternate_bases();
    let alleles = alleles.as_ref();
    if alleles.contains(&b',') {
        return Err(Error::MultipleAltAllelesNotSupported);
    }


    // FIXME: try retrieving from samples first.
    // 1. Get AF from INFO
    let allele_frequency = match record.info().get(header, "AF").transpose()?.flatten() {
        Some(Value::Float(af)) => af as f64,
        _ => return Err(Error::MissingOrInvalidAfTag),
    };
    let alt_a = alleles.to_vec();
    let ref_a = record.reference_bases().as_ref().to_vec();

    Ok(vec![Population { pos, ref_a, alt_a, allele_frequency }])
}
