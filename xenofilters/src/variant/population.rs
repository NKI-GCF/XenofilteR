use crate::Error;
use crate::variant::Variant;
use noodles::bcf::record::Record;
use noodles::vcf::{Header, variant::record::info::field::Value};

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
pub(crate) fn parse_population_record(
    record: &mut Record,
    header: &Header,
) -> Result<Vec<Population>, Error> {
    let pos    = record.variant_start().transpose()?.map(|p| p.get()).unwrap_or(0);
    let ref_a  = record.reference_bases().as_ref().to_vec();

    // Collect all ALT alleles (comma-separated in BCF alternate_bases).
    let alt_bytes = record.alternate_bases();
    let alts: Vec<Vec<u8>> = alt_bytes.as_ref()
        .split(|&b| b == b',')
        .filter(|a| !a.is_empty())
        .map(|a| a.to_vec())
        .collect();

    if alts.is_empty() {
        return Err(Error::MissingOrInvalidAfTag);
    }

    // Parse AF: scalar for single-ALT records, array for multi-ALT.
    // NEEDS VERIFICATION: Array::Float iterator API for noodles 0.111.0.
    let afs: Vec<f64> = match record.info().get(header, "AF").transpose()?.flatten() {
        Some(Value::Float(af)) => vec![af as f64],
        Some(Value::Array(arr)) => match arr {
            noodles::vcf::variant::record::info::field::value::Array::Float(fs) => {
                fs.iter()
                    .filter_map(|r| r.ok().flatten())
                    .map(|f| f as f64)
                    .collect()
            }
            _ => return Err(Error::MissingOrInvalidAfTag),
        },
        _ => return Err(Error::MissingOrInvalidAfTag),
    };

    // Pair each ALT with its per-allele AF.
    // If only one AF is provided (non-standard), apply it to all ALTs.
    alts.iter()
        .enumerate()
        .map(|(i, alt_a)| {
            let af = afs.get(i)
                .copied()
                .or_else(|| afs.first().copied())
                .ok_or(Error::MissingOrInvalidAfTag)?;
            Ok(Population {
                pos,
                ref_a: ref_a.clone(),
                alt_a: alt_a.clone(),
                allele_frequency: af,
            })
        })
        .collect()
}
