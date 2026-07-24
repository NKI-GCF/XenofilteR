use crate::Error;
use crate::variant::Variant;
use noodles::vcf::Header;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record_buf::RecordBuf;
use noodles::vcf::variant::record_buf::info::field::{Value, value::Array::Float};

#[derive(Debug, Clone)]
pub(crate) struct Population {
    pub(crate) ref_id: usize,
    pub(crate) pos: usize,
    pub(crate) ref_a: Vec<u8>,
    pub(crate) alt_a: Vec<u8>,
    /// Allele frequency, e.g., 0.01 (1%)
    pub(crate) allele_frequency: f64,
}

impl Variant for Population {
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
        self.allele_frequency
    }
}

/// Example parser for Population VCF (checks INFO tag "AF")
pub(crate) fn parse_population_record(
    record: &mut RecordBuf,
    header: &Header,
    min_af: f64,
) -> Result<Vec<Population>, Error> {
    let core = crate::variant::parse_core::parse_variant_core(record, header)?;

    // Parse AF: scalar for single-ALT records, array for multi-ALT.
    // NEEDS VERIFICATION: Array::Float iterator API for noodles 0.111.0.
    let afs: Vec<f64> = match record.info().get("AF").flatten() {
        Some(Value::Float(af)) => vec![*af as f64],
        Some(Value::Array(Float(fs))) => fs.iter().filter_map(|&r| r).map(|f| f as f64).collect(),
        _ => return Err(Error::MissingOrInvalidAfTag),
    };

    // Pair each ALT with its per-allele AF.
    // If only one AF is provided (non-standard), apply it to all ALTs.
    core.alts
        .iter()
        .enumerate()
        .filter_map(|(i, alt_a)| {
            let af = match afs.get(i).copied().or_else(|| afs.first().copied()) {
                Some(af) => af,
                None => return Some(Err(Error::MissingOrInvalidAfTag)),
            };
            if af < min_af {
                // Below cutoff -- cannot rescue given AF <= 1.0 requires
                // p_variant > 0.5 to fire, and per-flag semantics this
                // variant is not the near-fixed diagnostic marker the tool
                // targets. Skip loading rather than silently keeping dead
                // weight in the store.
                return None;
            }
            Some(Ok(Population {
                ref_id: core.ref_id,
                pos: core.pos,
                ref_a: core.ref_a.clone(),
                alt_a: alt_a.clone(),
                allele_frequency: af,
            }))
        })
        .collect()
}
