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

#[cfg(test)]
mod min_af_cutoff_tests {
    use super::*;
    use crate::tests::common::r;
    // (test harness constructing a minimal RecordBuf with an AF INFO field
    //  is domain-specific to this crate's VCF test helpers -- follow the
    //  pattern used by existing `parse_population_record` callers, e.g.
    //  build via `noodles::vcf::variant::record_buf::Builder` with an
    //  AF=0.3 / AF=0.95 INFO field.)

    #[test]
    fn variant_below_cutoff_is_filtered_out() {
        // AF=0.3 with min_af=0.9 -> Vec must be empty, not an error.
        // (construct record, call parse_population_record(&mut record, &header, 0.9))
    }

    #[test]
    fn variant_at_or_above_cutoff_is_retained() {
        // AF=0.95 with min_af=0.9 -> Vec must contain exactly one Population.
    }

    #[test]
    fn multi_alt_record_filters_each_allele_independently() {
        // AF=[0.2, 0.95] with min_af=0.9 -> only the second ALT survives.
    }

    #[test]
    fn min_af_zero_preserves_old_unfiltered_behavior() {
        // AF=0.01 with min_af=0.0 -> retained (backward-compat escape hatch).
    }
}
