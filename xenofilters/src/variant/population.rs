use crate::variant::Variant;
use crate::Error;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record_buf::info::field::{value::Array::Float, Value};
use noodles::vcf::variant::record_buf::RecordBuf;
use noodles::vcf::Header;

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

// src/variant/population.rs
#[cfg(test)]
mod min_af_cutoff_tests {
    use super::*;
    use noodles::core::Position;
    use noodles::vcf::header::record::value::{map::Contig, Map};
    use noodles::vcf::variant::record_buf::info::field::{
        value::Array as InfoArray, Value as InfoValue,
    };
    use noodles::vcf::variant::record_buf::RecordBuf;
    use noodles::vcf::Header;

    /// Minimal single-contig header sufficient for `parse_variant_core`'s
    /// `header.contigs().get_index_of(&chrom)` lookup. Not used for real
    /// file I/O, so contig length and other write-required fields are
    /// deliberately omitted.
    fn test_header() -> Header {
        Header::builder()
            .add_contig("chr1", Map::<Contig>::builder().build().unwrap())
            .build()
    }

    fn population_record(ref_a: &str, alt_a: &str, af: f32) -> RecordBuf {
        let mut record = RecordBuf::builder()
            .set_reference_sequence_name("chr1")
            .set_variant_start(Position::try_from(100).unwrap())
            .set_reference_bases(ref_a)
            .set_alternate_bases(vec![alt_a.to_string()].into())
            .build();
        record
            .info_mut()
            .insert(String::from("AF"), Some(InfoValue::Float(af)));
        record
    }

    fn population_record_multi_alt(ref_a: &str, alts: &[&str], afs: &[f32]) -> RecordBuf {
        let mut record = RecordBuf::builder()
            .set_reference_sequence_name("chr1")
            .set_variant_start(Position::try_from(100).unwrap())
            .set_reference_bases(ref_a)
            .set_alternate_bases(
                alts.iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>()
                    .into(),
            )
            .build();
        let af_values: Vec<Option<f32>> = afs.iter().map(|&f| Some(f)).collect();
        record.info_mut().insert(
            String::from("AF"),
            Some(InfoValue::Array(InfoArray::Float(af_values))),
        );
        record
    }

    #[test]
    fn variant_below_cutoff_is_filtered_out() {
        let header = test_header();
        let mut record = population_record("A", "G", 0.3);
        let result = parse_population_record(&mut record, &header, 0.9).unwrap();
        assert!(
            result.is_empty(),
            "AF=0.3 must be filtered out at min_af=0.9"
        );
    }

    #[test]
    fn variant_at_or_above_cutoff_is_retained() {
        let header = test_header();
        let mut record = population_record("A", "G", 0.95);
        let result = parse_population_record(&mut record, &header, 0.9).unwrap();
        assert_eq!(result.len(), 1, "AF=0.95 must be retained at min_af=0.9");
        assert!((result[0].p_variant() - 0.95).abs() < 1e-6);
    }

    #[test]
    fn multi_alt_record_filters_each_allele_independently() {
        let header = test_header();
        let mut record = population_record_multi_alt("A", &["G", "T"], &[0.2, 0.95]);
        let result = parse_population_record(&mut record, &header, 0.9).unwrap();
        assert_eq!(
            result.len(),
            1,
            "only the AF=0.95 allele should survive min_af=0.9"
        );
        assert_eq!(result[0].alt_a, b"T");
    }

    #[test]
    fn min_af_zero_preserves_old_unfiltered_behavior() {
        let header = test_header();
        let mut record = population_record("A", "G", 0.01);
        let result = parse_population_record(&mut record, &header, 0.0).unwrap();
        assert_eq!(
            result.len(),
            1,
            "min_af=0.0 must retain even very rare variants"
        );
    }
}
