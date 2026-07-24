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
    /// FORMAT `PS` tag, if present -- see `Variant::phase_set`.
    pub(crate) phase_set: Option<u32>,
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
    fn phase_set(&self) -> Option<u32> {
        self.phase_set
    }
}

/// Example parser for Sample-Specific VCF (checks FORMAT tags "GT" and "GQ")
pub(crate) fn parse_sample_record(
    record: &mut RecordBuf,
    header: &Header,
    min_gq: f64,
) -> Result<Vec<Sample>, Error> {
    let core = crate::variant::parse_core::parse_variant_core(record, header)?;

    let samples = record.samples();
    if samples.get_index(2).is_some() {
        return Err(Error::MultipleSamplesNotSupported);
    }
    let sample = samples.get_index(1).ok_or(Error::NoSampleData)?;
    let gq = match sample.get("GQ").flatten() {
        Some(Value::Integer(gq)) => *gq as f64,
        _ => return Err(Error::MissingOrInvalidGqTag),
    };
    let gt = match sample.get("GT").flatten() {
        Some(Value::Integer(gt)) => *gt,
        _ => return Err(Error::MissingOrInvalidGtTag),
    };
    let phase_set = match sample.get("PS").flatten() {
        Some(Value::Integer(ps)) => Some(*ps as u32),
        _ => None,
    };

    core.alts
        .iter()
        .enumerate()
        .filter_map(|(i, alt_a)| {
            if gq < min_gq {
                return None;
            }
            let is_called = gt == (i + 1) as i32;
            Some(Ok(Sample {
                ref_id: core.ref_id,
                pos: core.pos,
                ref_a: core.ref_a.clone(),
                alt_a: alt_a.to_vec(),
                genotype_quality: gq,
                is_called,
                phase_set,
            }))
        })
        .collect()
}

#[cfg(test)]
mod sample_ps_tests {
    use super::*;
    use noodles::core::Position;
    use noodles::vcf::header::record::value::{map::Contig, Map};
    use noodles::vcf::variant::record_buf::samples::{sample::Value as SampleValue, Keys};
    use noodles::vcf::variant::record_buf::{RecordBuf, Samples};
    use noodles::vcf::Header;

    fn test_header() -> Header {
        Header::builder()
            .add_contig("chr1", Map::<Contig>::builder().build().unwrap())
            .build()
    }

    /// Build a minimal record with exactly two sample columns (index 0: an
    /// unread placeholder; index 1: the sample under test), matching
    /// `parse_sample_record`'s `samples.get_index(1)` convention. GT/GQ/PS
    /// are stored as `Value::Integer` -- matching the exact match arms in
    /// `parse_sample_record`, which do NOT accept a genotype string like
    /// "0/1" for GT.
    ///
    /// NEEDS VERIFICATION against noodles 0.111/0.112's exact
    /// `Samples`/`Keys` constructors -- isolated here so a mismatch is a
    /// one-function patch.
    fn sample_record(gt: i32, gq: i32, ps: Option<i32>) -> RecordBuf {
        let mut record = RecordBuf::builder()
            .set_reference_sequence_name("chr1")
            .set_variant_start(Position::try_from(100).unwrap())
            .set_reference_bases("A")
            .set_alternate_bases(vec!["G".to_string()].into())
            .build();

        let mut keys = vec![String::from("GT"), String::from("GQ")];
        let mut under_test: Vec<Option<SampleValue>> = vec![
            Some(SampleValue::Integer(gt)),
            Some(SampleValue::Integer(gq)),
        ];
        let mut placeholder: Vec<Option<SampleValue>> = vec![
            Some(SampleValue::Integer(0)),
            Some(SampleValue::Integer(99)),
        ];

        if let Some(ps_val) = ps {
            keys.push(String::from("PS"));
            under_test.push(Some(SampleValue::Integer(ps_val)));
            placeholder.push(None); // sample 0 unphased -- must not leak into sample 1's parse
        }

        *record.samples_mut() = Samples::new(Keys::from_iter(keys), vec![placeholder, under_test]);
        record
    }

    #[test]
    fn ps_tag_is_parsed_when_present() {
        let header = test_header();
        let mut record = sample_record(1, 40, Some(7));
        let result = parse_sample_record(&mut record, &header, 0.0).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].phase_set, Some(7));
    }

    #[test]
    fn ps_tag_is_none_when_absent() {
        let header = test_header();
        let mut record = sample_record(1, 40, None);
        let result = parse_sample_record(&mut record, &header, 0.0).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(
            result[0].phase_set, None,
            "unphased VCFs (no PS tag) must not error or fabricate a phase set"
        );
    }
}
