use crate::bam_format::{path_unicode_ok, out_from_file, out_stdout};
use crate::variant::{
    PopulationVariant, SampleVariant, parse_population_record, parse_sample_record,
};
use crate::vcf_format::vcf_reader;
use crate::{Config, StripReadSuffix};
use anyhow::{Result, anyhow, ensure};
use rust_htslib::bam::{self, Read, record::Record};
use std::collections::HashMap;

pub trait AlignmentStream {
    fn next_qname(&self) -> &[u8];
    fn un_next(&mut self, rec: Record) -> Result<()>;
    fn next_rec(&mut self) -> Result<Option<Record>>;
    fn write_record(&mut self, rec: Record, is_best: Option<bool>) -> Result<()>;
    fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<()>;
}

pub struct AlnStream {
    ambiguous: Option<bam::Writer>,
    pub bam: Option<bam::Reader>,
    filt: Option<bam::Writer>,
    next: Option<Record>,
    output: Option<bam::Writer>,
    #[allow(dead_code)]
    sample_variants: Option<Vec<HashMap<i64, Vec<SampleVariant>>>>,
    #[allow(dead_code)]
    population_variants: Option<Vec<HashMap<i64, Vec<PopulationVariant>>>>,
}

impl AlnStream {
    pub fn new(opt: &mut Config, i: usize) -> Result<Self> {
        let bam_str = opt.alignment[i].as_str();
        let mut bam = bam::Reader::from_path(bam_str)?;

        let mut test_record = Record::new();
        bam.read(&mut test_record)
            .transpose()?
            .ok_or_else(|| anyhow!("{bam_str} has no records"))?;

        let qname = test_record.qname();
        let strip_read_suffix = match opt.strip_read_suffix {
            StripReadSuffix::True => {
                ensure!(
                    qname.ends_with(b"/1") || qname.ends_with(b"/2"),
                    "Input read names do not have /1 or /2 suffixes, but strip_read_suffix is true."
                );
                StripReadSuffix::True
            }
            StripReadSuffix::False => {
                ensure!(
                    !qname.ends_with(b"/1") && !qname.ends_with(b"/2"),
                    "Input read names have /1 or /2 suffixes, but strip_read_suffix is false."
                );
                StripReadSuffix::False
            }
            StripReadSuffix::Auto => {
                // Auto-detect based on first read
                if qname.ends_with(b"/1") || qname.ends_with(b"/2") {
                    StripReadSuffix::True
                } else {
                    StripReadSuffix::False
                }
            }
            StripReadSuffix::Variable => StripReadSuffix::Variable,
        };
        let is_paired = if i == 0 && opt.is_paired.is_none() {
            Some(test_record.is_paired())
        } else {
            ensure!(
                opt.is_paired == Some(test_record.is_paired()),
                "All input BAMs must be either paired-end or single-end."
            );
            opt.is_paired
        };

        let next = Some(test_record);
        let sample_variants = opt
            .sample_variants
            .get(i)
            .map(|_| vcf_reader(&opt.sample_variants[i], parse_sample_record))
            .transpose()?;
        let population_variants = opt
            .population_variants
            .get(i)
            .map(|_| vcf_reader(&opt.population_variants[i], parse_population_record))
            .transpose()?;

        // check output paths are unicode here, so we hopefully only create files once all are ok.
        opt.output
            .get(i)
            .map(|f| path_unicode_ok(f))
            .transpose()?;
        opt.filtered_output
            .get(i)
            .map(|f| path_unicode_ok(f))
            .transpose()?;
        opt.ambiguous_output
            .get(i)
            .map(|f| path_unicode_ok(f))
            .transpose()?;

        let stream = AlnStream {
            ambiguous: None,
            bam: Some(bam),
            filt: None,
            next,
            output: None,
            sample_variants,
            population_variants,
        };
        let bam = stream.bam.as_ref().expect("no in");

        let parts: Vec<Vec<&[u8]>> = bam
            .header()
            .as_bytes()
            .split(|&b| b == b'\n')
            .map(|s| s.split(|&b| b == b'\t').collect())
            .collect();

        // Should be first line in header, but in some version does not comply to specs.
        for p in &parts {
            ensure!(
                p.len() < 3
                    || p[0] != b"@HD"
                    || (p[2] != b"SO:coordinate" && p[2] != b"GO:reference"),
                "Coordinate sorted input, would require hashmap lookup."
            );
        }

        Ok(stream)
    }
}

impl AlignmentStream for AlnStream {
    fn next_qname(&self) -> &[u8] {
        self.next.as_ref().map_or(b"", |r| r.qname())
    }

    fn un_next(&mut self, rec: Record) -> Result<()> {
        if self.next.is_some() {
            return Err(anyhow!("Cannot un-next more than one record"));
        }
        self.next = Some(rec);
        Ok(())
    }

    fn next_rec(&mut self) -> Result<Option<Record>> {
        self.next
            .take()
            .map(Ok)
            .or_else(|| {
                let mut record = Record::new();
                match self.bam.as_mut().and_then(|b| b.read(&mut record)) {
                    Some(Err(e)) => Some(Err(anyhow!(e))),
                    Some(Ok(())) => Some(Ok(record)),
                    None => None,
                }
            })
            .transpose()
    }

    fn write_record(&mut self, rec: Record, is_best: Option<bool>) -> Result<()> {
        let output = match is_best {
            Some(true) => self.output.as_mut(),
            Some(false) => self.filt.as_mut(),
            None => self.ambiguous.as_mut(),
        };
        output.map(|output| output.write(&rec)).transpose()?;
        Ok(())
    }
    fn init_writers(&mut self, opt: &Config, i: usize) -> Result<()> {

        #[cfg(test)]
        let allow_stdout = false;
        #[cfg(not(test))]
        let allow_stdout = true;

        let bam = self.bam.as_ref().ok_or_else(|| anyhow!("No BAM reader"))?;
        self.output = opt
            .output
            .get(i)
            .map(|f| out_from_file(f, bam.header()))
            .transpose()?;

        if i == 0 && self.output.is_none() && allow_stdout {
            self.output = Some(out_stdout(bam.header(), opt.stdout_format.into())?);
        }
        self.filt = opt
            .filtered_output
            .get(i)
            .map(|f| out_from_file(f, bam.header()))
            .transpose()?;
        self.ambiguous = opt
            .ambiguous_output
            .get(i)
            .map(|f| out_from_file(f, bam.header()))
            .transpose()?;
        Ok(())
    }
}

// Minimal Mock to satisfy the AlnStream trait/interface

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::StripReadSuffix;
    use crate::tests::{BamFormat, create_record};

    pub struct MockStream {
        pub reads: Vec<Record>,
        written: Vec<(Record, Option<bool>)>,
        aln_stream: AlnStream,
        i: usize,
    }

    impl AlignmentStream for MockStream {
        fn next_qname(&self) -> &[u8] {
            self.aln_stream.next_qname()
        }

        fn un_next(&mut self, rec: Record) -> Result<()> {
            self.un_next(rec)
        }

        fn next_rec(&mut self) -> Result<Option<Record>> {
            self.next_rec()
        }

        fn write_record(&mut self, rec: Record, is_best: Option<bool>) -> Result<()> {
            self.write_record(rec, is_best)
        }
        fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<()> {
            Ok(())
        }
    }

    impl MockStream {
        pub fn new(i: usize, reads: Vec<Record>) -> Self {
            let aln_stream = AlnStream {
                ambiguous: None,
                bam: None,
                filt: None,
                next: None,
                output: None,
                sample_variants: None,
                population_variants: None,
            };
            Self {
                reads,
                written: Vec::new(),
                aln_stream,
                i,
            }
        }
        fn next_rec(&mut self) -> Result<Option<Record>> {
            match self.aln_stream.next_rec()? {
                Some(rec) => {
                    /*eprintln!(
                        "re-nexted({}): {}",
                        self.i,
                        std::str::from_utf8(rec.qname()).unwrap_or("Invalid UTF-8")
                    );*/
                    return Ok(Some(rec));
                }
                None => {}
            }
            if self.reads.is_empty() {
                return Ok(None);
            }
            let rec = self.reads.remove(0);
            self.aln_stream.un_next(rec)?;
            self.aln_stream.next_rec()
        }
        fn un_next(&mut self, rec: Record) -> Result<()> {
            /*eprintln!(
                "Un-next({}) read: {}",
                self.i,
                std::str::from_utf8(rec.qname()).unwrap_or("Invalid UTF-8")
            );*/
            self.aln_stream.un_next(rec)
        }
        fn write_record(&mut self, rec: Record, state: Option<bool>) -> Result<()> {
            self.written.push((rec, state));
            Ok(())
        }
    }

    /*#[test]
    fn test_aln_stream_new_ok() -> Result<()> {
        let mut config = Config {
            alignment: vec!["tests/data/test_input_1_a.bam".to_string()],
            stdout_format: BamFormat::Sam,
            ..Default::default()
        };

        let aln_stream = AlnStream::new(&mut config, 0)?;
        let qname = std::str::from_utf8(aln_stream.next_qname())?;
        assert_eq!(qname, "r000");
        Ok(())
    }*/

    #[test]
    fn test_aln_stream_new_mismatch_strip_suffix_true_instead_of_false() {
        let mut config = Config {
            alignment: vec!["tests/data/test_input_1_a.bam".to_string()],
            stdout_format: BamFormat::Sam,
            strip_read_suffix: StripReadSuffix::True,
            ..Default::default()
        };

        let result = AlnStream::new(&mut config, 0);
        assert!(result.is_err());
    }
    #[test]
    fn test_aln_stream_next_rec() -> Result<()> {
        let records = vec![
            create_record(b"read1/1", "10M", &[], &vec![30; 10], "10", false)?,
            create_record(b"read2/1", "10M", &[], &vec![30; 10], "10", false)?,
        ];
        let mut mock_stream = MockStream::new(0, records.clone());

        for expected in records {
            let rec = mock_stream.next_rec()?.unwrap();
            assert_eq!(rec.qname(), expected.qname());
        }
        assert!(mock_stream.next_rec()?.is_none());
        Ok(())
    }
    #[test]
    fn test_aln_stream_un_next() -> Result<()> {
        let records = vec![
            create_record(b"read1/1", "10M", &[], &vec![30; 10], "10", false)?,
            create_record(b"read2/1", "10M", &[], &vec![30; 10], "10", false)?,
        ];
        let mut mock_stream = MockStream::new(0, records);

        let rec1 = mock_stream.next_rec()?.unwrap();
        assert_eq!(rec1.qname(), b"read1/1");

        mock_stream.un_next(rec1)?;

        let rec2 = mock_stream.next_rec()?.unwrap();
        assert_eq!(rec2.qname(), b"read1/1");

        let rec3 = mock_stream.next_rec()?.unwrap();
        assert_eq!(rec3.qname(), b"read2/1");

        Ok(())
    }
    /*#[test]
    fn test_aln_stream_new_force_strip_suffix_ok() -> Result<()> {
        // None of the test BAMs have read suffixes..
        let mut config = Config {
            alignment: vec!["tests/data/test_input_1_b.bam".to_string()],
            stdout_format: BamFormat::Sam,
            strip_read_suffix: StripReadSuffix::True,
            ..Default::default()
        };

        let aln_stream = AlnStream::new(&mut config, 0)?;
        //assert!(aln_stream.bam.is_paired());
        assert_eq!(config.strip_read_suffix, StripReadSuffix::True);
        Ok(())
    }
    #[test]
    fn test_aln_stream_new_auto_detect_true_ok() -> Result<()> {
        // None of the test BAMs have read suffixes..
        let mut config = Config {
            alignment: vec!["tests/data/test_input_1_b.bam".to_string()],
            stdout_format: BamFormat::Sam,
            strip_read_suffix: StripReadSuffix::Auto,
            ..Default::default()
        };

        let aln_stream = AlnStream::new(&mut config, 0)?;
        //assert!(aln_stream.bam.is_paired());
        assert_eq!(config.strip_read_suffix, StripReadSuffix::True);
        Ok(())
    }
    fn test_aln_stream_new_mismatch_strip_suffix_false_instead_of_true() {
        // None of the test BAMs have read suffixes..
        let mut config = Config {
            alignment: vec!["tests/data/test_input_1_b.bam".to_string()],
            stdout_format: BamFormat::Sam,
            strip_read_suffix: StripReadSuffix::False,
            ..Default::default()
        };

        let result = AlnStream::new(&mut config, 0);
        assert!(result.is_err());
    }*/
}
