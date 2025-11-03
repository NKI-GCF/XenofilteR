use std::collections::HashMap;

use anyhow::{Result, anyhow, ensure};
use rust_htslib::bam::{self, Read, record::Record};

use crate::Config;
use crate::bam_format::{out_from_file, out_stdout};
use crate::vcf_format::vcf_reader;
use crate::variant::{parse_sample_record, parse_population_record, SampleVariant, PopulationVariant};

pub struct AlnStream {
    ambiguous: Option<bam::Writer>,
    pub bam: bam::Reader,
    filt: Option<bam::Writer>,
    next: Option<Record>,
    output: Option<bam::Writer>,
    sample_variants: Option<Vec<HashMap<i64, Vec<SampleVariant>>>>,
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
        match opt.strip_read_suffix {
            Some(true) => {
                ensure!(
                    qname.ends_with(b"/1") || qname.ends_with(b"/2"),
                    "Input read names do not have /1 or /2 suffixes, but strip_read_suffix is true."
                );
            }
            Some(false) => {
                ensure!(
                    !qname.ends_with(b"/1") && !qname.ends_with(b"/2"),
                    "Input read names have /1 or /2 suffixes, but strip_read_suffix is false."
                );
            }
            None => {
                // auto-detect
                if qname.ends_with(b"/1") || qname.ends_with(b"/2") {
                    opt.strip_read_suffix = Some(true);
                } else {
                    opt.strip_read_suffix = Some(false);
                }
            }
        }

        let next = Some(test_record);
        let output = opt
            .output
            .get(i)
            .map(|f| out_from_file(f, bam.header()))
            .transpose()?;
        let filt = opt
            .filtered_output
            .get(i)
            .map(|f| out_from_file(f, bam.header()))
            .transpose()?;
        let ambiguous = opt
            .ambiguous_output
            .get(i)
            .map(|f| out_from_file(f, bam.header()))
            .transpose()?;
        let sample_variants = opt.sample_variants.get(i).map(|_| {
            vcf_reader(
                &opt.sample_variants[i],
                parse_sample_record,
            )
        }).transpose()?;
        let population_variants = opt.population_variants.get(i).map(|_| {
            vcf_reader(
                &opt.population_variants[i],
                parse_population_record,
            )
        }).transpose()?;

        let mut stream = AlnStream {
            ambiguous,
            bam,
            filt,
            next,
            output,
            sample_variants,
            population_variants,
        };

        if i == 0 && stream.output.is_none() {
            stream.output = Some(out_stdout(stream.bam.header(), opt.stdout_format.into())?);
        }
        let parts: Vec<Vec<&[u8]>> = stream
            .bam
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
    pub fn next_qname(&self) -> &[u8] {
        self.next.as_ref().map_or(b"", |r| r.qname())
    }

    pub fn un_next(&mut self, rec: Record) {
        self.next = Some(rec);
    }

    pub fn next_rec(&mut self) -> Result<Option<Record>> {
        self.next
            .take()
            .map(Ok)
            .or_else(|| {
                let mut record = Record::new();
                match self.bam.read(&mut record) {
                    Some(Err(e)) => Some(Err(anyhow!(e))),
                    Some(Ok(())) => Some(Ok(record)),
                    None => None,
                }
            })
            .transpose()
    }

    pub fn write_record(&mut self, rec: Record, is_best: Option<bool>) -> Result<()> {
        let output = match is_best {
            Some(true) => self.output.as_mut(),
            Some(false) => self.filt.as_mut(),
            None => self.ambiguous.as_mut(),
        };
        output.map(|output| output.write(&rec)).transpose()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn skip_filtered_write_first_bit_differs() {
        // in progress_state() this is used to skip unmapped, if in config
        assert_eq!(SKIP_FILTERED | 1, WRITE_FILTERED);
    }
}
