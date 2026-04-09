use crate::bam_format::{out_from_file, out_stdout, path_unicode_ok};
use crate::variant::{
    VariantStoreTrait, VariantStore, PopulationVariant, SampleVariant, parse_population_record, parse_sample_record, vcf_reader
};
use crate::{Config, StripReadSuffix};
use anyhow::{Result, anyhow, ensure};
use rust_htslib::bam::{self, Read, record::Record};

pub trait AlignmentStream {
    fn next_qname(&self) -> &[u8];
    fn un_next(&mut self, rec: Record) -> Result<()>;
    fn next_rec(&mut self) -> Result<Option<Record>>;
    fn write_record(&mut self, rec: Record, is_best: Option<bool>) -> Result<()>;
    fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<()>;
    fn variant_store(&self) -> Option<&dyn VariantStoreTrait>;
}

pub struct AlnStream {
    ambiguous: Option<bam::Writer>,
    pub bam: Option<bam::Reader>,
    filt: Option<bam::Writer>,
    next: Option<Record>,
    output: Option<bam::Writer>,
    sample_variants: Option<VariantStore<SampleVariant>>,
    population_variants: Option<VariantStore<PopulationVariant>>,
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
        opt.strip_read_suffix = match opt.strip_read_suffix {
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
        opt.is_paired = if i == 0 && opt.is_paired.is_none() {
            Some(test_record.is_paired())
        } else {
            ensure!(
                opt.is_paired == Some(test_record.is_paired()),
                "All input BAMs must be either paired-end or single-end."
            );
            opt.is_paired
        };

        let sample_variants = opt
            .sample_variants
            .get(i)
            .map(|p| vcf_reader(p, parse_sample_record))
            .transpose()?;
        let population_variants = opt
            .population_variants
            .get(i)
            .map(|p| vcf_reader(p, parse_population_record))
            .transpose()?;

        // check output paths are unicode here, so we hopefully only create files once all are ok.
        opt.output.get(i).map(path_unicode_ok).transpose()?;
        opt.filtered_output
            .get(i)
            .map(path_unicode_ok)
            .transpose()?;
        opt.ambiguous_output
            .get(i)
            .map(path_unicode_ok)
            .transpose()?;

        let stream = AlnStream {
            ambiguous: None,
            bam: Some(bam),
            filt: None,
            next: Some(test_record),
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
        let add_pg_line = !opt.no_program_line;

        let bam = self.bam.as_ref().ok_or_else(|| anyhow!("No BAM reader"))?;
        self.output = opt
            .output
            .get(i)
            .map(|f| out_from_file(f, bam.header(), add_pg_line))
            .transpose()?;

        if i == 0 && self.output.is_none() && allow_stdout {
            self.output = Some(out_stdout(
                bam.header(),
                opt.stdout_format.into(),
                add_pg_line,
            )?);
        }
        self.filt = opt
            .filtered_output
            .get(i)
            .map(|f| out_from_file(f, bam.header(), add_pg_line))
            .transpose()?;
        self.ambiguous = opt
            .ambiguous_output
            .get(i)
            .map(|f| out_from_file(f, bam.header(), add_pg_line))
            .transpose()?;
        Ok(())
    }
    fn variant_store(&self) -> Option<&dyn VariantStoreTrait> {
        self.sample_variants
            .as_ref()
            .map(|s| s as &dyn VariantStoreTrait)
            .or_else(|| self.population_variants.as_ref().map(|p| p as &dyn VariantStoreTrait))
    }
}

// Minimal Mock to satisfy the AlnStream trait/interface

#[cfg(test)]
pub mod tests;
