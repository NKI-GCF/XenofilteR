use crate::bam::{out_from_file, path_unicode_ok};
use crate::variant::{
    StoreTrait, Store, Population, Sample, parse_population_record, parse_sample_record
};
use crate::config::{Config, StripReadSuffix};
use anyhow::{Result, anyhow, ensure};
use noodles::bam::{io::{Reader as BamReader, Writer as BamWriter}, record::Record};
use std::io::Read as ioRead;
use noodles::bgzf::io::{Reader as BgzfReader, Writer as BgzfWriter};
use std::fs::File;
use noodles::sam::header::Header;
use noodles::sam::alignment::io::Write;
use crate::alignment::SimpleRec;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::io::{Error, ErrorKind};


pub(crate) trait AlignmentStream<R: SimpleRec> {
    fn next_qname(&self) -> &[u8];
    fn un_next(&mut self, rec: R) -> Result<()>;
    fn next_rec(&mut self) -> Result<Option<R>>;
    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<()>;
    fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<()>;
    fn variant_store(&self) -> Option<&dyn StoreTrait>;
    fn header(&self) -> &Header;
}

pub(crate) struct AlnStream<R> {
    ambiguous: Option<BamWriter<BgzfWriter<File>>>,
    pub(crate) bam: Option<BamReader<BgzfReader<File>>>,
    filt: Option<BamWriter<BgzfWriter<File>>>,
    next: Option<R>,
    output: Option<BamWriter<BgzfWriter<File>>>,
    sample_variants: Option<Store<Sample>>,
    population_variants: Option<Store<Population>>,
    header: Header,
}

impl AlnStream<Record> {
    pub(crate) fn new(opt: &mut Config, i: usize) -> Result<Self> {
        let bam_str = opt.alignment[i].as_str();
        let file = File::open(bam_str)?;
        let mut bam = BamReader::new(file);
        let header = bam.read_header()?;

        // XXX: workaround for sort order.
        let mut header_reader = bam.header_reader();
        header_reader.read_magic_number()?;
        let mut raw_sam_header_reader = header_reader.raw_sam_header_reader()?;

        let mut buf = Vec::new();
        raw_sam_header_reader.read_to_end(&mut buf)?;

        for parts in buf
            .split(|&b| b == b'\n')
            .map(|s| s.split(|&b| b == b'\t').collect::<Vec<_>>())
        {
            ensure!(
                parts.len() < 3
                    || parts[0] != b"@HD"
                    || (parts[2] != b"SO:coordinate" && parts[2] != b"GO:reference"),
                "Coordinate sorted input, would require hashmap lookup."
            );
        }

        let test_record = match bam.records().next() {
            Some(Ok(rec)) => rec,
            Some(Err(e)) => return Err(anyhow!(e)),
            None => return Err(anyhow!("{bam_str} has no records")),
        };

        let name = test_record.name().ok_or_else(|| anyhow!("Record has no name"))?;
        opt.strip_read_suffix = match opt.strip_read_suffix {
            StripReadSuffix::True => {
                ensure!(
                    name.ends_with(b"/1") || name.ends_with(b"/2"),
                    "Input read names do not have /1 or /2 suffixes, but strip_read_suffix is true."
                );
                StripReadSuffix::True
            }
            StripReadSuffix::False => {
                ensure!(
                    !name.ends_with(b"/1") && !name.ends_with(b"/2"),
                    "Input read names have /1 or /2 suffixes, but strip_read_suffix is false."
                );
                StripReadSuffix::False
            }
            StripReadSuffix::Auto => {
                // Auto-detect based on first read
                if name.ends_with(b"/1") || name.ends_with(b"/2") {
                    StripReadSuffix::True
                } else {
                    StripReadSuffix::False
                }
            }
            StripReadSuffix::Variable => StripReadSuffix::Variable,
        };
        opt.is_paired = if i == 0 && opt.is_paired.is_none() {
            Some(test_record.flags().is_segmented())
        } else {
            ensure!(
                opt.is_paired == Some(test_record.flags().is_segmented()),
                "All input BAMs must be either paired-end or single-end."
            );
            opt.is_paired
        };

        let sample_variants = opt
            .sample_variants
            .get(i)
            .map(|p| Store::new(p, parse_sample_record))
            .transpose()?;
        let population_variants = opt
            .population_variants
            .get(i)
            .map(|p| Store::new(p, parse_population_record))
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

        Ok(AlnStream {
            ambiguous: None,
            bam: Some(bam),
            filt: None,
            next: Some(test_record),
            output: None,
            sample_variants,
            population_variants,
            header
        })
    }
}

trait FromBamRecord: Sized {
    fn from_bam_record(header: &Header, rec: Record) -> std::io::Result<Self>;
}

impl FromBamRecord for Record {
    fn from_bam_record(_header: &Header, rec: Record) -> std::io::Result<Self> {
        Record::try_from(rec).map_err(|_| Error::new(ErrorKind::Other, "next_rec"))
    }
}

impl FromBamRecord for RecordBuf {
    fn from_bam_record(header: &Header, rec: Record) -> std::io::Result<Self> {
        RecordBuf::try_from_alignment_record(header, &rec)
    }
}

impl<R> AlignmentStream<R> for AlnStream<R>
where
    R: SimpleRec + FromBamRecord,
{
    fn next_qname(&self) -> &[u8] {
        self.next.as_ref().and_then(|r| r.name()).map_or(b"", |n| n.as_ref())
    }

    fn un_next(&mut self, rec: R) -> Result<()> {
        if self.next.is_some() {
            return Err(anyhow!("Cannot un-next more than one record"));
        }
        self.next = Some(rec);
        Ok(())
    }

    fn next_rec(&mut self) -> Result<Option<R>> {
        let header = &self.header;   // disjoint field borrow before &mut self.bam
        Ok(self.next
            .take()
            .map(Ok)
            .or_else(|| {
                self.bam.as_mut().and_then(|b| {
                    b.records().next().map(|r| {
                        r.and_then(|r| R::from_bam_record(header, r))
                    })
                })
            })
            .transpose()?)
    }

    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<()> {
        let output = match is_best {
            Some(true) => self.output.as_mut(),
            Some(false) => self.filt.as_mut(),
            None => self.ambiguous.as_mut(),
        };
        if let Some(o) = output {
            o.write_alignment_record(&self.header, &rec)?;
        }
        Ok(())
    }

    fn init_writers(&mut self, opt: &Config, i: usize) -> Result<()> {
        let add_pg_line = !opt.no_program_line;
        self.output = opt.output.get(i)
            .map(|f| out_from_file(f, &self.header, add_pg_line))
            .transpose()?;
        self.filt = opt.filtered_output.get(i)
            .map(|f| out_from_file(f, &self.header, add_pg_line))
            .transpose()?;
        self.ambiguous = opt.ambiguous_output.get(i)
            .map(|f| out_from_file(f, &self.header, add_pg_line))
            .transpose()?;
        Ok(())
    }

    fn variant_store(&self) -> Option<&dyn StoreTrait> {
        self.sample_variants.as_ref().map(|s| s as &dyn StoreTrait)
            .or_else(|| self.population_variants.as_ref().map(|p| p as &dyn StoreTrait))
    }

    fn header(&self) -> &Header {
        &self.header
    }
}

#[cfg(test)]
pub(crate) mod tests;
