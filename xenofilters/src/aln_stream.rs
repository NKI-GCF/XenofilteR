//! Alignment stream abstraction.
//!
//! [`AlnStream`] wraps a name-sorted BAM reader together with its optional
//! output writers and variant stores. The [`AlignmentStream`] trait is the
//! interface used by the filter algorithm so that tests can substitute a
//! lightweight mock without touching any file I/O.

use crate::alignment::SimpleRec;
use crate::bam::{out_from_file, path_unicode_ok, BamOutput};
use crate::config::{Config, StripReadSuffix};
use crate::variant::{
    parse_population_record, parse_sample_record, Population, Sample, Store, StoreTrait,
};
use anyhow::{anyhow, ensure, Result};
use noodles::bam::{io::Reader as BamReader, record::Record};
use noodles::bgzf::io::Reader as BgzfReader;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::header::Header;
use std::fs::File;
use std::io::Read as ioRead;
use std::path::PathBuf;

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
    ambiguous: Option<BamOutput>,
    pub(crate) bam: Option<BamReader<BgzfReader<File>>>,
    filt: Option<BamOutput>,
    /// One-record look-ahead buffer (supports `un_next`).
    next: Option<R>,
    output: Option<BamOutput>,
    sample_variants: Option<Store<Sample>>,
    population_variants: Option<Store<Population>>,
    header: Header,
    /// Cached from config so writers can be opened later with the correct count.
    threads: usize,
}

impl AlnStream<Record> {
    pub(crate) fn new(opt: &mut Config, i: usize) -> Result<Self> {
        let bam_str = opt.alignment[i].as_str();
        tracing::debug!(stream = i, path = bam_str, "Opening BAM reader");

        let file = File::open(bam_str)
            .map_err(|e| anyhow!("Cannot open alignment file '{bam_str}': {e}"))?;
        let mut bam = BamReader::new(file);
        let header = bam.read_header()?;

        // -- Reject coordinate-sorted input --------------------------------
        // Coordinate-sorted files would require a hash map of all in-flight
        // fragments, defeating the low-memory streaming design.
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
                "Input file '{bam_str}' appears to be coordinate-sorted \
                 (SO:coordinate / GO:reference). xenofilters requires \
                 name-sorted BAM files. Re-sort with `samtools sort -n`."
            );
        }

        // -- Peek at the first record --------------------------------------
        // We need it to (a) detect whether reads carry /1 /2 suffixes and
        // (b) synchronise the name check in main.
        let test_record = match bam.records().next() {
            Some(Ok(rec)) => rec,
            Some(Err(e)) => return Err(anyhow!("Error reading from '{bam_str}': {e}")),
            None => return Err(anyhow!("'{bam_str}' contains no records")),
        };

        let name = test_record
            .name()
            .ok_or_else(|| anyhow!("First record in '{bam_str}' has no read name"))?;

        // -- Auto-detect / validate read-name suffix stripping -------------
        opt.strip_read_suffix = match opt.strip_read_suffix {
            StripReadSuffix::True => {
                ensure!(
                    name.ends_with(b"/1") || name.ends_with(b"/2"),
                    "Stream {i} ('{bam_str}'): --strip-read-suffix=true requested but \
                     first read name '{name:?}' has no /1 or /2 suffix."
                );
                StripReadSuffix::True
            }
            StripReadSuffix::False => {
                ensure!(
                    !name.ends_with(b"/1") && !name.ends_with(b"/2"),
                    "Stream {i} ('{bam_str}'): --strip-read-suffix=false requested but \
                     first read name '{name:?}' ends with a /1 or /2 suffix."
                );
                StripReadSuffix::False
            }
            StripReadSuffix::Auto => {
                // Auto-detect based on first read
                if name.ends_with(b"/1") || name.ends_with(b"/2") {
                    tracing::debug!(stream = i, "Auto-detected /1 /2 read-name suffixes");
                    StripReadSuffix::True
                } else {
                    StripReadSuffix::False
                }
            }
            StripReadSuffix::Variable => StripReadSuffix::Variable,
        };

        // -- Paired-end consistency check ----------------------------------
        opt.is_paired = if i == 0 && opt.is_paired.is_none() {
            let paired = test_record.flags().is_segmented();
            tracing::debug!(paired, "Auto-detected paired-end status from stream 0");
            Some(paired)
        } else {
            ensure!(
                opt.is_paired == Some(test_record.flags().is_segmented()),
                "Stream {i} ('{bam_str}') has different paired-end status than stream 0. \
                 All inputs must be either all paired-end or all single-end."
            );
            opt.is_paired
        };

        // -- Variant stores ------------------------------------------------
        let sample_variants = opt
            .sample_variants
            .get(i)
            .filter(|s| !s.is_empty())
            .map(|s| {
                tracing::debug!(stream = i, path = s, "Loading sample variants");
                Store::new(&PathBuf::try_from(s.as_str())?, parse_sample_record)
            })
            .transpose()?;

        let population_variants = opt
            .population_variants
            .get(i)
            .filter(|s| !s.is_empty())
            .map(|s| {
                tracing::debug!(stream = i, path = s, "Loading population variants");
                Store::new(&PathBuf::try_from(s.as_str())?, parse_population_record)
            })
            .transpose()?;

        // -- Validate output paths early -----------------------------------
        // Check before opening writers so no partial files are created if one
        // path is invalid.
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
            header,
            threads: opt.threads,
        })
    }
}

trait FromBamRecord: Sized {
    fn from_bam_record(header: &Header, rec: Record) -> std::io::Result<Self>;
}

impl FromBamRecord for Record {
    fn from_bam_record(_header: &Header, rec: Record) -> std::io::Result<Self> {
        Ok(rec)
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
        self.next
            .as_ref()
            .and_then(|r| r.name())
            .map_or(b"", |n| n.as_ref())
    }

    fn un_next(&mut self, rec: R) -> Result<()> {
        if self.next.is_some() {
            return Err(anyhow!(
                "un_next called while look-ahead buffer is already occupied \
                 (this is a bug — please report it)"
            ));
        }
        self.next = Some(rec);
        Ok(())
    }

    fn next_rec(&mut self) -> Result<Option<R>> {
        // Return the buffered record first, then read from the BAM reader.
        let header = &self.header;
        Ok(self
            .next
            .take()
            .map(Ok)
            .or_else(|| {
                self.bam.as_mut().and_then(|b| {
                    b.records()
                        .next()
                        .map(|r| r.and_then(|r| R::from_bam_record(header, r)))
                })
            })
            .transpose()?)
    }

    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<()> {
        let output: Option<&mut BamOutput> = match is_best {
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
        let add_pg = !opt.no_program_line;
        let threads = self.threads;

        self.output = opt
            .output
            .get(i)
            .map(|f| out_from_file(f, &self.header, add_pg, threads))
            .transpose()?;
        self.filt = opt
            .filtered_output
            .get(i)
            .map(|f| out_from_file(f, &self.header, add_pg, threads))
            .transpose()?;
        self.ambiguous = opt
            .ambiguous_output
            .get(i)
            .map(|f| out_from_file(f, &self.header, add_pg, threads))
            .transpose()?;

        tracing::debug!(
            stream = i,
            output = opt.output.get(i).map(|p| p.display().to_string()),
            filtered = opt.filtered_output.get(i).map(|p| p.display().to_string()),
            ambiguous = opt.ambiguous_output.get(i).map(|p| p.display().to_string()),
            "Writers initialised"
        );

        Ok(())
    }

    fn variant_store(&self) -> Option<&dyn StoreTrait> {
        self.sample_variants
            .as_ref()
            .map(|s| s as &dyn StoreTrait)
            .or_else(|| {
                self.population_variants
                    .as_ref()
                    .map(|p| p as &dyn StoreTrait)
            })
    }

    fn header(&self) -> &Header {
        &self.header
    }
}

#[cfg(test)]
pub(crate) mod tests;
