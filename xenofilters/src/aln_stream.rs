use crate::alignment::SimpleRec;
use crate::bam::{out_from_file, path_unicode_ok, BamOutput, OutputMode};
use crate::config::{Config, StripReadSuffix};
use crate::variant::{
    parse_population_record, parse_sample_record, Population, Sample, Store, StoreTrait,
};
use anyhow::{anyhow, ensure, Result};
use noodles::bam::{
    io::{Reader as BamReader},
    record::Record,
};
use noodles::bgzf::io::{Reader as BgzfReader};
use noodles::bgzf::VirtualPosition;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::Header;
use std::fs::File;
use std::io::{Read as ioRead, Seek as ioSeek};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Arc;

pub(crate) trait AlignmentStream<R: SimpleRec> {
    fn next_qname(&self) -> &[u8];
    fn un_next(&mut self, rec: R) -> Result<()>;
    fn next_rec(&mut self) -> Result<Option<R>>;
    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<()>;
    fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<()>;
    fn variant_store(&self) -> Option<Arc<dyn StoreTrait>>;
    fn header(&self) -> &Header;
    /// Seek to a BGZF virtual offset and read one full record for pass-2 output.
    /// Returns `Err` for stream types that do not support seeking (e.g. mock streams).
    fn fetch_by_virtual_offset(&mut self, _virtual_offset: u64) -> Result<RecordBuf> {
        Err(anyhow!(
            "fetch_by_virtual_offset not supported for this stream type"
        ))
    }
}

pub(crate) trait BamStreamReader {
    fn next_record(&mut self) -> Option<std::io::Result<Record>>;
    fn seek_vpos(&mut self, pos: VirtualPosition) -> std::io::Result<VirtualPosition>;
}

impl<T: ioRead + ioSeek> BamStreamReader for BamReader<BgzfReader<T>> {
    fn next_record(&mut self) -> Option<std::io::Result<Record>> {
        self.records().next()
    }

    fn seek_vpos(&mut self, pos: VirtualPosition) -> std::io::Result<VirtualPosition> {
        // Bypass the BAM wrapper and seek directly on the underlying BGZF reader
        self.get_mut().seek(pos)
        /* Depending on the exact micro-version/patch-level resolution of your noodles crates,
         * BgzfReader might expect standard std::io::SeekFrom positioning where the packed u64
         * is treated as the virtual position. If the block above throws a type mismatch error,
         * use this variant instead:
        self.get_mut()
            .seek(std::io::SeekFrom::Start(u64::from(pos)))
            .map(VirtualPosition::from)
        */
    }
}

pub(crate) struct AlnStream<R, B = BamReader<BgzfReader<File>>> {
    ambiguous: Option<BamOutput>,
    pub(crate) bam: Option<B>,
    filt: Option<BamOutput>,
    next: Option<R>,
    output: Option<BamOutput>,
    sample_variants: Option<Arc<Store<Sample>>>,
    population_variants: Option<Arc<Store<Population>>>,
    pub(crate) header: Header,
    output_mode: OutputMode,
    threads: NonZeroUsize,
}

impl<R> AlnStream<R>
where
    R: SimpleRec + FromBamRecord,
{
    pub(crate) fn new(opt: &mut Config, i: usize) -> Result<Self> {
        let bam_str = opt.alignment[i].as_str();
        let file = File::open(bam_str)?;
        let mut bam = BamReader::new(file);
        let header = bam.read_header()?;

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

        let name = test_record
            .name()
            .ok_or_else(|| anyhow!("Record has no name"))?;
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
            .filter(|s| !s.is_empty())
            .map(|s| Store::new(&PathBuf::try_from(s)?, parse_sample_record).map(Arc::new))
            .transpose()?;
        let population_variants = opt
            .population_variants
            .get(i)
            .filter(|s| !s.is_empty())
            .map(|s| Store::new(&PathBuf::try_from(s)?, parse_population_record).map(Arc::new))
            .transpose()?;

        opt.output.get(i).map(path_unicode_ok).transpose()?;
        opt.filtered_output
            .get(i)
            .map(path_unicode_ok)
            .transpose()?;
        opt.ambiguous_output
            .get(i)
            .map(path_unicode_ok)
            .transpose()?;

        let threads = NonZeroUsize::new(opt.threads).unwrap_or(NonZeroUsize::MIN);
        let next_rec = R::from_bam_record(&header, test_record)?;

        Ok(AlnStream {
            ambiguous: None,
            bam: Some(bam),
            filt: None,
            next: Some(next_rec),
            output: None,
            sample_variants,
            population_variants,
            header,
            output_mode: OutputMode::default(),
            threads,
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
            return Err(anyhow!("Cannot un-next more than one record"));
        }
        self.next = Some(rec);
        Ok(())
    }

    fn next_rec(&mut self) -> Result<Option<R>> {
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
        let add_pg = !opt.no_program_line;
        let threads = self.threads.into();
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
        Ok(())
    }

    fn variant_store(&self) -> Option<Arc<dyn StoreTrait>> {
        self.sample_variants
            .as_ref()
            .map(|s| Arc::clone(s) as Arc<dyn StoreTrait>)
            .or_else(|| {
                self.population_variants
                    .as_ref()
                    .map(|p| Arc::clone(p) as Arc<dyn StoreTrait>)
            })
    }

    fn header(&self) -> &Header {
        &self.header
    }

    fn fetch_by_virtual_offset(&mut self, virtual_offset: u64) -> Result<RecordBuf> {
        let vpos = VirtualPosition::try_from(virtual_offset)
            .map_err(|_| anyhow!("Invalid virtual position {virtual_offset}"))?;
        let bam = self
            .bam
            .as_mut()
            .ok_or_else(|| anyhow!("No BAM reader available for seek"))?;
        bam.seek_vpos(vpos)?;
        let rec = bam
            .records()
            .next()
            .ok_or_else(|| anyhow!("No record at virtual offset {virtual_offset}"))?
            .map_err(|e| anyhow!("BAM read error: {e}"))?;
        RecordBuf::try_from_alignment_record(&self.header, &rec)
            .map_err(|e| anyhow!("RecordBuf conversion: {e}"))
    }
}

#[cfg(test)]
pub(crate) mod tests;
