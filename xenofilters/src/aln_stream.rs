mod cram;
mod sam_stdin;

use crate::alignment::SimpleRec;
use crate::bam::{out_from_file, path_unicode_ok, OutputMode};
use crate::config::MatchingAlgorithm;
use crate::config::{Config, StripReadSuffix};
use crate::variant::{
    parse_population_record, parse_sample_record, Population, Sample, Store, StoreTrait,
};
use crate::Error;
use noodles::bam::{io::Reader as BamReader, record::Record};
use noodles::bgzf::{self, io::MultithreadedReader, VirtualPosition};
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::Header;
use std::fs::File;
use std::io::Read as ioRead;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Arc;

pub(crate) use cram::CramStream;
pub(crate) use sam_stdin::SamStdinStream;

pub(crate) trait AlignmentStream<R: SimpleRec> {
    fn next_qname(&self) -> &[u8];
    fn un_next(&mut self, rec: R) -> Result<(), Error>;
    fn next_rec(&mut self) -> Result<Option<R>, Error>;
    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<(), Error>;
    fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<(), Error> {
        Ok(())
    }
    fn variant_store(&self) -> Option<Arc<dyn StoreTrait>> {
        None
    }
    fn header(&self) -> &Header;
    /// Seek to a BGZF virtual offset and read one full record for pass-2 output.
    /// Returns `Err` for stream types that do not support seeking (e.g. mock streams).
    fn fetch_by_virtual_offset(&mut self, _virtual_offset: u64) -> Result<RecordBuf, Error> {
        Err(Error::FetchByVirtualOffsetNotSupported)
    }
}

/// Wraps both bgzf reader variants in a zero-cost internal enum.
/// `Single` is seekable (required by HashLookup pass-2).
/// `Multi` is not seekable but parallelises decompression.
pub(crate) enum BgzfBamReader {
    Single(BamReader<bgzf::io::Reader<File>>),
    Multi(BamReader<MultithreadedReader<File>>),
}

impl BgzfBamReader {
    fn read_header(&mut self) -> std::io::Result<Header> {
        match self {
            Self::Single(r) => r.read_header(),
            Self::Multi(r) => r.read_header(),
        }
    }

    /// Returns an error for the multithreaded variant.
    /// HashLookup must always be constructed with `Single`.
    fn seek_vpos(&mut self, pos: VirtualPosition) -> std::io::Result<VirtualPosition> {
        match self {
            Self::Single(r) => r.get_mut().seek(pos),
            Self::Multi(_) => Err(std::io::Error::new(
                std::io::ErrorKind::Unsupported,
                "seek_vpos unavailable on MultithreadedReader; \
                 use --matching-algorithm namesorted or collated with --threads",
            )),
        }
    }

    fn next_record(&mut self) -> Option<std::io::Result<Record>> {
        match self {
            Self::Single(r) => r.records().next(),
            Self::Multi(r) => r.records().next(),
        }
    }

    /// Raw header bytes for SO: / GO: sort-order validation.
    /// Both variants can read sequentially from position 0.
    fn read_raw_header_bytes(&mut self) -> std::io::Result<Vec<u8>> {
        // Re-read the raw SAM header via the inner bgzf reader.
        // After `read_header()` the file cursor is at the first record;
        // this secondary read is only used during `new()` before any
        // records are consumed.
        match self {
            Self::Single(r) => {
                let mut hr = r.header_reader();
                hr.read_magic_number()?;
                let mut rhr = hr.raw_sam_header_reader()?;
                let mut buf = Vec::new();
                rhr.read_to_end(&mut buf)?;
                Ok(buf)
            }
            Self::Multi(r) => {
                let mut hr = r.header_reader();
                hr.read_magic_number()?;
                let mut rhr = hr.raw_sam_header_reader()?;
                let mut buf = Vec::new();
                rhr.read_to_end(&mut buf)?;
                Ok(buf)
            }
        }
    }
}

pub(crate) struct AlnStream<R> {
    pub(crate) bam: Option<BgzfBamReader>,
    next: Option<R>,
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
    pub(crate) fn new(opt: &mut Config, i: usize) -> Result<Self, Error> {
        let path = opt.alignment[i].as_str();
        let file = File::open(path)?;

        // MultithreadedReader parallelises bgzf block decompression.
        // Requires threads > 1 AND a non-seeking backend (namesorted / collated).
        // HashLookup pass-2 uses seek_vpos → must use Single.
        let seekable_required = opt.matching_algorithm == MatchingAlgorithm::Hashlookup;
        let threads = NonZeroUsize::new(opt.threads).unwrap_or(NonZeroUsize::MIN);

        let mut bam = if !seekable_required && opt.threads > 1 {
            tracing::debug!(
                stream = i,
                threads = opt.threads,
                "Using MultithreadedReader for bgzf decompression"
            );
            BgzfBamReader::Multi(BamReader::from(MultithreadedReader::with_worker_count(
                threads, file,
            )))
        } else {
            BgzfBamReader::Single(BamReader::new(file))
        };

        let header = bam.read_header()?;

        // Sort-order check (namesorted only).
        let raw = bam.read_raw_header_bytes()?;
        if opt.matching_algorithm == MatchingAlgorithm::Namesorted {
            for parts in raw
                .split(|&b| b == b'\n')
                .map(|s| s.split(|&b| b == b'\t').collect::<Vec<_>>())
            {
                if parts.len() >= 3
                    && parts[0] == b"@HD"
                    && (parts[2] == b"SO:coordinate" || parts[2] == b"GO:reference")
                {
                    return Err(Error::CoordinateSortedInputDetected);
                }
            }
        }

        let test_record = match bam.next_record() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(e.into()),
            None => {
                return Err(Error::BamHasNoRecords {
                    bam_str: path.to_string(),
                });
            }
        };

        let name = test_record.name().ok_or(Error::RecordHasNoReadName)?;
        opt.strip_read_suffix = match opt.strip_read_suffix {
            StripReadSuffix::True => {
                if !name.ends_with(b"/1") && !name.ends_with(b"/2") {
                    return Err(Error::ReadNamesMissingSuffixes);
                }
                StripReadSuffix::True
            }
            StripReadSuffix::False => {
                if name.ends_with(b"/1") || name.ends_with(b"/2") {
                    return Err(Error::ReadNamesHaveSuffixes);
                }
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
            if opt.is_paired != Some(test_record.flags().is_segmented()) {
                return Err(Error::MixedPairedAndSingleEnd);
            }
            opt.is_paired
        };

        let sample_variants = opt
            .sample_variants
            .get(i)
            .filter(|s| !s.is_empty())
            .map(|s| Store::new(&PathBuf::from(s), parse_sample_record).map(Arc::new))
            .transpose()?;
        let population_variants = opt
            .population_variants
            .get(i)
            .filter(|s| !s.is_empty())
            .map(|s| Store::new(&PathBuf::from(s), parse_population_record).map(Arc::new))
            .transpose()?;

        opt.output.get(i).map(path_unicode_ok).transpose()?;
        opt.discarded_output
            .get(i)
            .map(path_unicode_ok)
            .transpose()?;
        opt.ambiguous_output
            .get(i)
            .map(path_unicode_ok)
            .transpose()?;

        let next_rec = R::from_bam_record(&header, test_record)?;

        Ok(AlnStream {
            bam: Some(bam),
            next: Some(next_rec),
            sample_variants,
            population_variants,
            header,
            output_mode: OutputMode::MultiFile {
                output: None,
                filt: None,
                ambiguous: None,
            },
            threads,
        })
    }
}

pub(crate) trait FromBamRecord: Sized {
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
    fn next_rec(&mut self) -> Result<Option<R>, Error> {
        let header = &self.header;
        Ok(self
            .next
            .take()
            .map(Ok)
            .or_else(|| {
                self.bam.as_mut().and_then(|b| {
                    b.next_record()
                        .map(|r| r.and_then(|r| R::from_bam_record(header, r)))
                })
            })
            .transpose()?)
    }

    fn fetch_by_virtual_offset(&mut self, virtual_offset: u64) -> Result<RecordBuf, Error> {
        let vpos = VirtualPosition::from(virtual_offset);
        let bam = self.bam.as_mut().ok_or(Error::NoBamReaderForSeek)?;
        bam.seek_vpos(vpos).map_err(|e| {
            // Provide a clear message when misused with the wrong backend.
            if e.kind() == std::io::ErrorKind::Unsupported {
                Error::FetchByVirtualOffsetNotSupported
            } else {
                Error::from(e)
            }
        })?;
        let rec = bam
            .next_record()
            .ok_or(Error::NoRecordAtVirtualOffset { virtual_offset })??;
        RecordBuf::try_from_alignment_record(&self.header, &rec)
            .map_err(|e| Error::RecordBufConversion(e.to_string()))
    }

    fn next_qname(&self) -> &[u8] {
        self.next
            .as_ref()
            .and_then(|r| r.name())
            .map_or(b"", |n| n.as_ref())
    }

    fn un_next(&mut self, rec: R) -> Result<(), Error> {
        if self.next.is_some() {
            return Err(Error::CannotUnNext);
        }
        self.next = Some(rec);
        Ok(())
    }

    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<(), Error> {
        self.output_mode.write(rec, is_best, &self.header)
    }

    fn init_writers(&mut self, opt: &Config, i: usize) -> Result<(), Error> {
        let add_pg = !opt.no_program_line;
        let threads = self.threads.into();

        if let Some(ref path) = opt.merged_output {
            if i == 0 {
                // Stream 0 owns the single merged writer.
                // KNOWN LIMITATION: MergedOutput is owned per-stream (stream 0 only).
                // Correct multi-stream merged output requires lifting MergedOutput into
                // LineByLine and routing records from all streams through a shared writer
                // (Arc<Mutex<MergedOutput>> or caller-level aggregation).
                // Tracked in TODO_rust.md.
                use crate::bam::{expand_header, MergedOutput};
                let expanded = expand_header(self.header.clone());
                self.output_mode =
                    OutputMode::Merged(MergedOutput::new(path, expanded, add_pg, threads)?);
            }
            // i > 0: no-op writer; all records written through stream 0's MergedOutput.
            // This is a known limitation; correct multi-stream merged output requires
            // an Arc<Mutex<MergedOutput>> or caller-level aggregation.
        } else {
            self.output_mode = OutputMode::MultiFile {
                output: opt
                    .output
                    .get(i)
                    .map(|f| out_from_file(f, &self.header, add_pg, threads))
                    .transpose()?,
                filt: opt
                    .discarded_output
                    .get(i)
                    .map(|f| out_from_file(f, &self.header, add_pg, threads))
                    .transpose()?,
                ambiguous: opt
                    .ambiguous_output
                    .get(i)
                    .map(|f| out_from_file(f, &self.header, add_pg, threads))
                    .transpose()?,
            };
        }
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
}

#[cfg(test)]
pub(crate) mod tests;
