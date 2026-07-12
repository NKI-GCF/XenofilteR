use crate::alignment::SimpleRec;
use crate::bam::{
    expand_header, out_from_file, path_unicode_ok, rewrite_rg, BamOutput, SUFFIX_AMBIGUOUS,
    SUFFIX_FILTERED,
};
use crate::config::MatchingAlgorithm;
use crate::config::{Config, StripReadSuffix};
use crate::region::{PositiveRegions, ScoreFn};
use crate::variant::{
    indel_equiv::{
        build_diagnostic_store_expanded, build_population_store_expanded,
        build_sample_store_expanded,
    },
    name_to_id::header_name_to_id,
    population::Population,
    sample::Sample,
    store::Store,
    StoreTrait,
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
    fn positive_regions(&self) -> Option<(&PositiveRegions, ScoreFn)> {
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
    output: Option<BamOutput>,
    ambiguous: Option<BamOutput>,
    write_discarded: bool,
    threads: NonZeroUsize,
    positive_regions: Option<(PositiveRegions, ScoreFn)>,
}

impl<R> AlnStream<R>
where
    R: SimpleRec + FromBamRecord,
{
    pub(crate) fn new(
        opt: &mut Config,
        i: usize,
        positive_regions: Option<(PositiveRegions, ScoreFn)>,
    ) -> Result<Self, Error> {
        let path = opt.alignment[i].as_str();
        let seekable_required = opt.matching_algorithm == MatchingAlgorithm::Hashlookup;

        let file = File::open(path)?;

        // MultithreadedReader parallelises bgzf block decompression.
        // Requires threads > 1 AND a non-seeking backend (namesorted / collated).
        // HashLookup pass-2 uses seek_vpos → must use Single.
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

        let is_pass2 = header
            .read_groups()
            .keys()
            .any(|id| id.ends_with(SUFFIX_AMBIGUOUS.as_bytes()));

        if is_pass2 {
            opt.is_pass2 = true; // set on Config; overrides threshold selection
        }

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

        let (sample_variants, population_variants) = build_variant_stores(opt, &header, i)?;

        opt.output.get(i).map(path_unicode_ok).transpose()?;
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
            output: None,
            ambiguous: None,
            positive_regions,
            write_discarded: opt.write_discarded,
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

    fn write_record(&mut self, mut rec: RecordBuf, is_best: Option<bool>) -> Result<(), Error> {
        let (dest, suffix) = match is_best {
            Some(true) => (self.output.as_mut(), None),
            Some(false) => match self.write_discarded {
                true => (self.ambiguous.as_mut(), Some(SUFFIX_FILTERED)),
                false => (None, None),
            },
            None => (self.ambiguous.as_mut(), Some(SUFFIX_AMBIGUOUS)),
        };
        if let Some(sfx) = suffix {
            rewrite_rg(&mut rec, sfx)?;
        }
        if let Some(w) = dest {
            w.write_alignment_record(&self.header, &rec)?;
        }
        Ok(())
    }

    fn init_writers(&mut self, opt: &Config, i: usize) -> Result<(), Error> {
        let add_pg = !opt.no_program_line;
        let threads = self.threads.into();
        self.output = opt
            .output
            .get(i)
            .map(|f| out_from_file(f, &self.header, add_pg, threads))
            .transpose()?;
        let expanded = expand_header(self.header.clone(), opt);
        self.ambiguous = opt
            .ambiguous_output
            .get(i)
            .map(|f| out_from_file(f, &expanded, add_pg, threads))
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

    fn positive_regions(&self) -> Option<(&PositiveRegions, ScoreFn)> {
        self.positive_regions.as_ref().map(|(pr, sf)| (pr, *sf))
    }

    fn header(&self) -> &Header {
        &self.header
    }
}

/// Build all three variant stores for alignment stream `i`, choosing between
/// plain loading and indel-equivalence-expanded loading based on
/// `config.expand_indels`.
///
/// Called from `AlnStream::new()` after the BAM header has been read.
///
/// # Arguments
/// - `config`    : parsed CLI args (CommonArgs + subcommand-specific fields)
/// - `header`    : SAM/BAM header of stream `i` (for name→id mapping)
/// - `stream_idx`: which stream we are building stores for (0-based)
///
/// # Returns
/// `(sample_store, population_store, diagnostic_vcf_path_if_any)`
///
/// The diagnostic store for HashLookup is built as an in-memory
/// `AmbiguousRegions`-style struct; see `build_diagnostic_store_expanded`.
pub(crate) fn build_variant_stores(
    config: &CommonArgs,
    header: &Header,
    stream_idx: usize,
) -> Result<(Option<Arc<dyn StoreTrait>>, Option<Arc<dyn StoreTrait>>)> {
    // Resolve the VCF paths for this stream index.
    let sample_vcf_path = variant_path_for_stream(&config.sample_variants, stream_idx);
    let population_vcf_path = variant_path_for_stream(&config.population_variants, stream_idx);

    if !config.expand_indels {
        // -- Plain path: existing behavior, unchanged -------------------------
        let sample_store: Option<Arc<dyn StoreTrait>> = sample_vcf_path
            .map(|p| {
                Store::<Sample>::from_vcf(p, header).map(|s| Arc::new(s) as Arc<dyn StoreTrait>)
            })
            .transpose()?;

        let population_store: Option<Arc<dyn StoreTrait>> = population_vcf_path
            .map(|p| {
                Store::<Population>::from_vcf(p, header).map(|s| Arc::new(s) as Arc<dyn StoreTrait>)
            })
            .transpose()?;

        return Ok((sample_store, population_store));
    }

    // -- Expanded path: requires --reference ----------------------------------
    let fasta_path = config
        .reference
        .as_deref()
        .ok_or(crate::Error::ExpandIndelsRequiresReference)?;

    // Validate that the .fai sidecar exists before opening the expander.
    let fai_path = fasta_path.with_extension(
        fasta_path
            .extension()
            .map(|e| format!("{}.fai", e.to_string_lossy()))
            .unwrap_or_else(|| "fai".to_string()),
    );
    if !fai_path.exists() {
        return Err(crate::Error::FastaIndexMissing(fasta_path.to_path_buf()).into());
    }

    let name_to_id = header_name_to_id(header);

    let sample_store: Option<Arc<dyn StoreTrait>> = sample_vcf_path
        .map(|p| {
            build_sample_store_expanded(
                p,
                fasta_path,
                &name_to_id,
                config.gamete_streams.contains(&stream_idx),
            )
            .map(|s| Arc::new(s) as Arc<dyn StoreTrait>)
        })
        .transpose()?;

    let population_store: Option<Arc<dyn StoreTrait>> = population_vcf_path
        .map(|p| {
            build_population_store_expanded(p, fasta_path, &name_to_id)
                .map(|s| Arc::new(s) as Arc<dyn StoreTrait>)
        })
        .transpose()?;

    Ok((sample_store, population_store))
}

/// Build the diagnostic variant store for one stream.
///
/// Separated from `build_variant_stores` because the diagnostic store has
/// a different type (`DiagnosticVariants`) and is consumed by different
/// code paths (HashLookup BED/VCF region check).
pub(crate) fn build_diagnostic_store_for_stream(
    config: &CommonArgs,
    header: &Header,
    stream_idx: usize,
) -> Result<Option<crate::region::diagnostic::DiagnosticVariants>> {
    use crate::region::diagnostic::DiagnosticVariants;
    use crate::variant::name_to_id::header_name_to_id;

    let diag_vcf_path = variant_path_for_stream(&config.diagnostic_variants, stream_idx);
    let Some(path) = diag_vcf_path else {
        return Ok(None);
    };

    if !config.expand_indels {
        // Plain diagnostic loading: existing behavior.
        return DiagnosticVariants::from_vcf(path, &header_name_to_id(header)).map(Some);
    }

    let fasta_path = config
        .reference
        .as_deref()
        .ok_or(crate::Error::ExpandIndelsRequiresReference)?;

    use crate::variant::indel_equiv::IndelEquivalenceExpander;
    use noodles::fasta;
    use std::{fs::File, io::BufReader};

    let fasta_reader = fasta::IndexedReader::new(BufReader::new(File::open(fasta_path)?));
    let mut expander = IndelEquivalenceExpander::new(fasta_reader);
    let name_to_id = header_name_to_id(header);

    // Read the VCF header separately so we can pass it by reference.
    let vcf_header = crate::variant::indel_equiv::read_vcf_or_bcf_header(path)?;

    build_diagnostic_store_expanded(path, &mut expander, &name_to_id, &vcf_header).map(Some)
}

// -- Private helper ------------------------------------------------------------

/// Resolve the VCF/BCF path for stream `stream_idx` from a list of
/// `[IDX:]FILE` strings.
///
/// Accepts two forms:
/// - `"FILE"` — positional: stream 0 gets index 0, stream 1 gets index 1, …
/// - `"IDX:FILE"` — explicit stream index
fn variant_path_for_stream(specs: &[String], stream_idx: usize) -> Option<&std::path::Path> {
    for spec in specs {
        if let Some((idx_str, path_str)) = spec.split_once(':') {
            if let Ok(idx) = idx_str.trim().parse::<usize>() {
                if idx == stream_idx {
                    return Some(std::path::Path::new(path_str.trim()));
                }
            }
        }
    }
    // Positional fallback: no IDX: prefix → assign by position in the list.
    let positional: Vec<&str> = specs
        .iter()
        .filter(|s| !s.contains(':') || s.starts_with('/') || s.starts_with('.'))
        .map(|s| s.as_str())
        .collect();
    positional.get(stream_idx).map(|s| std::path::Path::new(s))
}

#[cfg(test)]
pub(crate) mod tests;
