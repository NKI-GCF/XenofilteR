use crate::alignment::SimpleRec;
use crate::bam::reader::BgzfBamReader;
use crate::bam::{
    expand_header, out_from_file, path_unicode_ok, rewrite_rg, BamOutput, SUFFIX_AMBIGUOUS,
    SUFFIX_FILTERED,
};
use crate::config::run_config::RunConfig;
use crate::config::{MatchingAlgorithm, StripReadSuffix};
use crate::file_spec::path_for_stream;
use crate::region::{diagnostic::SegregateVariants, PositiveRegions, ScoreFn};
use crate::variant::{
    build_diagnostic_store_expanded,
    indel_equiv::corrected::read_vcf_or_bcf_header,
    indel_equiv::{
        build_population_store_expanded, build_sample_store_expanded, IndelEquivalenceExpander,
    },
    name_to_id::header_name_to_id,
    population::{parse_population_record, Population},
    sample::{parse_sample_record, Sample},
    store::Store,
    StoreTrait,
};
use crate::Error;
use noodles::bam::{io::Reader as BamReader, record::Record};
use noodles::bgzf::{io::MultithreadedReader, VirtualPosition};
use noodles::fasta::io::indexed_reader::Builder;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::Header;
use std::fs::File;
use std::num::NonZeroUsize;
use std::sync::Arc;

pub(crate) trait AlignmentStream<R: SimpleRec> {
    fn next_qname(&self) -> &[u8];
    fn un_next(&mut self, rec: R) -> Result<(), Error>;
    fn next_rec(&mut self) -> Result<Option<R>, Error>;
    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<(), Error>;
    fn init_writers(&mut self, _opt: &RunConfig, _i: usize) -> Result<(), Error> {
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

pub(crate) struct AlnStream<R> {
    pub(crate) bam: Option<BgzfBamReader>,
    next: Option<R>,
    sample_variants: Option<Arc<dyn StoreTrait>>,
    population_variants: Option<Arc<dyn StoreTrait>>,
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
        opt: &mut RunConfig,
        algorithm: &MatchingAlgorithm,
        i: usize,
        threads: NonZeroUsize,
        positive_regions: Option<(PositiveRegions, ScoreFn)>,
    ) -> Result<Self, Error> {
        let path = opt.io.alignment[i].to_string_lossy();
        let seekable_required = MatchingAlgorithm::Hashlookup == *algorithm;

        let file = File::open(&opt.io.alignment[i])?;

        let mut bam = if !seekable_required && usize::from(threads) > 1 {
            tracing::debug!(
                stream = i,
                threads = usize::from(threads),
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
            opt.is_pass2 = true; // set on RunConfig; overrides threshold selection
        }

        // Sort-order check (namesorted only).
        let raw = bam.read_raw_header_bytes()?;
        if MatchingAlgorithm::Namesorted == *algorithm {
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
        opt.io.strip_read_suffix = match opt.io.strip_read_suffix {
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

        let (sample_variants, population_variants) = build_variant_stores(opt, i)?;

        opt.output.output.get(i).map(path_unicode_ok).transpose()?;
        opt.output
            .ambiguous_output
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
            write_discarded: opt.io.write_discarded,
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

    fn init_writers(&mut self, opt: &RunConfig, i: usize) -> Result<(), Error> {
        let add_pg = !opt.io.no_program_line;
        let threads = self.threads.into();
        self.output = opt
            .output
            .output
            .get(i)
            .map(|f| out_from_file(f, &self.header, add_pg, threads))
            .transpose()?;
        let expanded = expand_header(self.header.clone(), opt.io.write_discarded);
        self.ambiguous = opt
            .output
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
fn build_variant_stores(
    config: &RunConfig,
    stream_idx: usize,
) -> Result<(Option<Arc<dyn StoreTrait>>, Option<Arc<dyn StoreTrait>>), Error> {
    // Resolve the VCF paths for this stream index.
    let sample_vcf_path = path_for_stream(&config.variants.sample_variants, stream_idx);
    let population_vcf_path = path_for_stream(&config.variants.population_variants, stream_idx);

    if !config.variants.expand_indels {
        // -- Plain path: existing behavior, unchanged -------------------------
        let sample_store: Option<Arc<dyn StoreTrait>> = sample_vcf_path
            .map(|p| {
                Store::<Sample>::new_from_path(p, parse_sample_record)
                    .map(|s| Arc::new(s) as Arc<dyn StoreTrait>)
            })
            .transpose()?;

        let population_store: Option<Arc<dyn StoreTrait>> = population_vcf_path
            .map(|p| {
                Store::<Population>::new_from_path(p, parse_population_record)
                    .map(|s| Arc::new(s) as Arc<dyn StoreTrait>)
            })
            .transpose()?;

        return Ok((sample_store, population_store));
    }

    // -- Expanded path: requires --reference ----------------------------------
    let fasta_path = path_for_stream(&config.io.reference, stream_idx);
    let fasta_path = fasta_path
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
        return Err(crate::Error::FastaIndexMissing(fasta_path.to_path_buf()));
    }

    let sample_store: Option<Arc<dyn StoreTrait>> = sample_vcf_path
        .map(|p| {
            build_sample_store_expanded(p, fasta_path).map(|s| Arc::new(s) as Arc<dyn StoreTrait>)
        })
        .transpose()?;

    let population_store: Option<Arc<dyn StoreTrait>> = population_vcf_path
        .map(|p| {
            build_population_store_expanded(p, fasta_path)
                .map(|s| Arc::new(s) as Arc<dyn StoreTrait>)
        })
        .transpose()?;

    Ok((sample_store, population_store))
}

/// Build the diagnostic variant store for one stream.
///
/// Separated from `build_variant_stores` because the diagnostic store has
/// a different type (`SegregateVariants`) and is consumed by different
/// code paths (HashLookup BED/VCF region check).
pub(crate) fn build_diagnostic_store_for_stream(
    config: &RunConfig,
    header: &Header,
    stream_idx: usize,
) -> Result<Option<crate::region::diagnostic::SegregateVariants>, Error> {
    if let Some (segregate) = &config.segregate {
        let diag_vcf_path = path_for_stream(&segregate.distinct_variants, stream_idx);
        let Some(path) = diag_vcf_path else {
            return Ok(None);
        };
        if !config.variants.expand_indels {
            // Plain diagnostic loading: existing behavior.
            return SegregateVariants::from_vcf(path, &header_name_to_id(header)).map(Some);
        }
        let fasta_path = path_for_stream(&config.io.reference, stream_idx);
        let fasta_path = fasta_path
            .as_deref()
            .ok_or(crate::Error::ExpandIndelsRequiresReference)?;

        let fasta_reader = Builder::default().build_from_path(fasta_path)?;

        let mut expander = IndelEquivalenceExpander::new(fasta_reader);
        let name_to_id = header_name_to_id(header);

        // Read the VCF header separately so we can pass it by reference.
        let vcf_header = read_vcf_or_bcf_header(path)?;

        build_diagnostic_store_expanded(path, &mut expander, &name_to_id, &vcf_header).map(Some)
    } else {
        return Ok(None);
    }


}

#[cfg(test)]
pub(crate) mod tests;
