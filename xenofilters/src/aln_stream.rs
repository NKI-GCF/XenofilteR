use crate::Error;
use crate::alignment::SimpleRec;
use crate::bam::reader::BgzfBamReader;
use crate::bam::{
    BamOutput, SUFFIX_AMBIGUOUS, SUFFIX_FILTERED, expand_header, out_from_file, path_unicode_ok,
    rewrite_rg,
};
use crate::config::MatchingAlgorithm;
use crate::config::run_config::RunConfig;
use crate::file_spec::path_for_stream;
use crate::region::{PositiveRegions, ScoreFn, ScoredRegions, diagnostic::SegregateVariants};
use crate::variant::{
    StoreTrait, build_diagnostic_store_expanded,
    indel_equiv::IndelEquivalenceExpander,
    indel_equiv::corrected::{
        build_population_store_expanded, build_sample_store_expanded, read_vcf_or_bcf_header,
    },
    name_to_id::header_name_to_id,
    population::{Population, parse_population_record},
    sample::{Sample, parse_sample_record},
    store::Store,
};
use noodles::bam::{io::Reader as BamReader, record::Record};
use noodles::bgzf::{VirtualPosition, io::MultithreadedReader};
use noodles::fasta::io::indexed_reader::Builder;
use noodles::sam::Header;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::header::record::value::map::program::tag;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
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
    ) -> Result<Self, Error> {
        let path = opt.io.alignment[i].to_string_lossy();
        let seekable_required = MatchingAlgorithm::Hashlookup == *algorithm;

        let file = File::open(&opt.io.alignment[i])?;

        let mut bam = if !seekable_required && usize::from(threads) > 1 {
            let n_threads = usize::from(threads);
            tracing::debug!(
                stream = i,
                threads = n_threads,
                "Using rayon-backed MultithreadedReader for bgzf decompression"
            );
            let pool = ThreadPoolBuilder::new()
                .num_threads(n_threads)
                .thread_name(move |idx| format!("xenofilters-bgzf-{i}-{idx}"))
                .build()
                .map_err(|e| {
                    Error::InvalidInput(format!("failed to build bgzf thread pool: {e}"))
                })?;
            let pool = Arc::new(pool);
            let reader = pool.install(|| BamReader::from(MultithreadedReader::new(file)));
            BgzfBamReader::Multi { reader, pool }
        } else {
            BgzfBamReader::Single(BamReader::new(file))
        };

        let header = bam.read_header()?;
        let name_to_id = header_name_to_id(&header);
        let positive = &opt.variants.positive_regions;
        let score_fn = opt.variants.region_score_fn;
        let positive_regions = path_for_stream(positive, i)
            .map(|p| ScoredRegions::from_bed(p.as_path(), &name_to_id).map(|s| (s, score_fn)))
            .transpose()?;

        let is_pass2 = header
            .read_groups()
            .keys()
            .any(|id| id.ends_with(SUFFIX_AMBIGUOUS.as_bytes()));

        if is_pass2 {
            opt.is_pass2 = true; // set on RunConfig; overrides threshold selection
            if opt.scoring.ambiguous_threshold == u32::MAX && !header_indicates_bqsr(&header) {
                return Err(Error::Pass2RequiresExplicitThresholdWithoutBqsr { stream: i });
            }
        }

        // Sort-order check (namesorted only).
        if MatchingAlgorithm::Namesorted == *algorithm {
            use noodles::sam::header::record::value::map::header::tag as hd_tag;
            let sorted_or_grouped = header.header().is_some_and(|hd| {
                let of = hd.other_fields();
                let so = of.get(&hd_tag::SORT_ORDER).map(|v| v.to_string());
                let go = of.get(&hd_tag::GROUP_ORDER).map(|v| v.to_string());
                so.as_deref() == Some("coordinate") || go.as_deref() == Some("reference")
            });
            if sorted_or_grouped {
                return Err(Error::CoordinateSortedInputDetected);
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

        opt.is_paired = if i == 0 && opt.is_paired.is_none() {
            Some(test_record.flags().is_segmented())
        } else {
            if opt.is_paired != Some(test_record.flags().is_segmented()) {
                return Err(Error::MixedPairedAndSingleEnd);
            }
            opt.is_paired
        };

        let (sample_variants, population_variants) = build_variant_stores(opt, i, name_to_id)?;

        opt.io.output.get(i).map(path_unicode_ok).transpose()?;
        opt.io
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
            .io
            .output
            .get(i)
            .map(|f| out_from_file(f, &self.header, add_pg, threads))
            .transpose()?;
        let expanded = expand_header(self.header.clone(), opt.io.write_discarded);
        self.ambiguous = opt
            .io
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
/// - `header`    : SAM/BAM header of stream `i` (for name->id mapping)
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
    name_to_id: HashMap<String, usize>, // BAM-header-derived, correct for overlapping_multi's tid
) -> Result<(Option<Arc<dyn StoreTrait>>, Option<Arc<dyn StoreTrait>>), Error> {
    let sample_vcf_path = path_for_stream(&config.variants.sample_variants, stream_idx);
    let population_vcf_path = path_for_stream(&config.variants.population_variants, stream_idx);
    let min_gq = config.variants.min_sample_gq;
    let min_af = config.variants.min_population_af;

    match path_for_stream(&config.io.reference, stream_idx) {
        None => {
            let sample_store = sample_vcf_path
                .map(|p| Store::<Sample>::new_from_path(p, parse_sample_record, min_gq)
                    .map(|s| Arc::new(s) as Arc<dyn StoreTrait>))
                .transpose()?;
            let population_store = population_vcf_path
                .map(|p| Store::<Population>::new_from_path(p, parse_population_record, min_af)
                    .map(|s| Arc::new(s) as Arc<dyn StoreTrait>))
                .transpose()?;
            Ok((sample_store, population_store))
        }
        Some(fasta_path) => {
            let fai_path = fasta_path.with_extension(
                fasta_path.extension().map(|e| format!("{}.fai", e.to_string_lossy())).unwrap_or_else(|| "fai".to_string()),
            );
            if !fai_path.exists() {
                return Err(Error::FastaIndexMissing(fasta_path.to_path_buf()));
            }
            let sample_store = sample_vcf_path
                .map(|p| build_sample_store_expanded(p, fasta_path, &name_to_id, min_gq)
                    .map(|s| Arc::new(s) as Arc<dyn StoreTrait>))
                .transpose()?;
            let population_store = population_vcf_path
                .map(|p| build_population_store_expanded(p, fasta_path, &name_to_id, min_af)
                    .map(|s| Arc::new(s) as Arc<dyn StoreTrait>))
                .transpose()?;
            Ok((sample_store, population_store))
        }
    }
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
) -> Result<Option<SegregateVariants>, Error> {
    let Some(segregate) = &config.segregate else { return Ok(None) };
    let Some(path) = path_for_stream(&segregate.distinct_variants, stream_idx) else { return Ok(None) };
    let name_to_id = header_name_to_id(header);

    match path_for_stream(&config.io.reference, stream_idx) {
        None => SegregateVariants::from_vcf(path, &name_to_id).map(Some),
        Some(fasta_path) => {
            let fasta_reader = Builder::default().build_from_path(fasta_path)?;
            let mut expander = IndelEquivalenceExpander::new(fasta_reader);
            let vcf_header = read_vcf_or_bcf_header(path)?;
            build_diagnostic_store_expanded(path, &mut expander, &name_to_id, &vcf_header).map(Some)
        }
    }
}

/// Heuristically detects whether BQSR (or an equivalent base-quality
/// recalibration step) appears in the BAM's @PG chain, by scanning program
/// IDs, @PG NAME, and @PG CO/command-line fields for known tool signatures.
///
/// This is inherently a heuristic, not a proof: a BAM recalibrated by an
/// unlisted tool will be (safely) treated as "not recalibrated." The
/// failure direction is deliberately conservative -- worst case, a user
/// with a legitimately recalibrated BAM has to pass `--ambiguous-threshold`
/// explicitly once for pass2; the alternative (silently trusting
/// `auto` -> threshold 0 on non-recalibrated data) risks resolving
/// ambiguous reads that a fair per-base error model would not have
/// resolved.
fn header_indicates_bqsr(header: &Header) -> bool {
    const SIGNATURES: &[&str] = &["baserecalibrator", "applybqsr", "bqsr", "recalibrat"];

    header.programs().roots().any(|(id, map)| {
        let id_lower = id.to_string().to_ascii_lowercase();
        if SIGNATURES.iter().any(|s| id_lower.contains(s)) {
            return true;
        }
        let of = map.other_fields();
        [tag::NAME, tag::COMMAND_LINE].iter().any(|t| {
            of.get(t)
                .map(|v| {
                    let v_lower = v.to_string().to_ascii_lowercase();
                    SIGNATURES.iter().any(|s| v_lower.contains(s))
                })
                .unwrap_or(false)
        })
    })
}

#[cfg(any(test, feature = "bench-internals"))]
pub(crate) mod tests;
