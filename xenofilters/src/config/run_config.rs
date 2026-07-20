use crate::{
    aln_stream::{AlignmentStream, AlnStream},
    config::{
        CommonArgs,
        args::{
            IoArgs, OutputArgsMulti, RegionArgsTabix, ScoringArgs, VariantArgs, ChimericArgs
        },
        MatchingAlgorithm, NameEncoderKind,
    },
    file_spec::path_for_stream,
    filter_algorithm::line_by_line::MAX_STREAMS,
    Error,
    variant::name_to_id::header_name_to_id,
    region::{ScoreFn, ScoredRegions},
    ensure,
};
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{smallvec, SmallVec};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::ops::RangeInclusive;

/// Single flat struct consumed by all three engines. No Args-struct
/// indirection survives past `into_run_config()`. This is what
/// `LineByLine::new`, `HashLookup::new`, `CollatedMatcher::new` accept —
/// identical to before subcommands were introduced.
#[derive(Debug, Default)]
pub(crate) struct RunConfig {
    pub(crate) io: IoArgs,
    pub(crate) scoring: ScoringArgs,
    pub(crate) variants: VariantArgs,
    pub(crate) output: OutputArgsMulti,
    pub(crate) threads: usize,
    pub(crate) chimeric_pairs: Vec<String>,
    pub(crate) stream_labels: Vec<String>,
    pub(crate) region_tabix: RegionArgsTabix,
    pub(crate) name_encoder: Option<NameEncoderKind>,
    pub(crate) is_pass2: bool,
    pub(crate) is_paired: Option<bool>,
}

impl RunConfig {
    pub(crate) fn new(
        common: CommonArgs,
        chimeric: Option<ChimericArgs>,
        name_encoder: Option<NameEncoderKind>,
        threads: usize,
        streams: RangeInclusive<usize>,
        max_stdin: usize,
    ) -> Result<Self, Error> {
        let io = common.io;

        // validate number of streams
        let n = io.validate(max_stdin)?;
        common.scoring.validate()?;
        let max = *streams.end();
        ensure!(streams.contains(&n), Error::InvalidStreamCount { n, min: *streams.start(), max });

        let output = common.output;
        output.validate(if max == 1 { 2 } else { n })?;

        Ok(RunConfig {
            io,
            scoring: common.scoring,
            variants: common.variants,
            output,
            threads,
            score_threads: common.score_threads,
            chimeric_pairs: chimeric.map_or(vec![], |c| c.chimeric_pairs),
            stream_labels: chimeric.map_or(vec![], |c| c.stream_labels),
            region_tabix: common.region_tabix,
            name_encoder,
            ..Default::default()
        })
    }

    /// Open all streams for namesorted/collated (unified reader path).
    /// `bgzf_threads` applies to BAM decompression; ignored for SAM/CRAM.
    pub(crate) fn open_streams_unified(
        &mut self,
        bgzf_threads: usize,
    ) -> Result<SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>, Error> {
        let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];

        // MultithreadedReader parallelises bgzf block decompression.
        // Requires threads > 1 AND a non-seeking backend (namesorted / collated).
        // HashLookup pass-2 uses seek_vpos → must use Single.
        let threads = NonZeroUsize::new(bgzf_threads).unwrap_or(NonZeroUsize::MIN);
        let score_fn = self.scoring.region_score_fn;

        let n = self.io.alignment.len();
        for i in 0..n {
            let path = &self.io.alignment[i];
            let path_str = path.to_string_lossy().to_string();
            tracing::debug!(stream = i, path_str, "Opening stream");
            let positive_regions = &self.variants.positive_regions;
            let name_to_id = header_name_to_id(aln[i].header());
            let positive_regions = path_for_stream(positive_regions, 0).map(|p| ScoredRegions::from_bed(p.as_path(), &name_to_id).map(|s| (s, score_fn))).transpose()?;

            let stream = AlnStream::<RecordBuf>::new(self, i, threads.clone(), positive_regions)?;
            aln.push(Box::new(stream));
            if i > 0 {
                if aln[i].next_qname() != aln[0].next_qname() {
                    return Err(Error::InvalidInput(format!(
                        "HashLookup requires all input streams to be namesorted/collated. \
                         Stream 0 next_qname: {:?}, stream {} next_qname: {:?}",
                        aln[0].next_qname(),
                        i,
                        aln[i].next_qname()
                    )));
                }
            }
        }
        Ok(aln)
    }

    /// Open exactly 2 streams as raw seekable BGZF BAM readers (HashLookup pass-2
    /// requires seek_vpos, which the unified/CRAM-capable reader cannot provide).
    pub(crate) fn open_indexed_streams(
        &mut self,
    ) -> Result<SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>, Error> {
        if !self.io.alignment.iter().any(|p| !p.ends_with(".bam")) {
            return Err(Error::InvalidInput(format!(
                "hashlookup requires BAM/CRAM input. SAM does not support seek: {:?}",
                self.io.alignment
            )));
        }
        self.open_streams()
    }

    pub(crate) fn open_streams(
        &mut self,
    ) -> Result<SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>, Error> {
        let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];
        for i in 0..2 {
            let stream = AlnStream::<RecordBuf>::new_raw_bam(i, self)?;
            aln.push(Box::new(stream));
        }
        Ok(aln)
    }
    pub(crate) fn print_routing_counters(&self, counters: &[u64], tag: &str) {
        use crate::filter_algorithm::line_by_line::core::COUNTER_STRIDE;
        let stream_count = counters.len() / COUNTER_STRIDE;
        for nr in 0..stream_count {
            let b = nr * COUNTER_STRIDE;
            tracing::info!(
                stream = nr,
                backend = tag,
                discard = counters[b],
                out = counters[b + 1],
                ambiguous = counters[b + 2],
                chimeric = counters[b + 3],
                "Stream summary"
            );
        }
        let total: u64 = counters.iter().sum();
        if total == 0 {
            return;
        }
        let ambiguous: u64 = counters.chunks_exact(COUNTER_STRIDE).map(|c| c[2]).sum();
        let frac = ambiguous as f64 / total as f64;
        if frac > self.scoring.warn_ambig_fraction {
            let threshold_phred = match self.scoring.ambiguous_threshold {
                u32::MAX => {
                    if self.io.is_pass2 {
                        0
                    } else {
                        10
                    }
                }
                p => p,
            };
            tracing::warn!(
                ambiguous_pct = format!("{:.1}", frac * 100.0),
                threshold_phred,
                "Ambiguous fraction exceeds warning level."
            );
        }
    }
}
