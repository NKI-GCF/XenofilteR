use crate::{
    Error,
    aln_stream::{AlignmentStream, AlnStream},
    config::{
        CommonArgs, MatchingAlgorithm, NameEncoderKind,
        args::{ChimericArgs, IoArgs, RelatedArgs, ScoringArgs, SegregateArgs},
    },
    file_spec::path_for_stream,
    region::ScoredRegions,
    variant::name_to_id::header_name_to_id,
};
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{SmallVec, smallvec};
use std::num::NonZeroUsize;
use std::ops::RangeInclusive;

/// Single flat struct consumed by all three engines. No Args-struct
/// indirection survives past `into_run_config()`. This is what
/// `LineByLine::new`, `HashLookup::new`, `CollatedMatcher::new` accept --
/// identical to before subcommands were introduced.
#[derive(Debug, Default)]
pub(crate) struct RunConfig {
    pub(crate) io: IoArgs,
    pub(crate) scoring: ScoringArgs,
    pub(crate) variants: RelatedArgs,
    pub(crate) threads: usize,
    pub(crate) chimeric_pairs: Vec<String>,
    pub(crate) stream_labels: Vec<String>,
    pub(crate) segregate: Option<SegregateArgs>,
    pub(crate) name_encoder: Option<NameEncoderKind>,
    pub(crate) is_pass2: bool,
    pub(crate) is_paired: Option<bool>,
}

impl RunConfig {
    pub(crate) fn new(
        mut common: CommonArgs,
        chimeric: Option<ChimericArgs>,
        name_encoder: Option<NameEncoderKind>,
        threads: usize,
        streams: RangeInclusive<usize>,
        max_stdin: usize,
        segregate: Option<SegregateArgs>,
    ) -> Result<Self, Error> {
        let io = common.io;

        // validate number of streams
        io.validate(max_stdin, streams)?;
        common.scoring.validate()?;

        let (chimeric_pairs, stream_labels) = match chimeric {
            Some(c) => (c.chimeric_pairs, c.stream_labels),
            None => (vec![], vec![]),
        };

        Ok(RunConfig {
            io,
            scoring: common.scoring,
            variants: common.variants,
            threads,
            chimeric_pairs,
            stream_labels,
            segregate,
            name_encoder,
            ..Default::default()
        })
    }

    /// Open all streams for namesorted/collated (unified reader path).
    /// `bgzf_threads` applies to BAM decompression; ignored for SAM.
    pub(crate) fn open_streams_unified(
        &mut self,
        algorithm: MatchingAlgorithm,
        bgzf_threads: usize,
    ) -> Result<SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>, Error> {
        let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];

        // MultithreadedReader parallelises bgzf block decompression.
        // Requires threads > 1 AND a non-seeking backend (namesorted / collated).
        // HashLookup pass-2 uses seek_vpos -> must use Single.
        let threads = NonZeroUsize::new(bgzf_threads).unwrap_or(NonZeroUsize::MIN);
        let score_fn = self.variants.region_score_fn;

        let n = self.io.alignment.len();
        for i in 0..n {
            let path = &self.io.alignment[i];
            let path_str = path.to_string_lossy().to_string();
            tracing::debug!(stream = i, path_str, "Opening stream");
            let positive_regions = &self.variants.positive_regions;
            let name_to_id = header_name_to_id(aln[i].header());
            let positive_regions = path_for_stream(positive_regions, 0)
                .map(|p| ScoredRegions::from_bed(p.as_path(), &name_to_id).map(|s| (s, score_fn)))
                .transpose()?;

            let stream =
                AlnStream::<RecordBuf>::new(self, &algorithm, i, threads, positive_regions)?;
            aln.push(Box::new(stream));
            if i > 0 && aln[i].next_qname() != aln[0].next_qname() {
                return Err(Error::InvalidInput(format!(
                    "HashLookup requires all input streams to be namesorted/collated. \
                         Stream 0 next_qname: {:?}, stream {} next_qname: {:?}",
                    aln[0].next_qname(),
                    i,
                    aln[i].next_qname()
                )));
            }
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
