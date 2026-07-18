use crate::{
    aln_stream::{AlignmentStream, AlnStream},
    config::{
        args::{
            IoArgs, OutputArgsMulti, RegionArgsMemory, RegionArgsTabix, ScoringArgs, VariantArgs,
        },
        run_config::RunConfig,
        MatchingAlgorithm, NameEncoderKind,
    },
    file_spec::path_for_stream,
    filter_algorithm::line_by_line::MAX_STREAMS,
    Error,
    ensure,
};
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{smallvec, SmallVec};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::ops::Range;

/// Single flat struct consumed by all three engines. No Args-struct
/// indirection survives past `into_run_config()`. This is what
/// `LineByLine::new`, `HashLookup::new`, `CollatedMatcher::new` accept —
/// identical to before subcommands were introduced.
#[derive(Debug, Default)]
pub(crate) struct RunConfig {
    pub(crate) single_alignment_mode: bool,
    pub(crate) io: IoArgs,
    pub(crate) scoring: ScoringArgs,
    pub(crate) variants: VariantArgs,
    pub(crate) output: OutputArgsMulti,
    pub(crate) threads: usize,
    pub(crate) score_threads: usize,
    pub(crate) chimeric_pairs: Vec<String>,
    pub(crate) stream_labels: Vec<String>,
    pub(crate) region_memory: RegionArgsMemory,
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
        rng: Range<usize>,
    } -> Result<Self, Error> {
        let n = common.io.alignment.len();
        let max = rng.end;
        ensure!(rng.contains(n), Error::InvalidStreamCount { n, min: rng.start, max });

        let output = common.output.clone();
        output.validate(if max == 1 { 2 } else { n })?;

        // CRAM sanity: any .cram input requires --reference.
        let cram_lacks_ref = common.io.alignment
            .iter()
            .enumerate()
            .any(|(i, p)| p.ends_with(".cram") && path_for_stream(&common.io.reference, i).is_none());
        ensure!(!cram_lacks_ref, Error::CramRequiresReference);

        let mut config = RunConfig {
            algorithm: common.algorithm,
            single_alignment_mode: common.single_alignment_mode,
            io: common.io,
            scoring: common.scoring,
            variants: common.variants,
            output,
            threads,
            score_threads: common.score_threads,
            chimeric_pairs: chimeric.map_or(vec![], |c| c.pairs),
            stream_labels: chimeric.map_or(vec![], |c| c.stream_labels),
            region_memory: common.region_memory,
            region_tabix: common.region_tabix,
            name_encoder,
            ...Default::default()
        };
        Ok(config)
    }

    pub(crate) fn validate_and_init(&mut self) -> Result<(), Error> {
        self.io.validate()?;
        self.scoring.validate()?;

        let n = self.alignment.len();

        match self.algorithm {
            MatchingAlgorithm::Namesorted => {
                if n < 1 || n > MAX_STREAMS {
                    return Err(Error::NamesortedStreamCount { got: n });
                }
                self.output.validate(n)?;
            }
            MatchingAlgorithm::Hashlookup => {
                if n != 2 {
                    return Err(Error::HashlookupStreamCount { got: n });
                }
                if self.score_threads > 1 {
                    return Err(Error::AlgorithmNotParallel {
                        algorithm: "hashlookup",
                    });
                }
                self.output.validate(2)?;
            }
            MatchingAlgorithm::Collated => {
                if n != 2 {
                    return Err(Error::CollatedStreamCount { got: n });
                }
                if self.score_threads > 1 {
                    return Err(Error::AlgorithmNotParallel {
                        algorithm: "collated",
                    });
                }
                self.output.validate(2)?;
            }
            MatchingAlgorithm::Strain => {
                if n != 1 {
                    return Err(Error::StrainStreamCount { got: n });
                }
                self.output.validate(2)?;
            }
            MatchingAlgorithm::ViralIntegration => {
                if n < 2 {
                    return Err(Error::ViralIntegrationStreamCount { got: n });
                }
                self.output.validate(n)?;
            }
        }

        // CRAM sanity: any .cram input requires --reference.
        if self
            .alignment
            .iter()
            .enumerate()
            .any(|(i, p)| p.ends_with(".cram") && path_for_stream(&self.io.reference, i).is_none())
        {
            return Err(Error::CramRequiresReference);
        }

        // stdin: at most one stream, namesorted only.
        let stdin_count = self
            .alignment
            .iter()
            .filter(|p| p.to_string_lossy() == "-")
            .count();
        if stdin_count > 1 {
            return Err(Error::MultipleStdinStreams);
        }
        if stdin_count == 1 && self.algorithm != MatchingAlgorithm::Namesorted {
            return Err(Error::StdinNamesortedOnly);
        }

        Ok(())
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

        let n = self.alignment.len();
        for i in 0..n {
            let path = &self.alignment[i];
            let path_str = path.to_string_lossy().to_string();
            tracing::debug!(stream = i, path_str, "Opening stream");
            let positive_regions = path_for_stream(specs, 0)
                .map(|s| 


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
    pub(crate) fn open_streams_raw_bam(
        &mut self,
    ) -> Result<SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>, Error> {
        if !self.alignment.iter().all(|p| p.ends_with(".bam")) {
            return Err(Error::InvalidInput(format!(
                "hashlookup requires BAM input (CRAM/SAM seeking not yet supported; \
                 see ROADMAP: CRAM index seeking): {:?}",
                self.alignment
            )));
        }
        let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];
        for i in 0..2 {
            let stream = AlnStream::<RecordBuf>::new_raw_bam(i, self)?;
            aln.push(Box::new(stream));
        }
        Ok(aln)
    }
}
