use crate::config::args::{
    IoArgs, OutputArgsMulti, RegionArgsMemory,
    RegionArgsTabix, ScoringArgs, VariantArgs,
};
use crate::config::MatchingAlgorithm;
use crate::{Error, filter_algorithm::line_by_line::MAX_STREAMS};

/// Single flat struct consumed by all three engines. No Args-struct
/// indirection survives past `into_run_config()`. This is what
/// `LineByLine::new`, `HashLookup::new`, `CollatedMatcher::new` accept —
/// identical to before subcommands were introduced.
#[derive(Debug, Default)]
pub(crate) struct RunConfig {
    pub(crate) algorithm: MatchingAlgorithm,
    pub(crate) alignment: Vec<String>,
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
    pub(crate) name_encoder: NameEncoderKind,
}

impl RunConfig {
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
        }

        // CRAM sanity: any .cram input requires --reference.
        if self.alignment.iter().any(|p| p.ends_with(".cram"))
            && self.io.reference.is_none()
        {
            return Err(Error::CramRequiresReference);
        }

        // stdin: at most one stream, namesorted only.
        let stdin_count =
            self.alignment.iter().filter(|p| p.as_str() == "-").count();
        if stdin_count > 1 {
            return Err(Error::MultipleStdinStreams);
        }
        if stdin_count == 1
            && self.algorithm != MatchingAlgorithm::Namesorted
        {
            return Err(Error::StdinNamesortedOnly);
        }

        Ok(())
    }
}
