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

// src/config/run_config.rs

impl RunConfig {
    pub(crate) fn validate_and_init(&mut self) -> Result<(), Error> {
        self.io.validate()?;
        self.scoring.validate()?;
        let n = self
            .alignment
            .len()
            .max(if self.single_alignment_mode { 2 } else { 0 });

        match self.algorithm {
            MatchingAlgorithm::Namesorted => {
                ensure!(
                    n >= 1 && n <= MAX_STREAMS,
                    "namesorted supports 1..={MAX_STREAMS} streams (got {n})"
                );
                self.output.validate(n)?;
            }
            MatchingAlgorithm::Hashlookup | MatchingAlgorithm::Collated => {
                ensure!(
                    self.alignment.len() == 2,
                    "{:?} requires exactly 2 streams (got {})",
                    self.algorithm,
                    self.alignment.len()
                );
                ensure!(
                    self.score_threads == 1,
                    "{:?} does not support --score-threads > 1",
                    self.algorithm
                );
                self.output.validate(2)?;
            }
        }

        // CRAM sanity: any .cram input requires --reference.
        if self.alignment.iter().any(|p| p.ends_with(".cram")) {
            ensure!(
                self.io.reference.is_some(),
                "CRAM input requires --reference <fasta>"
            );
        }

        // stdin: at most one stream, namesorted only.
        let stdin_count = self.alignment.iter().filter(|p| p.as_str() == "-").count();
        ensure!(
            stdin_count <= 1,
            "at most one stream may be read from stdin"
        );
        if stdin_count == 1 {
            ensure!(
                self.algorithm == MatchingAlgorithm::Namesorted,
                "stdin input is only supported with namesorted"
            );
        }

        Ok(())
    }
}
