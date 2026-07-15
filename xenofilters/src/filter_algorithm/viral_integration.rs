use clap::Args;
use crate::{
    config::args::{IoArgs, ScoringArgs, VariantArgs, OutputArgsMulti},
    config::run_config::RunConfig,
    config::MatchingAlgorithm,
    Error,
};

/// Preconfigured chimeric-pair detection for viral integration studies.
/// Requires exactly 2 or 3 alignment streams: host, virus, [optional xenograft host].
#[derive(Args, Clone, Debug)]
pub(crate) struct ViralIntegrationArgs {
    /// Host and viral reference alignments (2 streams), or host, virus,
    /// xenograft-host (3 streams, e.g. HPV+human tumour grafted in mouse).
    #[arg(required = true, num_args = 2..=3, value_hint = clap::ValueHint::FilePath,
          help_heading = "Input")]
    pub(crate) alignment: Vec<String>,

    /// Labels for host/virus[/xenograft] streams. REQUIRED — used in XC:Z tags.
    #[arg(long, required = true, num_args = 2..=3, help_heading = "Chimeric")]
    pub(crate) stream_labels: Vec<String>,

    #[command(flatten)]
    pub(crate) io: IoArgs,

    #[command(flatten)]
    pub(crate) scoring: ScoringArgs,

    #[command(flatten)]
    pub(crate) variants: VariantArgs,

    #[command(flatten)]
    pub(crate) output: OutputArgsMulti,

    #[arg(short = 't', long, default_value = "4", help_heading = "Parallelism")]
    pub(crate) threads: usize,

    #[arg(short = 'S', long, default_value = "1", help_heading = "Parallelism")]
    pub(crate) score_threads: usize,
}

impl ViralIntegrationArgs {
    pub(crate) fn into_run_config(self) -> Result<RunConfig, Error> {
        Ok(RunConfig {
            algorithm: MatchingAlgorithm::Namesorted,
            alignment: self.alignment,
            io: self.io, scoring: self.scoring, variants: self.variants,
            output: self.output,
            threads: self.threads, score_threads: self.score_threads,
            // Preset: streams 0 and 1 are always the host↔virus chimeric pair.
            chimeric_pairs: vec!["0:1".to_string()],
            stream_labels: self.stream_labels,
            ..Default::default()
        })
    }
}
