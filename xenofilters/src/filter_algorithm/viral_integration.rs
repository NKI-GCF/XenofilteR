use crate::{
    Error,
    config::args::{IoArgs, RelatedArgs, ScoringArgs},
    config::run_config::RunConfig,
};
use clap::Args;

/// Preconfigured chimeric-pair detection for viral integration studies.
/// Requires exactly 2 or 3 alignment streams: host, virus, [optional xenograft host].
#[derive(Args, Clone, Debug)]
pub(crate) struct ViralIntegrationArgs {
    /// Labels for host/virus[/xenograft] streams. REQUIRED -- used in XC:Z tags.
    #[arg(long, required = true, num_args = 2..=3, help_heading = "Chimeric")]
    pub(crate) stream_labels: Vec<String>,

    #[command(flatten)]
    pub(crate) io: IoArgs,

    #[command(flatten)]
    pub(crate) scoring: ScoringArgs,

    #[command(flatten)]
    pub(crate) variants: RelatedArgs,

    #[arg(short = 'S', long, default_value = "1", help_heading = "Parallelism")]
    pub(crate) score_threads: usize,
}

impl ViralIntegrationArgs {
    pub(crate) fn into_run_config(self) -> Result<RunConfig, Error> {
        Ok(RunConfig {
            io: self.io,
            scoring: self.scoring,
            variants: self.variants,
            ..Default::default()
        })
    }
    pub(crate) fn validate_and_init(&mut self) -> Result<(), Error> {
        let got = self.io.alignment.len();
        if !(2..=3).contains(&got) {
            return Err(Error::ViralIntegrationStreamCount { got });
        }
        if self.stream_labels.len() != got {
            return Err(Error::ViralIntegrationLabelCount {
                alignments: got,
                labels: self.stream_labels.len(),
            });
        }
        Ok(())
    }
}
