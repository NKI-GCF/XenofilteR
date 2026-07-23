use crate::{
    Error,
    config::args::{IoArgs, RelatedArgs, ScoringArgs},
    config::run_config::RunConfig,
};
use clap::Args;

/// Within-species strain disambiguation: single alignment, two variant
/// profiles (one per strain). No chimeric pairs, no multi-stream regions,
/// no parallelism flags beyond bgzf threads (NW scoring is inherently
/// single-pair here; --score-threads would only help with many fragments
/// in flight, which the underlying engine already supports -- exposed here
/// as --threads only, kept simple).
#[derive(Args, Clone, Debug)]
pub(crate) struct StrainArgs {
    #[command(flatten)]
    pub(crate) io: IoArgs,

    #[command(flatten)]
    pub(crate) scoring: ScoringArgs,

    /// Strain A / strain B variant profiles. At least one of
    /// --sample-variants or --population-variants must supply both
    /// indices 0 and 1 (one profile per strain).
    #[command(flatten)]
    pub(crate) variants: RelatedArgs,

    #[command(flatten)]
    pub(crate) parallel: crate::config::args::ParallelArgs,
}

impl StrainArgs {
    pub(crate) fn into_run_config(self) -> Result<RunConfig, Error> {
        if !self.variants.has_index(0) || !self.variants.has_index(1) {
            return Err(Error::StrainMissingVariantProfile);
        }
        let mut io = self.io;
        let alignment = io.alignment[0].clone();
        io.alignment.push(alignment);
        Ok(RunConfig {
            io,
            scoring: self.scoring,
            variants: self.variants,
            chimeric_pairs: vec![],
            stream_labels: vec![],
            ..Default::default()
        })
    }
    pub(crate) fn validate_and_init(&mut self) -> Result<(), Error> {
        let got = self.io.alignment.len();
        if got != 1 {
            return Err(Error::StrainStreamCount { got });
        }
        if !self.variants.has_index(0) || !self.variants.has_index(1) {
            return Err(Error::StrainMissingVariantProfile);
        }
        Ok(())
    }
}
