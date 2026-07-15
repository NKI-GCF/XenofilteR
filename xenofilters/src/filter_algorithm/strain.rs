/// Within-species strain disambiguation: single alignment, two variant
/// profiles (one per strain). No chimeric pairs, no multi-stream regions,
/// no parallelism flags beyond bgzf threads (NW scoring is inherently
/// single-pair here; --score-threads would only help with many fragments
/// in flight, which the underlying engine already supports — exposed here
/// as --threads only, kept simple).
#[derive(Args, Debug)]
pub(crate) struct StrainArgs {
    /// Single alignment BAM/CRAM (duplicated internally into two logical
    /// streams, one scored against each strain's variant profile).
    #[arg(required = true, value_hint = clap::ValueHint::FilePath,
          help_heading = "Input")]
    pub(crate) alignment: String,

    #[command(flatten)]
    pub(crate) io: IoArgs,

    #[command(flatten)]
    pub(crate) scoring: ScoringArgs,

    /// Strain A / strain B variant profiles. At least one of
    /// --sample-variants or --population-variants must supply both
    /// indices 0 and 1 (one profile per strain).
    #[command(flatten)]
    pub(crate) variants: VariantArgs,

    #[command(flatten)]
    pub(crate) output: OutputArgsPair,

    #[arg(short = 't', long, default_value = "4", help_heading = "Parallelism")]
    pub(crate) threads: usize,
}

impl StrainArgs {
    pub(crate) fn into_run_config(self) -> Result<RunConfig, Error> {
        ensure!(
            has_profile_for(&self.variants, 0) && has_profile_for(&self.variants, 1),
            "strain: requires a variant profile for both strain 0 and strain 1 \
             (via --sample-variants or --population-variants, e.g. '0:a.vcf 1:b.vcf')"
        );
        Ok(RunConfig {
            algorithm: MatchingAlgorithm::Namesorted,
            alignment: vec![self.alignment.clone(), self.alignment], // duplicated stream
            single_alignment_mode: true,
            io: self.io,
            scoring: self.scoring,
            variants: self.variants,
            output: self.output.into(),
            threads: self.threads,
            score_threads: 1,
            chimeric_pairs: vec![],
            stream_labels: vec![],
            ..Default::default()
        })
    }
}
