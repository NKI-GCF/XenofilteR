pub(crate) mod args;
pub(crate) mod run_config;

use crate::{
    filter_algorithm::{strain::StrainArgs, viral_integration::ViralIntegrationArgs},
    region::ScoreFn,
};
use clap::{Args, Parser, Subcommand, ValueEnum};
use clap_complete::Shell;

const ARG_MAX: usize = 4;

#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Default)]
pub enum StripReadSuffix {
    #[default]
    Auto,
    True,
    False,
    Variable,
}

/// Fragment-matching algorithm.
#[derive(Subcommand, Debug, Default, PartialEq)]
pub enum MatchingAlgorithm {
    /// Streaming merge. Both BAMs must be in identical query-name order.
    /// Lowest memory. Fastest. Supports multithreading. Default.
    #[default]
    Namesorted,

    /// Hash-table matching. Each BAM must be collated (all records for a read
    /// name contiguous) but the two BAMs may present fragments in different
    /// orders. Memory proportional to name-order skew. Single-threaded only.
    /// Output order is not guaranteed.
    Collated,

    /// Hash-table matching. Works on arbitrary (non-sorted, non-collated) BAM
    /// input. High memory usage — proportional to the number of in-flight
    /// fragments. Single-threaded only. Preserves driving-stream (stream 0)
    /// output order.
    Hashlookup,

    /// Within-species strain/genotype disambiguation from a single alignment.
    ///
    /// Compares one BAM against two strain-specific variant profiles
    /// (e.g. C57BL/6 vs CAST/EiJ) rather than two separate reference
    /// alignments. Equivalent to `namesorted --single-alignment-mode`
    /// with a restricted flag surface.
    Strain,

    /// Cross-species viral integration detection (HPV, HBV, HTLV-1, …).
    ///
    /// Preconfigures chimeric-pair detection between the first two streams.
    /// Equivalent to `namesorted --chimeric-pairs 0:1` with --stream-labels
    /// required and output arity restricted to the streams involved.
    ViralIntegration,

}

// -- Top-level CLI -------------------------------------------------------------

/// Fast multi-stream read classifier for xenograft, PDX, viral integration,
/// and cross-species sequencing experiments.
#[derive(Parser, Debug, Default, Clone)]
#[command(
    author, version,
    about = "Disambiguate reads from multiple sequence alignments",
    long_about = None,
    subcommand_required = true,
    arg_required_else_help = true,
)]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub(crate) command: AlgorithmCommand,
}

#[derive(Copy, Clone, Debug, ValueEnum, Default)]
pub(crate) enum NameEncoderKind {
    #[default]
    Illumina,     // strip machine:run:flowcell prefix
    Passthrough,  // full name (PacBio, ONT, custom)
}

#[derive(Subcommand, Clone, Debug)]
pub(crate) enum AlgorithmCommand {
    /// Disambiguate name-sorted BAM/CRAM streams (default for most pipelines).
    ///
    /// All input files must be sorted by query name in identical order.
    /// Supports up to 32 simultaneous alignment streams, parallel scoring,
    /// chimeric-pair detection, and stdin BAM input.
    Namesorted(NamesortedArgs),

    /// Disambiguate coordinate-sorted BAMs without re-sorting.
    ///
    /// Two-pass algorithm: pass 1 indexes fragment names via BGZF virtual
    /// offsets; pass 2 retrieves and emits records in driving-stream order.
    /// Requires exactly 2 BAM streams. Single-threaded. In-memory BED/VCF
    /// region files force full NW scoring on overlapping reads.
    Hashlookup(HashlookupArgs),

    /// Disambiguate collated BAM streams (records grouped by name, any order).
    ///
    /// Each input BAM must be collated (all records for a fragment contiguous)
    /// but the two streams need not present fragments in the same order.
    /// Requires exactly 2 streams. Supports tabix-indexed BED/VCF regions.
    Collated(CollatedArgs),

    /// Within-species strain/genotype disambiguation from a single alignment.
    ///
    /// Compares one BAM against two strain-specific variant profiles
    /// (e.g. C57BL/6 vs CAST/EiJ) rather than two separate reference
    /// alignments. Equivalent to `namesorted --single-alignment-mode`
    /// with a restricted flag surface.
    Strain(StrainArgs),

    /// Cross-species viral integration detection (HPV, HBV, HTLV-1, …).
    ///
    /// Preconfigures chimeric-pair detection between the first two streams.
    /// Equivalent to `namesorted --chimeric-pairs 0:1` with --stream-labels
    /// required and output arity restricted to the streams involved.
    ViralIntegration(ViralIntegrationArgs),

    /// Generate shell completion script and print to stdout.
    ///
    /// Pipe the output to the appropriate file for your shell:
    ///
    ///   bash:   xenofilters completion bash > ~/.bash_completion.d/xenofilters
    ///   zsh:    xenofilters completion zsh  > ~/.zfunc/_xenofilters
    ///   fish:   xenofilters completion fish > ~/.config/fish/completions/xenofilters.fish
    ///   elvish: xenofilters completion elvish > ~/.config/elvish/completions/xenofilters.elv
    ///
    /// For zsh, ensure 'autoload -Uz compinit && compinit' is in ~/.zshrc
    /// and that ~/.zfunc is in $fpath.
    #[command(name = "completion")]
    Completion(CompletionArgs),
}

impl AlgorithmCommand {
    /// Extract the common args regardless of which subcommand was chosen.
    pub(crate) fn common(&self) -> &CommonArgs {
        match self {
            AlgorithmCommand::Namesorted(a) => &a.common,
            AlgorithmCommand::Hashlookup(a) => &a.common,
            AlgorithmCommand::Collated(a) => &a.common,
            AlgorithmCommand::Strain(a) => &a.common,
            AlgorithmCommand::ViralIntegration(a) => &a.common,
        }
    }
    pub(crate) fn common_mut(&mut self) -> &mut CommonArgs {
        match self {
            AlgorithmCommand::Namesorted(a) => &mut a.common,
            AlgorithmCommand::Hashlookup(a) => &mut a.common,
            AlgorithmCommand::Collated(a) => &mut a.common,
            AlgorithmCommand::Strain(a) => &mut a.common,
            AlgorithmCommand::ViralIntegration(a) => &mut a.common,
        }
    }
}

impl Default for AlgorithmCommand {
    fn default() -> Self {
        AlgorithmCommand::Namesorted(NamesortedArgs {
            common: CommonArgs::default(),
            threads: 4,
            name_encoder: NameEncoderKind::Illumina,
            ..Default::default()
        })
    }
}

// -- Shared argument groups (flattened into each subcommand) ------------------

/// Arguments shared by all three algorithms.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct CommonArgs {
    #[command(flatten)]
    pub(crate) io: crate::config::args::IoArgs,
    #[command(flatten)]
    pub(crate) scoring: crate::config::args::ScoringArgs,
    #[command(flatten)]
    pub(crate) variants: crate::config::args::VariantArgs,
    #[command(flatten)]
    pub(crate) output: crate::config::args::OutputArgsMulti,
}

// -- Per-algorithm arg structs -------------------------------------------------

#[derive(Args, Clone, Debug)]
pub(crate) struct NamesortedArgs {
    #[command(flatten)]
    pub(crate) common: CommonArgs,

    /// bgzf decompression worker threads.
    #[arg(
        short = 't',
        long,
        default_value = "4",
        env = "XENOFILTERS_THREADS",
        help_heading = "Parallelism"
    )]
    pub(crate) threads: usize,

    #[command(flatten)]
    pub(crate) parallel: crate::config::args::ParallelArgs,

    #[command(flatten)]
    pub(crate) chimeric: crate::config::args::ChimericArgs,

    /// Name encoder for FragmentTable key compression.
    /// illumina (strips machine:run:flowcell prefix) | passthrough
    #[arg(long, default_value = "illumina", help_heading = "Advanced")]
    pub(crate) name_encoder: NameEncoderKind,
}

#[derive(Args, Clone, Debug)]
pub(crate) struct HashlookupArgs {
    #[command(flatten)]
    pub(crate) common: CommonArgs,

    /// In-memory BED file(s) forcing full NW scoring on overlap.
    /// One per stream (positional: stream 0, then 1).
    /// Strand-aware when BED column 6 present.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) ambiguous_regions: Vec<String>,

    /// In-memory VCF/BCF of diagnostic positions per stream.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) diagnostic_variants: Vec<String>,

    /// BED file(s) of regions giving reads a positive score bonus.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) positive_regions: Vec<String>,

    /// Score function for --positive-regions BED score column.
    /// Format: fn[:weight]  fn ∈ {linear, log, constant, overlap_fraction}
    #[arg(long, default_value = "linear:1.0", help_heading = "Regions")]
    pub(crate) region_score_fn: ScoreFn,

    /// Name encoder for key compression.
    #[arg(long, default_value = "illumina", help_heading = "Advanced")]
    pub(crate) name_encoder: NameEncoderKind,
}

#[derive(Args, Clone, Debug)]
pub(crate) struct CollatedArgs {
    #[command(flatten)]
    pub(crate) common: CommonArgs,

    /// Tabix-indexed BED.gz file(s) forcing full NW on overlap.
    /// One per stream. Strand-aware when BED column 6 present.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) ambiguous_regions: Vec<String>,

    /// Tabix-indexed VCF/BCF of diagnostic positions per stream.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) diagnostic_variants: Vec<String>,

    /// BED.gz file(s) of positive-score regions. Tabix-indexed.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) positive_regions: Vec<String>,

    #[arg(long, default_value = "linear:1.0", help_heading = "Regions")]
    pub(crate) region_score_fn: ScoreFn,
}

#[derive(Args, Debug, Clone, PartialEq)]
pub(crate) struct CompletionArgs {
    /// Shell to generate completions for.
    #[arg(value_enum)]
    pub(crate) shell: ShellChoice,
}

/// Supported shells. Wraps `clap_complete::Shell` plus Nu.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum ShellChoice {
    Bash,
    Zsh,
    Fish,
    Elvish,
    PowerShell,
}

impl ShellChoice {
    pub(crate) fn to_clap_shell(self) -> Option<Shell> {
        match self {
            ShellChoice::Bash => Some(Shell::Bash),
            ShellChoice::Zsh => Some(Shell::Zsh),
            ShellChoice::Fish => Some(Shell::Fish),
            ShellChoice::Elvish => Some(Shell::Elvish),
            ShellChoice::PowerShell => Some(Shell::PowerShell),
        }
    }
}

impl Default for NamesortedArgs {
    fn default() -> Self {
        Self {
            common: CommonArgs::default(),
            threads: 4,
            parallel: crate::config::args::ParallelArgs::default(),
            chimeric: crate::config::args::ChimericArgs::default(),
            name_encoder: NameEncoderKind::Illumina,
        }
    }
}

#[cfg(test)]
pub(crate) mod tests;
