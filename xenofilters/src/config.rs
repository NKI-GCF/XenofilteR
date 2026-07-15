use crate::{
    filter_algorithm::line_by_line::{core::COUNTER_STRIDE, MAX_STREAMS},
    Error,
    regions::ScoreFn,
    bam::AlnFormat,
    penalty::{ErrorModel, Penalty},
};
use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

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
#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Default)]
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
}

// -- Top-level CLI --------------------------------------------------------------

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

#[derive(Subcommand, Debug)]
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
}

impl AlgorithmCommand {
    /// Extract the common args regardless of which subcommand was chosen.
    pub(crate) fn common(&self) -> &CommonArgs {
        match self {
            AlgorithmCommand::Namesorted(a) => &a.common,
            AlgorithmCommand::Hashlookup(a) => &a.common,
            AlgorithmCommand::Collated(a) => &a.common,
        }
    }
    pub(crate) fn common_mut(&mut self) -> &mut CommonArgs {
        match self {
            AlgorithmCommand::Namesorted(a) => &mut a.common,
            AlgorithmCommand::Hashlookup(a) => &mut a.common,
            AlgorithmCommand::Collated(a) => &mut a.common,
        }
    }
}

// ── Shared argument groups (flattened into each subcommand) ──────────────────

/// Arguments shared by all three algorithms.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct CommonArgs {
    // ── Input ────────────────────────────────────────────────────────────────
    /// Input alignment files (BAM, CRAM, or - for stdin BAM).
    #[arg(required = true, num_args = 1..=32, help_heading = "Input")]
    pub(crate) alignment: Vec<String>,

    /// Reference FASTA for CRAM decoding. Required when any input is CRAM.
    #[arg(long, help_heading = "Input")]
    pub(crate) reference: Option<PathBuf>,

    /// Strip read-name suffix (/1, /2). auto | true | false | variable
    #[arg(short = 'R', long, default_value = "auto", help_heading = "Input")]
    pub(crate) strip_read_suffix: StripReadSuffix,

    // ── Output ───────────────────────────────────────────────────────────────
    /// Winner output, one file per stream.
    #[arg(short = 'o', long, num_args = 1..=32, help_heading = "Output")]
    pub(crate) output: Vec<PathBuf>,

    /// Ambiguous and discarded output, one file per stream.
    /// Both ambiguous and discarded reads go here; tagged _xenoambig /
    /// _xenodiscard in RG:Z when configured.
    #[arg(short = 'a', long, num_args = 0..=32, help_heading = "Output")]
    pub(crate) ambiguous_output: Vec<PathBuf>,

    /// Human-readable labels for each stream (positional: stream 0, 1, …).
    /// Used in XC:Z chimeric tags and JSON stats output.
    #[arg(long, num_args = 0..=32, help_heading = "Output")]
    pub(crate) stream_labels: Vec<String>,

    /// Write JSON summary statistics (MultiQC-compatible).
    #[arg(long, help_heading = "Output")]
    pub(crate) stats_output: Option<PathBuf>,

    /// Add XF:C (confidence) or XR:C (variant-rescued) aux tags to records.
    #[arg(short = 'A', long, default_value = "false", help_heading = "Output")]
    pub(crate) add_decision_tag: bool,

    /// Suppress @PG header line in output BAMs.
    #[arg(short = 'P', long, default_value = "false", help_heading = "Output")]
    pub(crate) no_program_line: bool,

    /// Write JSON summary to this path (MultiQC generalstats format).
    #[arg(long, env = "XENOFILTERS_STATS_OUTPUT", help_heading = "Output")]
    pub(crate) stats_path: Option<PathBuf>,

    // ── Scoring ───────────────────────────────────────────────────────────────
    /// Platform error model. Sets gap/quality defaults.
    /// illumina (default) | hifi | ont
    #[arg(long, default_value = "illumina", help_heading = "Scoring")]
    pub(crate) error_model: ErrorModel,

    /// Mismatch penalty (positive; internally negated).
    #[arg(short = 'm', long, default_value = "4.0", help_heading = "Scoring")]
    pub(crate) mismatch_penalty: f64,

    /// Gap-open penalty (positive; internally negated).
    #[arg(short = 'g', long, default_value = "6.0", help_heading = "Scoring")]
    pub(crate) gap_open: f64,

    /// Gap-extend penalty (positive; internally negated).
    #[arg(short = 'e', long, default_value = "1.0", help_heading = "Scoring")]
    pub(crate) gap_extend: f64,

    /// Chimeric-junction structural penalty constant (bases).
    #[arg(short = 'J', long, default_value = "20",
          value_parser = clap::value_parser!(u32).range(0..=10000),
          help_heading = "Scoring")]
    pub(crate) chimeric_junction_bases: u32,

    /// Ambiguous threshold (Phred). u32::MAX = auto (10 for pass 1, 0 for pass 2).
    #[arg(long, default_value_t = u32::MAX, value_name = "PHRED|auto",
          help_heading = "Scoring")]
    pub(crate) ambiguous_threshold: u32,

    /// Warn when ambiguous fraction exceeds this fraction. Default 0.05 (5%).
    #[arg(long, default_value = "0.05", help_heading = "Scoring")]
    pub(crate) warn_ambig_fraction: f64,

    /// Enable bisulfite/WGBS scoring: C→T (forward) and G→A (reverse)
    /// mismatches scored as zero-penalty conversions.
    #[arg(long, default_value = "false", help_heading = "Scoring")]
    pub(crate) bisulfite: bool,

    // ── Variants ──────────────────────────────────────────────────────────────
    /// Sample VCF/BCF per stream ([IDX:]FILE). FORMAT/GT + FORMAT/GQ required.
    #[arg(short = 's', long, num_args = 0..=32,
          value_name = "[IDX:]FILE", help_heading = "Variants")]
    pub(crate) sample_variants: Vec<String>,

    /// Population VCF/BCF per stream ([IDX:]FILE). INFO/AF required.
    #[arg(short = 'p', long, num_args = 0..=32,
          value_name = "[IDX:]FILE", help_heading = "Variants")]
    pub(crate) population_variants: Vec<String>,

    // ── Misc ──────────────────────────────────────────────────────────────────
    /// Suppress progress output to stderr.
    #[arg(long, default_value = "false", help_heading = "Misc")]
    pub(crate) quiet: bool,

    /// Verbosity: -v = INFO, -vv = DEBUG. Respects RUST_LOG.
    #[arg(short, long, action = clap::ArgAction::Count, help_heading = "Misc")]
    pub(crate) verbose: u8,

    /// Write unmapped reads to output. For testing only.
    #[arg(long, default_value = "false", hide = true)]
    pub(crate) write_unmapped: bool,

    // ── Internal (skip) ───────────────────────────────────────────────────────
    #[arg(skip)]
    pub(crate) is_pass2: bool,

    #[arg(skip)]
    pub(crate) parsed_chimeric_pairs: Vec<[usize; 2]>,
}

// ── Per-algorithm arg structs ─────────────────────────────────────────────────

#[derive(Args, Debug)]
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

    /// Parallel fragment-scoring worker threads (0 = all logical CPUs, max 16).
    #[arg(
        short = 'S',
        long,
        default_value = "1",
        env = "XENOFILTERS_SCORE_THREADS",
        help_heading = "Parallelism"
    )]
    pub(crate) score_threads: usize,

    /// Chimeric stream-index pairs (format: A:B).
    /// Reads spanning species boundaries are tagged XC:Z:<other_label>
    /// and written to both streams' outputs.
    #[arg(long, num_args = 0.., value_name = "A:B",
          help_heading = "Chimeric")]
    pub(crate) chimeric_pairs: Vec<String>,

    /// Name encoder for FragmentTable key compression.
    /// illumina (strips machine:run:flowcell prefix) | passthrough
    #[arg(long, default_value = "illumina", help_heading = "Advanced")]
    pub(crate) name_encoder: NameEncoderKind,
}

#[derive(Args, Debug)]
pub(crate) struct HashlookupArgs {
    #[command(flatten)]
    pub(crate) common: CommonArgs,

    /// In-memory BED file(s) forcing full NW scoring on overlap.
    /// One per stream (positional: stream 0, then 1).
    /// Strand-aware when BED column 6 present.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions")]
    pub(crate) ambiguous_regions: Vec<String>,

    /// In-memory VCF/BCF of diagnostic positions per stream.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions")]
    pub(crate) diagnostic_variants: Vec<String>,

    /// BED file(s) of regions giving reads a positive score bonus.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions")]
    pub(crate) positive_regions: Vec<String>,

    /// Score function for --positive-regions BED score column.
    /// Format: fn[:weight]  fn ∈ {linear, log, constant, overlap_fraction}
    #[arg(long, default_value = "linear:1.0", help_heading = "Regions")]
    pub(crate) region_score_fn: ScoreFn,

    /// Name encoder for key compression.
    #[arg(long, default_value = "illumina", help_heading = "Advanced")]
    pub(crate) name_encoder: NameEncoderKind,
}

#[derive(Args, Debug)]
pub(crate) struct CollatedArgs {
    #[command(flatten)]
    pub(crate) common: CommonArgs,

    /// Tabix-indexed BED.gz file(s) forcing full NW on overlap.
    /// One per stream. Strand-aware when BED column 6 present.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions")]
    pub(crate) ambiguous_regions: Vec<String>,

    /// Tabix-indexed VCF/BCF of diagnostic positions per stream.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions")]
    pub(crate) diagnostic_variants: Vec<String>,

    /// BED.gz file(s) of positive-score regions. Tabix-indexed.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions")]
    pub(crate) positive_regions: Vec<String>,

    #[arg(long, default_value = "linear:1.0", help_heading = "Regions")]
    pub(crate) region_score_fn: ScoreFn,
}

#[cfg(test)]
pub(crate) mod tests;
