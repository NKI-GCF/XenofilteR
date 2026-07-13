pub(crate) mod args;
pub(crate) mod run_config;

use crate::{
    bam::AlnFormat,
    filter_algorithm::line_by_line::{core::COUNTER_STRIDE, MAX_STREAMS},
    penalty::{ErrorModel, Penalty},
    regions::ScoreFn,
    Error,
};
use clap::{Args, Parser, Subcommand, ValueEnum};
use clap_complete::Shell;
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
#[derive(Clone, Debug, ValueEnum, PartialEq, Default)]
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

// -- Shared argument groups (flattened into each subcommand) ------------------

/// Arguments shared by all three algorithms.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct CommonArgs {
    // -- Input ----------------------------------------------------------------
    /// Input alignment files (BAM, CRAM, or - for stdin BAM).
    #[arg(required = true, num_args = 1..=32, help_heading = "Input", value_hint = clap::ValueHint::FilePath)]
    pub(crate) alignment: Vec<PathBuf>,

    /// Reference FASTA for CRAM decoding. Required when any input is CRAM.
    #[arg(long, help_heading = "Input", value_hint = clap::ValueHint::FilePath)]
    pub(crate) reference: Option<PathBuf>,

    /// Strip read-name suffix (/1, /2). auto | true | false | variable
    #[arg(short = 'R', long, default_value = "auto", help_heading = "Input")]
    pub(crate) strip_read_suffix: StripReadSuffix,

    // -- Output ---------------------------------------------------------------
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
    pub(crate) stream_labels: Vec<PathBuf>,

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

    // -- Scoring ---------------------------------------------------------------
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

    // -- Variants --------------------------------------------------------------
    /// Sample VCF/BCF per stream ([IDX:]FILE). FORMAT/GT + FORMAT/GQ required.
    #[arg(short = 's', long, num_args = 0..=32,
          value_name = "[IDX:]FILE", help_heading = "Variants", value_hint = clap::ValueHint::FilePath)]
    pub(crate) sample_variants: Vec<PathBuf>,

    /// Population VCF/BCF per stream ([IDX:]FILE). INFO/AF required.
    #[arg(short = 'p', long, num_args = 0..=32,
          value_name = "[IDX:]FILE", help_heading = "Variants", value_hint = clap::ValueHint::FilePath)]
    pub(crate) population_variants: Vec<PathBuf>,

    // -- Misc ------------------------------------------------------------------
    /// Suppress progress output to stderr.
    #[arg(long, default_value = "false", help_heading = "Misc")]
    pub(crate) quiet: bool,

    /// Verbosity: -v = INFO, -vv = DEBUG. Respects RUST_LOG.
    #[arg(short, long, action = clap::ArgAction::Count, help_heading = "Misc")]
    pub(crate) verbose: u8,

    /// Write unmapped reads to output. For testing only.
    #[arg(long, default_value = "false", hide = true)]
    pub(crate) write_unmapped: bool,

    // -- Internal (skip) -------------------------------------------------------
    #[arg(skip)]
    pub(crate) is_pass2: bool,

    #[arg(skip)]
    pub(crate) parsed_chimeric_pairs: Vec<[usize; 2]>,
}

impl CommonArgs {
    pub(super) fn validate_and_init(&mut self) -> Result<(), Error> {
        // Detect pass 2: any input BAM whose header contains _xenoambig RG.
        // Detection is done at stream-open time; store result for threshold selection.
        // self.is_pass2 is set in AlnStream::new() after reading the header.

        // Resolve threshold.
        self.resolved_ambiguous_log_threshold = self.ambiguous_threshold_to_log_likelihood();

        if self.is_pass2 {
            tracing::info!(
                threshold_phred = match self.ambiguous_threshold {
                    u32::MAX =>
                        if self.is_pass2 {
                            0
                        } else {
                            10
                        },
                    phred => phred,
                },
                "Pass 2 detected via _xenoambig read groups"
            );
        }

        let aln_count = self.alignment.len();
        // -- Chimeric pair parsing --------------------------------------------
        let mut parsed_chimeric_pairs: Vec<[usize; 2]> = Vec::new();
        for raw in &self.chimeric_pairs {
            let (a_str, b_str) = raw
                .split_once(':')
                .ok_or(Error::ChimericPairsInvalidFormat { raw: raw.clone() })?;
            let a =
                a_str
                    .trim()
                    .parse::<usize>()
                    .map_err(|_| Error::ChimericPairsInvalidIndex {
                        index_str: a_str.to_string(),
                    })?;
            let b =
                b_str
                    .trim()
                    .parse::<usize>()
                    .map_err(|_| Error::ChimericPairsInvalidIndex {
                        index_str: b_str.to_string(),
                    })?;

            if a == b {
                return Err(Error::ChimericPairsIdenticalIndices { raw: raw.clone() });
            }
            if self.matching_algorithm != MatchingAlgorithm::Namesorted {
                return Err(Error::ChimericPairsRequiresNamesorted);
            }
            // Canonical order: lower index first.
            parsed_chimeric_pairs.push([a.min(b), a.max(b)]);
        }
        // Deduplicate.
        parsed_chimeric_pairs.sort_unstable();
        parsed_chimeric_pairs.dedup();
        self.parsed_chimeric_pairs = parsed_chimeric_pairs;

        if !self.parsed_chimeric_pairs.is_empty() {
            tracing::info!(
                pairs = ?self.parsed_chimeric_pairs,
                "Chimeric pair detection enabled"
            );
        }
        if self.warn_ambig_fraction < 0.0 || self.warn_ambig_fraction > 1.0 {
            return Err(Error::WarnAmbigFractionOutOfRange {
                value: self.warn_ambig_fraction,
            });
        }

        if self.input_format == AlnFormat::Cram && self.reference.is_none() {
            return Err(Error::CramRequiresReference);
        }
        if self.input_format == AlnFormat::Cram
            && self.matching_algorithm != MatchingAlgorithm::Namesorted
        {
            return Err(Error::CramNamesortedOnly);
        }
        // Stdin: path "-" implies SAM unless explicitly overridden.
        if self.alignment.iter().any(|p| p == "-") {
            if self.alignment.iter().filter(|p| p.as_str() == "-").count() > 1 {
                return Err(Error::MultipleStdinStreams);
            }
            if self.matching_algorithm != MatchingAlgorithm::Namesorted {
                return Err(Error::StdinRequiresNamesorted);
            }
        }

        // Reject multi-threaded modes for non-namesorted algorithms.
        if self.matching_algorithm == MatchingAlgorithm::Namesorted {
            if !self.ambiguous_regions.is_empty() || !self.diagnostic_variants.is_empty() {
                return Err(Error::NamesortedUnsupportedOptions);
            }
            // Resolve score_threads = 0 → available logical CPUs with a max of 16.
            if self.score_threads == 0 {
                self.score_threads = std::thread::available_parallelism()
                    .map(|n| n.get().min(16))
                    .unwrap_or(1);
                tracing::info!(
                    score_threads = self.score_threads,
                    "score_threads=0: using all available logical CPUs"
                );
            }
        } else {
            if self.score_threads != 1 {
                return Err(Error::MultiThreadedScoringRequiresNamesorted);
            }
        }

        // 1. Guard against single-stream niche execution accidents first
        if aln_count == 1 {
            if !self.single_alignment_mode {
                return Err(Error::SingleStreamMissingFlag);
            }
            if self.read_from_stdin {
                return Err(Error::SingleStreamStdinUnsupported);
            }
            if self.matching_algorithm != MatchingAlgorithm::Namesorted {
                return Err(Error::SingleStreamRequiresNamesorted);
            }
        } else {
            if self.single_alignment_mode {
                return Err(Error::SingleAlignmentModeExpectsOneStream { count: aln_count });
            }
            if aln_count < 2 {
                return Err(Error::InsufficientAlignmentStreams { count: aln_count });
            }
            if aln_count > 2 {
                if self.matching_algorithm != MatchingAlgorithm::Namesorted {
                    return Err(Error::MultiStreamRequiresNamesorted);
                }
                if aln_count > MAX_STREAMS {
                    return Err(Error::MaxStreamsExceeded {
                        count: aln_count,
                        max: MAX_STREAMS,
                    });
                }
            }
        }

        if self.ambiguous_regions.len() > 2 {
            return Err(Error::TooManyAmbiguousRegionsFiles {
                count: self.ambiguous_regions.len(),
            });
        }
        if self.diagnostic_variants.len() > 2 {
            return Err(Error::TooManyDiagnosticVariantsFiles {
                count: self.diagnostic_variants.len(),
            });
        }
        // Determine effective dimensions (logical comparisons)
        let logical_len = if aln_count == 1 { 2 } else { aln_count };

        // 2. Parse, validate, and normalize variant arrays to size matching comparisons
        let mut normalized_samples = vec![PathBuf::new(); logical_len];
        let mut normalized_populations = vec![PathBuf::new(); logical_len];

        let mut stream_has_variants = vec![false; logical_len];

        for (i, arg) in self.sample_variants.iter().enumerate() {
            let (idx, path) = Self::parse_variant_string(arg, i)?;
            if idx >= logical_len {
                return Err(Error::SampleVariantIndexOutOfBounds {
                    idx,
                    max: logical_len,
                });
            }
            normalized_samples[idx] = path;
            stream_has_variants[idx] = true;
        }

        for (i, arg) in self.population_variants.iter().enumerate() {
            let (idx, path) = Self::parse_variant_string(arg, i)?;
            if idx >= logical_len {
                return Err(Error::PopulationVariantIndexOutOfBounds {
                    idx,
                    max: logical_len,
                });
            }
            normalized_populations[idx] = path;
            stream_has_variants[idx] = true;
        }

        if aln_count == 1 && (!stream_has_variants[0] || !stream_has_variants[1]) {
            return Err(Error::SingleStreamMissingVariantProfiles);
        }

        self.sample_variants = normalized_samples
            .into_iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();
        self.population_variants = normalized_populations
            .into_iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();

        // 3. Output bounds validation
        if self.output.len() > logical_len {
            return Err(Error::TooManyOutputPaths {
                count: self.output.len(),
                max: logical_len,
            });
        }

        if aln_count == 1 {
            if self.ambiguous_output.len() > 1 {
                return Err(Error::SingleStreamTooManyAmbiguousOutputs {
                    count: self.ambiguous_output.len(),
                });
            }
        } else {
            if self.ambiguous_output.len() > logical_len {
                return Err(Error::TooManyAmbiguousOutputPaths {
                    count: self.ambiguous_output.len(),
                    max: logical_len,
                });
            }
        }

        if self.gap_open > 0.0 {
            self.gap_open = -self.gap_open;
        }
        if self.gap_extend > 0.0 {
            self.gap_extend = -self.gap_extend;
        }

        if self.gap_open == 0.0 || self.mismatch_penalty <= 0.0 {
            return Err(Error::InvalidPenalties);
        }

        if self.bisulfite {
            tracing::info!(
                "Bisulfite mode enabled: C→T (forward) and G→A (reverse) \
                 mismatches scored as zero-penalty conversions. \
                 Tier-2 perfect-match fast path suppressed."
            );
        }
        if self.error_model != ErrorModel::Illumina {
            tracing::info!(
                model = ?self.error_model,
                quality_calibration = self.error_model.quality_calibration(),
                gap_open = self.gap_open,
                gap_extend = self.gap_extend,
                "Non-default error model active"
            );
        }

        tracing::debug!(
            threads = self.threads,
            score_threads = self.score_threads,
            alignments = aln_count,
            gap_open = self.gap_open,
            gap_extend = self.gap_extend,
            mismatch = self.mismatch_penalty,
            "Configuration validated"
        );

        Ok(())
    }
}

// -- Per-algorithm arg structs -------------------------------------------------

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

    #[command(flatten)]
    pub(crate) parallel: crate::config::args::ParallelArgs,

    #[command(flatten)]
    pub(crate) chimeric: crate::config::args::ChimericArgs,

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
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) ambiguous_regions: Vec<PathBuf>,

    /// In-memory VCF/BCF of diagnostic positions per stream.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) diagnostic_variants: Vec<PathBuf>,

    /// BED file(s) of regions giving reads a positive score bonus.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) positive_regions: Vec<PathBuf>,

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
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) ambiguous_regions: Vec<PathBuf>,

    /// Tabix-indexed VCF/BCF of diagnostic positions per stream.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) diagnostic_variants: Vec<PathBuf>,

    /// BED.gz file(s) of positive-score regions. Tabix-indexed.
    #[arg(long, num_args = 0..=2, value_name = "FILE",
          help_heading = "Regions", value_hint = clap::ValueHint::FilePath)]
    pub(crate) positive_regions: Vec<PathBuf>,

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

#[cfg(test)]
pub(crate) mod tests;
