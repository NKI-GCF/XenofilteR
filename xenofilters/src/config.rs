// src/config.rs — diff from previous version:
//   - `score_threads: usize` field added (CLI: --score-threads / -S)
//   - validate_and_init clamps score_threads to available parallelism

use crate::{
    bam::BamFormat,
    penalty::{Penalty, MAX_Q},
};
use anyhow::{ensure, Result};
use clap::{Parser, ValueEnum};
use std::path::PathBuf;

const ARG_MAX: usize = 4;

#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Default)]
pub(crate) enum StripReadSuffix {
    #[default] Auto,
    True,
    False,
    Variable,
}

#[derive(Parser, Debug, Default, Clone)]
#[command(
    author, version,
    about = "Fast alignment-based read classifier for xenograft / PDX sequencing data",
    long_about = None,
)]
pub(crate) struct Config {
    /// Input alignments to compare (name-sorted BAM, at least two required).
    #[arg(required = true, num_args = 1..ARG_MAX)]
    pub(crate) alignment: Vec<String>,

    /// Output file for the winning alignment, one per stream.
    #[arg(short, long, num_args = 1..ARG_MAX)]
    pub(crate) output: Vec<PathBuf>,

    /// Output file for all alignments (winners, filtered, and ambiguous).
    /// If set, overrides --output, --filtered-output, and --ambiguous-output.
    #[arg(short, long)]
    pub(crate) merged_output: Option<PathBuf>,

    /// Output file for the losing alignment, one per stream.
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) filtered_output: Vec<PathBuf>,

    /// Output file for ambiguous fragments.
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) ambiguous_output: Vec<PathBuf>,

    /// Output format for stdout.
    #[arg(short = 'O', long, default_value = "sam")]
    pub(crate) stdout_format: BamFormat,

    /// Number of bgzf (de)compression threads per reader/writer.
    #[arg(short = 't', long, default_value = "4")]
    pub(crate) threads: usize,

    /// Number of parallel scoring worker threads.
    ///
    /// Each worker owns its own DP scratch space and scores fragments
    /// independently.  The IO thread (reading + writing) is always single-
    /// threaded; only the log-likelihood and variant-aware scoring is
    /// parallelised.
    ///
    /// Set to 1 (the default) for deterministic output order.
    /// Set to 0 to use all available logical CPUs.
    #[arg(short = 'S', long, default_value = "1")]
    pub(crate) score_threads: usize,

    /// Increase log verbosity (-v = INFO, -vv = DEBUG). Overridden by RUST_LOG.
    #[arg(short, long, action = clap::ArgAction::Count)]
    pub(crate) verbose: u8,

    /// Read first alignment from stdin; enforced with only one input alignment.
    #[arg(short, long, default_value = "false")]
    pub(crate) read_from_stdin: bool,

    /// Exclude read-pairs that are unmapped in both alignments, even from the
    /// filtered output.
    #[arg(short = 'U', long, default_value = "false")]
    pub(crate) discard_unmapped: bool,

    /// Mismatch penalty (Phred-scaled; must be positive).
    #[arg(short, long, default_value = "4", value_parser = clap::value_parser!(f64))]
    pub(crate) mismatch_penalty: f64,

    /// Gap-open penalty for deletions and insertions (must be positive).
    #[arg(short, long, default_value = "6", value_parser = clap::value_parser!(f64))]
    pub(crate) gap_open: f64,

    /// Gap-extend penalty per base inside a gap (must be positive).
    #[arg(short = 'e', long, default_value = "1", value_parser = clap::value_parser!(f64))]
    pub(crate) gap_extend: f64,

    /// Penalty applied to 5′- and 3′-end soft-clipped bases.
    #[arg(short = 'c', long, default_value = "5", value_parser = clap::value_parser!(f64))]
    pub(crate) clipping_penalty: f64,

    /// How to handle the `/1` / `/2` FASTQ-style read-name suffixes.
    /// `auto` detects from the first record; `true` always strips; `false`
    /// never strips; `variable` handles mixed files.
    #[arg(short = 'R', long, default_value = "auto")]
    pub(crate) strip_read_suffix: StripReadSuffix,

    /// Sample-specific VCF/BCF used for variant-aware scoring.
    /// Prefix with stream index to target a specific alignment
    /// (e.g. `0:file.vcf`); otherwise assigned positionally.
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) sample_variants: Vec<String>,

    /// Population-frequency VCF/BCF used for variant-aware scoring.
    /// Same index-prefix syntax as `--sample-variants`.
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) population_variants: Vec<String>,

    /// Attach an `XF:C:<phred>` aux tag to winning records encoding the
    /// Phred-scaled confidence of the assignment (range 0–255).
    #[arg(short = 'A', long, default_value = "false")]
    pub(crate) add_decision_tag: bool,

    /// Suppress the `@PG` line that xenofilters appends to every output BAM
    /// header.
    #[arg(short = 'P', long, default_value = "false")]
    pub(crate) no_program_line: bool,

    /// Minimum Phred-score difference required to call one alignment better
    /// than another. Pairs below this threshold are written to
    /// `--ambiguous-output`. Set to 0 to disable (default).
    #[arg(short, long, default_value = "0", value_parser = clap::value_parser!(u32).range(..=0x8000))]
    pub(crate) ambiguous_threshold: u32,

    /// Skip secondary mappings even when writing the primary mapping.
    #[arg(short, long, default_value = "false")]
    pub(crate) skip_secondary: bool,

    /// Explicitly mark the dataset as paired-end (overrides auto-detection).
    #[arg(short, long)]
    pub(crate) is_paired: Option<bool>,

    /// Allow running with a single alignment stream (requires two variant
    /// profiles, one per strain slot).
    #[arg(long)]
    pub(crate) single_alignment_mode: bool,

}

impl Config {
    pub(super) fn validate_and_init(&mut self) -> Result<()> {
        let aln_count = self.alignment.len();

        if aln_count == 1 {
            ensure!(self.single_alignment_mode,
                "Single alignment stream detected. Pass --single-alignment-mode if intentional.");
            ensure!(!self.read_from_stdin,
                "Cannot use --single-alignment-mode with stdin.");
        } else {
            ensure!(!self.single_alignment_mode,
                "--single-alignment-mode can only be used with exactly 1 alignment stream.");
            ensure!(aln_count >= 2,
                "At least two alignments required outside single-alignment mode.");
        }

        if self.merged_output.is_some() {
            ensure!(
                self.output.is_empty() && self.filtered_output.is_empty() && self.ambiguous_output.is_empty(),
                "Cannot use --merged-output in combination with --output, --filtered-output, or --ambiguous-output."
            );
        }

        let logical_len = if aln_count == 1 { 2 } else { aln_count };

        let mut normalized_samples     = vec![PathBuf::new(); logical_len];
        let mut normalized_populations = vec![PathBuf::new(); logical_len];
        let mut stream_has_variants    = vec![false; logical_len];

        for (i, arg) in self.sample_variants.iter().enumerate() {
            let (idx, path) = Self::parse_variant_string(arg, i)?;
            ensure!(idx < logical_len, "Sample variant index {idx} out of bounds");
            normalized_samples[idx] = path;
            stream_has_variants[idx] = true;
        }
        for (i, arg) in self.population_variants.iter().enumerate() {
            let (idx, path) = Self::parse_variant_string(arg, i)?;
            ensure!(idx < logical_len, "Population variant index {idx} out of bounds");
            normalized_populations[idx] = path;
            stream_has_variants[idx] = true;
        }

        if aln_count == 1 {
            ensure!(stream_has_variants[0] && stream_has_variants[1],
                "Single-alignment mode requires variant profiles for both strain slots.");
        }

        self.sample_variants = normalized_samples.into_iter()
            .map(|p| p.to_string_lossy().into_owned()).collect();
        self.population_variants = normalized_populations.into_iter()
            .map(|p| p.to_string_lossy().into_owned()).collect();

        ensure!(self.output.len() <= logical_len,
            "More --output paths than alignment streams.");
        ensure!(self.filtered_output.len() <= logical_len,
            "More --filtered-output paths than alignment streams.");
        if aln_count == 1 {
            ensure!(self.ambiguous_output.len() <= 1,
                "Only one --ambiguous-output allowed in single-alignment-stream mode.");
        } else {
            ensure!(self.ambiguous_output.len() <= logical_len,
                "More --ambiguous-output paths than alignment streams.");
        }

        if self.gap_open   > 0.0 { self.gap_open   = -self.gap_open; }
        if self.gap_extend > 0.0 { self.gap_extend  = -self.gap_extend; }

        if self.gap_open == 0.0 || self.mismatch_penalty <= 0.0 {
            return Err(anyhow::anyhow!(
                "gap-open and mismatch-penalty must both be positive."
            ));
        }

        // Resolve score_threads = 0 → all available logical CPUs.
        if self.score_threads == 0 {
            self.score_threads = std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1);
            tracing::info!(score_threads = self.score_threads,
                "score_threads=0: using all available logical CPUs");
        }

        tracing::debug!(
            threads      = self.threads,
            score_threads = self.score_threads,
            alignments   = aln_count,
            gap_open     = self.gap_open,
            gap_extend   = self.gap_extend,
            mismatch     = self.mismatch_penalty,
            "Configuration validated"
        );

        Ok(())
    }

    /// Parse `"<idx>:<path>"` or fall back to `(default_idx, path)`.
    fn parse_variant_string(arg: &str, default_idx: usize) -> Result<(usize, PathBuf)> {
        if let Some((idx_str, path_str)) = arg.split_once(':')
            && let Ok(idx) = idx_str.parse::<usize>() {
                return Ok((idx, PathBuf::from(path_str)));
            }
        Ok((default_idx, PathBuf::from(arg)))
    }

    /// Build a [`Penalty`] from the current penalty parameters.
    pub(super) fn to_penalties(&self) -> Penalty {
        let mut error_prob = [0.0_f64; MAX_Q];
        for (q, item) in error_prob.iter_mut().enumerate() {
            *item = 10f64.powf(-(q as f64) / 10.0);
        }
        let mut log_likelihood_match = [0.0_f64; MAX_Q];
        for (q, item) in log_likelihood_match.iter_mut().enumerate() {
            *item = (1.0 - error_prob[q]).log10();
        }
        let scaling_factor = self.mismatch_penalty / 4.0;
        let mut log_likelihood_mismatch = [0.0f64; MAX_Q];
        for (q, item) in log_likelihood_mismatch.iter_mut().enumerate() {
            *item = -(q as f64) / 10.0 * scaling_factor;
        }
        Penalty { gap_open: self.gap_open, gap_extend: self.gap_extend,
                  log_likelihood_mismatch, log_likelihood_match }
    }
}

#[cfg(test)]
pub(crate) mod tests;
