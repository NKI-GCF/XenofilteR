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
    #[default]
    Auto,
    True,
    False,
    Variable,
}

#[derive(Parser, Debug, Default, Clone)]
#[command(author, version, about, long_about=None)]
pub(crate) struct Config {
    /// Assign fragments matching alignment to these respective files. Writes first alignment to stdout when omitted
    #[arg(short, long, num_args = 1..ARG_MAX)]
    pub(crate) output: Vec<PathBuf>,

    /// Discard fragments distancing more in alignment to these files. Default: do not discard
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) filtered_output: Vec<PathBuf>,

    /// Write ambiguous reads (equally good mappings) to these files. Default: do not write
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) ambiguous_output: Vec<PathBuf>,

    /// ouput format of stdout
    #[arg(short = 'O', long, default_value = "sam")]
    pub(crate) stdout_format: BamFormat,

    /// Input alignments to compare. If the same readnames are consecutive and in the same order for
    /// all inputs, a low memory non-hashing strategy is adopted.
    #[arg(required = true, num_args = 1..ARG_MAX)]
    pub(crate) alignment: Vec<String>,

    /// Read first alignment from stdin; enforced with only one input alignment
    #[arg(short, long, default_value = "false")]
    pub(crate) read_from_stdin: bool,

    /// Exclude read(pair)s, unmapped in both alignments, even from the filter output.
    #[arg(short = 'U', long, default_value = "false")]
    pub(crate) discard_unmapped: bool,

    /// Mismatch penalty (affects mismatches)
    #[arg(short, long, default_value = "4", value_parser = clap::value_parser!(f64))]
    pub(crate) mismatch_penalty: f64,

    /// Gap open penalty for deletions and insertions
    #[arg(short, long, default_value = "6", value_parser = clap::value_parser!(f64))]
    pub(crate) gap_open: f64,

    /// Gap extend penalty (affects indels)
    #[arg(short = 'e', long, default_value = "1", value_parser = clap::value_parser!(f64))]
    pub(crate) gap_extend: f64,

    /// penalty for 5'- and 3'-end clipping
    #[arg(short = 'c', long, default_value = "5", value_parser = clap::value_parser!(f64))]
    pub(crate) clipping_penalty: f64,

    /// strip fastq-style /1 and /2 from read names when comparing
    #[arg(short = 'R', long, default_value = "auto")]
    pub(crate) strip_read_suffix: StripReadSuffix,

    /// Sample-specific variants used for variant-aware scoring.
    /// For single alignments, prefix with index (e.g., '0:file.vcf').
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) sample_variants: Vec<String>,

    /// Population variants used for variant-aware scoring.
    /// For single alignments, prefix with index (e.g., '1:file.vcf').
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) population_variants: Vec<String>,

    /// Add an XF tag to the records.
    #[arg(short='A', long, default_value = "false")]
    pub(crate) add_decision_tag: bool,

    /// Don't add a PG line to the output BAM header.
    #[arg(short = 'P', long, default_value = "false")]
    pub(crate) no_program_line: bool,

    /// Threshold (in phred scale) for considering two alignments equally good and thus ambiguous. Set to
    /// 0 to disable.
    #[arg(short, long, default_value = "0", value_parser = clap::value_parser!(u32).range(..=0x8000))]
    pub(crate) ambiguous_threshold: u32,

    /// Skip secondary mappings even if the primary mapping is written
    #[arg(short, long, default_value = "false")]
    pub(crate) skip_secondary: bool,

    /// Explicitly indicate that reads are paired-end
    #[arg(short, long)]
    pub(crate) is_paired: Option<bool>,

    /// Required explicit flag to allow running with only a single alignment stream
    /// using strain-specific variant profiles.
    #[arg(long)]
    pub(crate) single_alignment_mode: bool,
}

impl Config {
    pub(super) fn validate_and_init(&mut self) -> Result<()> {
        let aln_count = self.alignment.len();

        // 1. Guard against single-stream niche execution accidents first
        if aln_count == 1 {
            ensure!(
                self.single_alignment_mode,
                "Single alignment stream detected. If this is intentional for within-species disambiguation, please pass --single-alignment-mode."
            );
            ensure!(
                !self.read_from_stdin,
                "Cannot use single alignment mode with stdin because the stream must be duplicated via file system access."
            );
        } else {
            ensure!(
                !self.single_alignment_mode,
                "--single-alignment-mode can only be used with exactly 1 alignment stream."
            );
            ensure!(
                aln_count >= 2,
                "At least two alignments required when not running in single alignment mode."
            );
        }

        // Determine effective dimensions (logical comparisons)
        let logical_len = if aln_count == 1 { 2 } else { aln_count };

        // 2. Parse, validate, and normalize variant arrays to size matching comparisons
        let mut normalized_samples = vec![PathBuf::new(); logical_len];
        let mut normalized_populations = vec![PathBuf::new(); logical_len];

        // Tracker to ensure each logical stream has at least one variant track assigned to it
        let mut stream_has_variants = vec![false; logical_len];

        // Parse explicit prefixes (e.g., "0:file.vcf") or default to position zip index
        for (i, arg) in self.sample_variants.iter().enumerate() {
            let (idx, path) = Self::parse_variant_string(arg, i)?;
            ensure!(idx < logical_len, "Sample variant stream index {idx} out of bounds");
            normalized_samples[idx] = path;
            stream_has_variants[idx] = true;
        }

        for (i, arg) in self.population_variants.iter().enumerate() {
            let (idx, path) = Self::parse_variant_string(arg, i)?;
            ensure!(idx < logical_len, "Population variant stream index {idx} out of bounds");
            normalized_populations[idx] = path;
            stream_has_variants[idx] = true;
        }

        // Enforce structural strain variations requirement
        if aln_count == 1 {
            ensure!(
                stream_has_variants[0] && stream_has_variants[1],
                "Single alignment mode requires both strain slots (index 0 and 1) to have a variant profile. \
                An option like '--sample-variant 0:a.vcf --population-variant 0:b.vcf' is invalid because strain 1 has no variations."
            );
        }

        // Overwrite raw user input string vectors with fully populated, position-normalized PathBuf tracks
        // NOTE: If you change these to Vec<PathBuf> in the Config struct, update the struct definitions!
        // For compliance with remaining downstream code expecting strings, we map them back to strings.
        self.sample_variants = normalized_samples.into_iter().map(|p| p.to_string_lossy().into_owned()).collect();
        self.population_variants = normalized_populations.into_iter().map(|p| p.to_string_lossy().into_owned()).collect();

        // 3. Output bounds validation
        ensure!(
            self.output.len() <= logical_len,
            "More output paths than logical alignment processing streams specified"
        );
        ensure!(
            self.filtered_output.len() <= logical_len,
            "More filtered output paths than logical alignment processing streams specified"
        );

        if aln_count == 1 {
            ensure!(
                self.ambiguous_output.len() <= 1,
                "Only one ambiguous output file is allowed when operating on a single alignment stream."
            );
        } else {
            ensure!(
                self.ambiguous_output.len() <= logical_len,
                "More ambiguous output paths than input specified"
            );
        }
        if self.gap_open > 0.0 {
            self.gap_open = -self.gap_open;
        }
        if self.gap_extend > 0.0 {
            self.gap_extend = -self.gap_extend;
        }

        if self.gap_open == 0.0 || self.mismatch_penalty <= 0.0 {
            return Err(anyhow::anyhow!(
                "Gap open/mismatch penalties must be positive"
            ));
        }
        Ok(())
    }
    /// Internal string parser to resolve explicit vs implicit stream indices
    fn parse_variant_string(arg: &str, default_idx: usize) -> Result<(usize, PathBuf)> {
        if let Some((idx_str, path_str)) = arg.split_once(':') {
            if let Ok(idx) = idx_str.parse::<usize>() {
                return Ok((idx, PathBuf::from(path_str)));
            }
        }
        Ok((default_idx, PathBuf::from(arg)))
    }
    pub(super) fn to_penalties(&self) -> Penalty {
        let mut error_prob = [0.0_f64; MAX_Q];
        for (q, item) in error_prob.iter_mut().enumerate() {
            *item = 10f64.powf(-(q as f64) / 10.0); // x = 10^{-Q/10}
        }

        let mut log_likelihood_match = [0.0_f64; MAX_Q];
        for (q, item) in log_likelihood_match.iter_mut().enumerate() {
            *item = (1.0 - error_prob[q]).log10();
            //eprintln!("Q{} match log-likelihood: {}", q, *item);
        }

        let reference_penalty = 4.0;
        let scaling_factor = self.mismatch_penalty / reference_penalty;

        let mut log_likelihood_mismatch = [0.0f64; MAX_Q];
        for (q, item) in log_likelihood_mismatch.iter_mut().enumerate() {
            *item = -(q as f64) / 10.0 * scaling_factor;
        }

        Penalty {
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            log_likelihood_mismatch,
            log_likelihood_match,
        }
    }
}

#[cfg(test)]
pub(crate) mod tests;
