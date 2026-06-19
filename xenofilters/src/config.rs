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
#[command(
    author,
    version,
    about = "Fast alignment-based read classifier for xenograft / PDX sequencing data",
    long_about = None,
    after_help = "\
EXAMPLES:
  # Compare graft vs host, write winners to graft.bam, discards to host.bam:
  xenofilters graft_aln.bam host_aln.bam -o graft.bam -o host.bam

  # With filtered and ambiguous outputs, 8 decompression threads:
  xenofilters graft.bam host.bam -o best_g.bam -o best_h.bam \\
              -f filt_g.bam -f filt_h.bam -a ambig.bam -t 8

  # Variant-aware rescue using a sample VCF for graft (index 0):
  xenofilters graft.bam host.bam -o best.bam -s 0:graft_variants.vcf
"
)]
pub(crate) struct Config {
    /// Input alignments to compare. Reads must arrive in identical name order
    /// across all inputs. Low-memory streaming mode is used automatically.
    #[arg(required = true, num_args = 1..ARG_MAX)]
    pub(crate) alignment: Vec<String>,

    /// Assign fragments matching alignment to these respective files.
    /// Writes the first alignment to stdout when omitted.
    #[arg(short, long, num_args = 1..ARG_MAX)]
    pub(crate) output: Vec<PathBuf>,

    /// Discard fragments that align worse in a given stream to these files.
    /// Default: discard silently.
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) filtered_output: Vec<PathBuf>,

    /// Write ambiguous reads (equally good mappings) to these files.
    /// Default: do not write.
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) ambiguous_output: Vec<PathBuf>,

    /// Output format for stdout (sam | bam | cram).
    #[arg(short = 'O', long, default_value = "sam")]
    pub(crate) stdout_format: BamFormat,


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
            ensure!(
                self.single_alignment_mode,
                "Single alignment stream detected. If this is intentional for \
                 within-species disambiguation, pass --single-alignment-mode."
            );
            ensure!(
                !self.read_from_stdin,
                "Cannot use --single-alignment-mode with stdin: the stream \
                 must be duplicated via file-system access."
            );
        } else {
            ensure!(
                !self.single_alignment_mode,
                "--single-alignment-mode can only be used with exactly 1 \
                 alignment stream."
            );
            ensure!(
                aln_count >= 2,
                "At least two alignments are required when not running in \
                 single-alignment mode."
            );
        }

        let logical_len = if aln_count == 1 { 2 } else { aln_count };

        let mut normalized_samples     = vec![PathBuf::new(); logical_len];
        let mut normalized_populations = vec![PathBuf::new(); logical_len];
        let mut stream_has_variants    = vec![false; logical_len];

        for (i, arg) in self.sample_variants.iter().enumerate() {
            let (idx, path) = Self::parse_variant_string(arg, i)?;
            ensure!(idx < logical_len,
                "Sample-variant stream index {idx} is out of bounds \
                 (only {logical_len} streams configured).");
            normalized_samples[idx] = path;
            stream_has_variants[idx] = true;
        }

        for (i, arg) in self.population_variants.iter().enumerate() {
            let (idx, path) = Self::parse_variant_string(arg, i)?;
            ensure!(idx < logical_len,
                "Population-variant stream index {idx} is out of bounds \
                 (only {logical_len} streams configured).");
            normalized_populations[idx] = path;
            stream_has_variants[idx] = true;
        }

        if aln_count == 1 {
            ensure!(
                stream_has_variants[0] && stream_has_variants[1],
                "Single-alignment mode requires both strain slots (index 0 and 1) \
                 to have a variant profile. An option like \
                 `--sample-variants 0:a.vcf --population-variants 0:b.vcf` is \
                 invalid because strain 1 has no variants."
            );
        }

        self.sample_variants = normalized_samples
            .into_iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();
        self.population_variants = normalized_populations
            .into_iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();

        ensure!(self.output.len() <= logical_len,
            "More --output paths ({}) than alignment streams ({logical_len}).",
            self.output.len());
        ensure!(self.filtered_output.len() <= logical_len,
            "More --filtered-output paths ({}) than alignment streams ({logical_len}).",
            self.filtered_output.len());

        if aln_count == 1 {
            ensure!(self.ambiguous_output.len() <= 1,
                "Only one --ambiguous-output file is allowed in \
                 single-alignment-stream mode.");
        } else {
            ensure!(self.ambiguous_output.len() <= logical_len,
                "More --ambiguous-output paths ({}) than alignment streams ({logical_len}).",
                self.ambiguous_output.len());
        }

        // clap accepts positive values for ergonomics; internally we need negatives.
        if self.gap_open > 0.0 {
            self.gap_open = -self.gap_open;
        }
        if self.gap_extend > 0.0 {
            self.gap_extend = -self.gap_extend;
        }

        if self.gap_open == 0.0 || self.mismatch_penalty <= 0.0 {
            return Err(anyhow::anyhow!(
                "gap-open ({}) and mismatch-penalty ({}) must both be positive.",
                self.gap_open.abs(), self.mismatch_penalty
            ));
        }

        Ok(())
    }

    /// Parse `"<idx>:<path>"` or fall back to `(default_idx, path)`.
    fn parse_variant_string(arg: &str, default_idx: usize) -> Result<(usize, PathBuf)> {
        if let Some((idx_str, path_str)) = arg.split_once(':') {
            if let Ok(idx) = idx_str.parse::<usize>() {
                return Ok((idx, PathBuf::from(path_str)));
            }
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

        // Scale mismatch by the ratio to the reference penalty (4.0) so that
        // the user-visible `--mismatch-penalty` maps intuitively onto the
        // log-likelihood space.
        let scaling_factor = self.mismatch_penalty / 4.0;
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
