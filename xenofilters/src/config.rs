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
    #[arg(required = true, num_args = 2..ARG_MAX)]
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

    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) sample_variants: Vec<PathBuf>,

    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) population_variants: Vec<PathBuf>,

    /// Add an XF tag to the records.
    #[arg(short, long, default_value = "false")]
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

    #[arg(short, long)]
    pub(crate) is_paired: Option<bool>,
}

impl Config {
    pub(super) fn validate_and_init(&mut self) -> Result<()> {
        ensure!(
            self.output.len() <= self.alignment.len(),
            "More output than input specified"
        );
        ensure!(
            self.filtered_output.len() <= self.alignment.len(),
            "More filtered output than input specified"
        );
        ensure!(
            self.ambiguous_output.len() <= self.alignment.len(),
            "More ambiguous output than input specified"
        );
        ensure!(
            self.alignment.len() >= 2,
            "At least two alignments required"
        );
        ensure!(
            !self.read_from_stdin || self.alignment.len() == 1,
            "Cannot read from stdin with multiple input alignments"
        );
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
