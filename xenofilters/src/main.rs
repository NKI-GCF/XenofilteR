extern crate anyhow;
extern crate clap;
extern crate core;
extern crate rust_htslib;
extern crate smallvec;

mod alignment;
mod aln_stream;
mod bam_format;
mod filter_algorithm;
mod fragment;
mod variant;
mod vcf_format;

use std::path::PathBuf;

use aln_stream::{AlignmentStream, AlnStream};
use anyhow::{Result, ensure};
use bam_format::BamFormat;
use clap::{Parser, ValueEnum};
use filter_algorithm::line_by_line::LineByLine;
use smallvec::{SmallVec, smallvec};

pub use alignment::{PreparedAlignmentPair, PreparedAlignmentPairIter};
pub use fragment::FragmentState;

const ARG_MAX: usize = 4;
const MAX_Q: usize = 93;

#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Default)]
pub enum StripReadSuffix {
    #[default]
    Auto,
    True,
    False,
    Variable,
}

#[derive(Parser, Debug, Default, Clone)]
#[clap(author, version, about, long_about=None)]
pub struct Config {
    /// Assign fragments matching alignment to these respective files. Writes first alignment to stdout when omitted
    #[clap(short, long, num_args = 0..ARG_MAX)]
    pub output: Vec<PathBuf>,

    /// Discard fragments distancing more in alignment to these files. Default: do not discard
    #[clap(short, long, num_args = 0..ARG_MAX)]
    pub filtered_output: Vec<PathBuf>,

    /// Write ambiguous reads (equally good mappings) to these files. Default: do not write
    #[clap(short, long, num_args = 0..ARG_MAX)]
    pub ambiguous_output: Vec<PathBuf>,

    /// ouput format of stdout
    #[clap(short = 'O', long, default_value = "sam")]
    pub stdout_format: BamFormat,

    /// Input alignments to compare. If the same readnames are consecutive and in the same order for
    /// all inputs, a low memory non-hashing strategy is adopted.
    #[clap(required = true, num_args = 2..ARG_MAX)]
    pub alignment: Vec<String>,

    /// Read first alignment from stdin; enforced with only one input alignment
    #[clap(short, long, default_value = "false")]
    pub read_from_stdin: bool,

    /// Exclude read(pair)s, unmapped in both alignments, even from the filter output.
    #[clap(short = 'U', long, default_value = "false")]
    pub discard_unmapped: bool,

    /// Gap open penalty (affects indels)
    #[clap(short, long, default_value = "6", value_parser = clap::value_parser!(f64))]
    pub gap_open: f64,

    /// Gap extend penalty (affects indels)
    #[clap(short = 'e', long, default_value = "1", value_parser = clap::value_parser!(f64))]
    pub gap_extend: f64,

    /// Mismatch penalty (affects mismatches)
    #[clap(short, long, default_value = "4", value_parser = clap::value_parser!(f64))]
    pub mismatch_penalty: f64,

    /// strip fastq-style /1 and /2 from read names when comparing
    #[clap(short = 'R', long, default_value = "Auto")]
    pub strip_read_suffix: StripReadSuffix,

    #[clap(short, long, num_args = 0..ARG_MAX)]
    pub sample_variants: Vec<PathBuf>,

    #[clap(short, long, num_args = 0..ARG_MAX)]
    pub population_variants: Vec<PathBuf>,

    /*/// Number of mismatches allowed in the second alignment
    #[clap(short, long, default_value = "4")]
    pub mismatch_threshold: u32,

    /// Penalty given to unmapped reads in favor of the alternative alignment. Set to 0 to disable
    #[clap(short, long, default_value = "8", value_parser = clap::value_parser!(u32).range(..0x8000))]
    pub unmapped_penalty: u32,

    /// Same, for mate. Defaults to same as first mapping penalty
    #[clap(long, default_value = "32768", value_parser = clap::value_parser!(u32).range(..=0x8000))]
    pub second_unmapped_penalty: u32,*/
    /// Skip secondary mappings even if the primary mapping is written
    #[clap(short, long, default_value = "false")]
    pub skip_secondary: bool,

    #[clap(short, long)]
    pub is_paired: Option<bool>,
}

struct Penalties {
    pub gap_open: f64,
    pub gap_extend: f64,
    pub log_likelihood_mismatch: [f64; MAX_Q],
    pub log_likelihood_match: [f64; MAX_Q],
}

impl Config {
    pub fn nreads(&self) -> usize {
        match self.is_paired {
            Some(true) => 2,
            Some(false) => 1,
            None => panic!("is_paired not set"), // don't call before AlnStream::new
        }
    }
    fn validate_and_init(&mut self) -> Result<()> {
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
    fn to_penalties(&self) -> Penalties {
        let mut error_prob = [0.0_f64; MAX_Q];
        for (q, item) in error_prob.iter_mut().enumerate() {
            *item = 10f64.powf(-(q as f64) / 10.0); // x = 10^{-Q/10}
        }

        let mut log_likelihood_match = [0.0_f64; MAX_Q];
        for (q, item) in log_likelihood_match.iter_mut().enumerate() {
            *item = (1.0 - error_prob[q]).log10();
        }

        let reference_penalty = 4.0;
        let scaling_factor = self.mismatch_penalty / reference_penalty;

        let mut log_likelihood_mismatch = [0.0f64; MAX_Q];
        for (q, item) in log_likelihood_mismatch.iter_mut().enumerate() {
            *item = -(q as f64) / 10.0 * scaling_factor;
        }

        Penalties {
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            log_likelihood_mismatch,
            log_likelihood_match,
        }
    }
}

fn main() -> Result<()> {
    let mut config = Config::parse();
    config.validate_and_init()?;

    // first alignment to quick check readnames are in same name order
    let mut aln: SmallVec<[Box<dyn AlignmentStream>; 2]> = smallvec![];
    for i in 0..config.alignment.len() {
        aln.push(Box::new(AlnStream::new(&mut config, i)?));
        ensure!(
            aln[i].next_qname() == aln[0].next_qname(),
            "Input alignments must have the same read order."
        );
    }

    LineByLine::new(config, aln).process()
}

#[cfg(test)]
pub mod tests {
    use super::*;
    pub use alignment::tests::*;
    pub use aln_stream::tests::*;
    pub use bam_format::BamFormat;
    impl Config {
        pub fn default_no_strip() -> Self {
            let mut c = Config::default();
            c.strip_read_suffix = StripReadSuffix::False;
            c
        }
    }
}
