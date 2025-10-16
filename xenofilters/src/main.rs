extern crate anyhow;
extern crate clap;
extern crate once_cell;
extern crate rust_htslib;
extern crate core;
extern crate smallvec;

mod aln_stream;
mod filter_algorithm;
mod bam_format;
mod alignment;
mod fragment;

use std::path::PathBuf;

use anyhow::{ensure, Result};
use clap::Parser;
use once_cell::sync::Lazy;
use smallvec::{SmallVec, smallvec};

use bam_format::BamFormat;
use aln_stream::AlnStream;
//use filter_algorithm::hasher::Hasher;
use filter_algorithm::line_by_line::LineByLine;


pub use alignment::*;
pub use fragment::{FragmentState, Evaluation};

const ARG_MAX: usize = 4;
const MAX_Q: usize = 93;

static ERROR_PROB: Lazy<[f64; MAX_Q]> = Lazy::new(|| {
    let mut arr = [0.0_f64; MAX_Q];
    for (q, item) in arr.iter_mut().enumerate() {
        *item = 10f64.powf(-(q as f64) / 10.0); // x = 10^{-Q/10
    }
    arr
});

static LOG_LIKELIHOOD_MATCH: Lazy<[f64; MAX_Q]> = Lazy::new(|| {
    let mut arr = [0.0_f64; MAX_Q];
    for (q, item) in arr.iter_mut().enumerate() {
        *item = (1.0 - ERROR_PROB[q]).log10();  // ln(1 - x)
    }
    arr
});

#[derive(Parser, Debug)]
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
    #[clap(short, long, default_value = "bam")]
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
    #[clap(short, long, default_value = "1", value_parser = clap::value_parser!(f64))]
    pub gap_extend: f64,

    /// Mismatch penalty (affects mismatches)
    #[clap(short, long, default_value = "4", value_parser = clap::value_parser!(f64))]
    pub mismatch_penalty: f64,

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
}

fn main() -> Result<()> {
    let config = Config::parse();
    // hla: --favor-last-alignment --ignore-clips-beyond-contig -m 5 -u 0

    /*if config.second_unmapped_penalty == 0x8000 {
        config.second_unmapped_penalty = config.unmapped_penalty;
    }
    config.unmapped_penalty <<= 16;*/

    ensure!(config.output.len() <= config.alignment.len(), "More output than input specified");
    ensure!(config.filtered_output.len() <= config.alignment.len(), "More filtered output than input specified");
    ensure!(config.ambiguous_output.len() <= config.alignment.len(), "More ambiguous output than input specified");
    ensure!(config.alignment.len() >= 2, "At least two alignments required");
    ensure!(!config.read_from_stdin || config.alignment.len() == 1, "Cannot read from stdin with multiple input alignments");

    let mut log_likelihood_mismatch = [0.0f64; MAX_Q + 2];

    for q in 0..MAX_Q {
        // mismatch likelihood (includes penalty)
        log_likelihood_mismatch[q] = (ERROR_PROB[q] / 3.0).ln() - config.mismatch_penalty;
    }
    if config.gap_open <= 0.0 || config.gap_extend < 0.0 {
        return Err(anyhow::anyhow!("Gap open/extend penalties must be positive"));
    }
    log_likelihood_mismatch[MAX_Q] = -config.gap_open; // gap open penalty
    log_likelihood_mismatch[MAX_Q + 1] = -config.gap_extend; // gap extend penalty

    // first alignment to quick check readnames are in same name order
    let mut aln: SmallVec<[AlnStream; 2]> = smallvec![];
    for i in 0..config.alignment.len() {
        aln.push(AlnStream::new(&config, i)?);
    }
    let mut line_by_line = LineByLine::new(config, log_likelihood_mismatch, aln)?;
    line_by_line.process()
}


