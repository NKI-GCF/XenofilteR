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

use std::path::PathBuf;
use std::sync::LazyLock;

use aln_stream::AlnStream;
use anyhow::{Result, ensure};
use bam_format::BamFormat;
use clap::Parser;
use filter_algorithm::line_by_line::LineByLine;
use smallvec::{SmallVec, smallvec};

pub use alignment::*;
pub use fragment::{Evaluation, FragmentState};

const ARG_MAX: usize = 4;
const MAX_Q: usize = 93;
const REFERENCE_PENALTY: f64 = 4.0;

static ERROR_PROB: LazyLock<[f64; MAX_Q]> = LazyLock::new(|| {
    let mut arr = [0.0_f64; MAX_Q];
    for (q, item) in arr.iter_mut().enumerate() {
        *item = 10f64.powf(-(q as f64) / 10.0); // x = 10^{-Q/10}
    }
    arr
});

static LOG_LIKELIHOOD_MATCH: LazyLock<[f64; MAX_Q]> = LazyLock::new(|| {
    let mut arr = [0.0_f64; MAX_Q];
    for (q, item) in arr.iter_mut().enumerate() {
        *item = (1.0 - ERROR_PROB[q]).log10();
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
    #[clap(short = 'O', long, default_value = "bam")]
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
    #[clap(short = 'R', long)]
    pub strip_read_suffix: Option<bool>,

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
    let mut config = Config::parse();

    ensure!(
        config.output.len() <= config.alignment.len(),
        "More output than input specified"
    );
    ensure!(
        config.filtered_output.len() <= config.alignment.len(),
        "More filtered output than input specified"
    );
    ensure!(
        config.ambiguous_output.len() <= config.alignment.len(),
        "More ambiguous output than input specified"
    );
    ensure!(
        config.alignment.len() >= 2,
        "At least two alignments required"
    );
    ensure!(
        !config.read_from_stdin || config.alignment.len() == 1,
        "Cannot read from stdin with multiple input alignments"
    );

    if config.gap_open <= 0.0 || config.gap_extend < 0.0 || config.mismatch_penalty <= 0.0 {
        return Err(anyhow::anyhow!(
            "Gap open/extend penalties must be positive"
        ));
    }

    // first alignment to quick check readnames are in same name order
    let mut aln: SmallVec<[AlnStream; 2]> = smallvec![];
    for i in 0..config.alignment.len() {
        aln.push(AlnStream::new(&mut config, i)?);
        ensure!(
            aln[i].next_qname() == aln[0].next_qname(),
            "Input alignments must have the same read order."
        );
    }

    LineByLine::new(config, aln).process()
}
