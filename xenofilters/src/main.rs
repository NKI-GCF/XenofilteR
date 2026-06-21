extern crate anyhow;
extern crate clap;
extern crate core;
extern crate noodles;
extern crate smallvec;

mod alignment;
mod aln_stream;
mod bam;
mod config;
mod filter_algorithm;
mod penalty;
mod region;
mod variant;

use aln_stream::{AlignmentStream, AlnStream};
use anyhow::{ensure, Result};
use clap::Parser;
use config::Config;
use filter_algorithm::line_by_line::LineByLine;
use noodles::bam::record::Record as BamRecord;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{smallvec, SmallVec};
use tracing_subscriber::{fmt, EnvFilter};

fn main() -> Result<()> {
    let mut config = Config::parse();

    // -- Logging ---------------------------------------------------------------
    let default_level = match config.verbose {
        0 => "warn",
        1 => "info",
        _ => "debug",
    };
    fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(default_level)),
        )
        .with_target(false)
        .with_writer(std::io::stderr)
        .init();

    tracing::info!(version = env!("CARGO_PKG_VERSION"), "xenofilters starting");

    config.validate_and_init()?;

    let score_threads = config.score_threads;

    // -- Open alignment streams ------------------------------------------------
    // When score_threads > 1 we need RecordBuf (owned) throughout because
    // fragments cross thread boundaries.  When sequential we can use the lazy
    // bam::Record path.
    //
    // For simplicity we always use RecordBuf; the copy cost is negligible
    // compared to bgzf decompression.

    let logical_loops = if config.alignment.len() == 1 {
        2
    } else {
        config.alignment.len()
    };

    if score_threads > 1 {
        // Parallel path: LineByLine<RecordBuf>
        let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];

        for i in 0..logical_loops {
            let idx = if config.alignment.len() == 1 { 0 } else { i };
            tracing::debug!(stream = i, path = %config.alignment[idx], "Opening alignment stream");
            aln.push(Box::new(AlnStream::<RecordBuf>::new(&mut config, idx)?));
            ensure!(
                aln[i].next_qname() == aln[0].next_qname(),
                "Input alignments must have the same read order."
            );
        }

        tracing::info!(
            streams = logical_loops,
            score_threads,
            "Starting parallel scoring pipeline"
        );

        LineByLine::new(config, aln)?.process_parallel()
    } else {
        // Sequential path: LineByLine<BamRecord> (zero-copy lazy decoding)
        let mut aln: SmallVec<[Box<dyn AlignmentStream<BamRecord>>; 2]> = smallvec![];

        for i in 0..logical_loops {
            let idx = if config.alignment.len() == 1 { 0 } else { i };
            tracing::debug!(stream = i, path = %config.alignment[idx], "Opening alignment stream");
            aln.push(Box::new(AlnStream::<BamRecord>::new(&mut config, idx)?));
            ensure!(
                aln[i].next_qname() == aln[0].next_qname(),
                "Input alignments must have the same read order."
            );
        }

        tracing::info!(
            streams = logical_loops,
            "Starting sequential scoring pipeline"
        );

        LineByLine::new(config, aln)?.process_sequential()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    pub(crate) use alignment::tests::*;
    pub(crate) use aln_stream::tests::*;
}
