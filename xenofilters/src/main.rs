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
mod variant;

use aln_stream::{AlignmentStream, AlnStream};
use anyhow::{ensure, Result};
use clap::Parser;
use config::Config;
use filter_algorithm::line_by_line::LineByLine;
use noodles::bam::record::Record as BamRecord;
use smallvec::{smallvec, SmallVec};
use tracing_subscriber::{fmt, EnvFilter};

fn main() -> Result<()> {
    let mut config = Config::parse();

    // -- Logging setup ------------------------------------------------------
    // Priority: RUST_LOG env var > -v count.
    // -v  → INFO, -vv → DEBUG, default → WARN.
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

    // -- Open alignment streams ---------------------------------------------
    let mut aln: SmallVec<[Box<dyn AlignmentStream<BamRecord>>; 2]> = smallvec![];

    let logical_loops = if config.alignment.len() == 1 {
        2
    } else {
        config.alignment.len()
    };

    for i in 0..logical_loops {
        let target_config_idx = if config.alignment.len() == 1 { 0 } else { i };

        tracing::debug!(stream = i, path = %config.alignment[target_config_idx], "Opening alignment stream");

        aln.push(Box::new(AlnStream::new(&mut config, target_config_idx)?));

        ensure!(
            aln[i].next_qname() == aln[0].next_qname(),
            "Input alignments must have the same read order. \
             Stream 0 has '{}', stream {i} has '{}'.",
            std::str::from_utf8(aln[0].next_qname()).unwrap_or("<invalid UTF-8>"),
            std::str::from_utf8(aln[i].next_qname()).unwrap_or("<invalid UTF-8>"),
        );
    }

    tracing::info!(
        streams = logical_loops,
        "All alignment streams opened and synchronised"
    );

    LineByLine::new(config, aln)?.process()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    pub(crate) use alignment::tests::*;
    pub(crate) use aln_stream::tests::*;
}
