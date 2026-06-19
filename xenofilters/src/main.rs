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

fn main() -> Result<()> {
    let mut config = Config::parse();
    config.validate_and_init()?;

    let mut aln: SmallVec<[Box<dyn AlignmentStream<BamRecord>>; 2]> = smallvec![];

    let logical_loops = if config.alignment.len() == 1 {
        2
    } else {
        config.alignment.len()
    };

    for i in 0..logical_loops {
        let target_config_idx = if config.alignment.len() == 1 { 0 } else { i };

        aln.push(Box::new(AlnStream::new(&mut config, target_config_idx)?));
        ensure!(
            aln[i].next_qname() == aln[0].next_qname(),
            "Input alignments must have the same read order. \
             Stream 0 has '{}', stream {i} has '{}'.",
            std::str::from_utf8(aln[0].next_qname()).unwrap_or("<invalid UTF-8>"),
            std::str::from_utf8(aln[i].next_qname()).unwrap_or("<invalid UTF-8>"),
        );
    }

    LineByLine::new(config, aln)?.process()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    pub(crate) use alignment::tests::*;
    pub(crate) use aln_stream::tests::*;
}
