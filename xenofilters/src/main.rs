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
use anyhow::{Result, ensure};
use clap::Parser;
use filter_algorithm::line_by_line::LineByLine;
use smallvec::{SmallVec, smallvec};
use config::Config;
use noodles::bam::record::Record as BamRecord;


fn main() -> Result<()> {
    let mut config = Config::parse();
    config.validate_and_init()?;

    let mut aln: SmallVec<[Box<dyn AlignmentStream<BamRecord>>; 2]> = smallvec![];

    // Determine the actual loop boundary
    let logical_loops = if config.alignment.len() == 1 { 2 } else { config.alignment.len() };

    // first alignment to quick check readnames are in same name order
    for i in 0..logical_loops {
        // If single_alignment_mode is active, both iteration 0 and 1 instantiate
        // a reader pointing to config.alignment[0]
        let target_config_idx = if config.alignment.len() == 1 { 0 } else { i };

        aln.push(Box::new(AlnStream::new(&mut config, target_config_idx)?));
        ensure!(
            aln[i].next_qname() == aln[0].next_qname(),
            "Input alignments must have the same read order."
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
