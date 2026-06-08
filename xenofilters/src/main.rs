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

    // first alignment to quick check readnames are in same name order
    let mut aln: SmallVec<[Box<dyn AlignmentStream<BamRecord>>; 2]> = smallvec![];
    for i in 0..config.alignment.len() {
        aln.push(Box::new(AlnStream::new(&mut config, i)?));
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
    pub(crate) use bam_format::BamFormat;
}
