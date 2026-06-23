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
use config::{Config, MatchingAlgorithm};
use filter_algorithm::{
    collated::CollatedMatcher, hash_lookup::HashLookup, line_by_line::LineByLine,
};
use noodles::bam::record::Record as BamRecord;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::Header;
use region::{AmbiguousRegions, DiagnosticVariants};
use smallvec::{smallvec, SmallVec};
use std::collections::HashMap;
use std::path::Path;
use tracing_subscriber::{fmt, EnvFilter};

fn main() -> Result<()> {
    let mut config = Config::parse();

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

    match config.matching_algorithm {
        MatchingAlgorithm::Namesorted => run_namesorted(config),
        MatchingAlgorithm::Hashlookup => run_hashlookup(config),
        MatchingAlgorithm::Collated => run_collated(config),
    }
}

// ---------------------------------------------------------------------------
// Name-sorted — streaming merge, sequential or parallel
// ---------------------------------------------------------------------------

fn run_namesorted(mut config: Config) -> Result<()> {
    let score_threads = config.score_threads;
    let logical_loops = if config.alignment.len() == 1 {
        2
    } else {
        config.alignment.len()
    };

    if score_threads > 1 {
        let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];
        for i in 0..logical_loops {
            let idx = if config.alignment.len() == 1 { 0 } else { i };
            tracing::debug!(stream = i, path = %config.alignment[idx], "Opening stream");
            aln.push(Box::new(AlnStream::<RecordBuf>::new(&mut config, idx)?));
            ensure!(
                aln[i].next_qname() == aln[0].next_qname(),
                "Input alignments must have the same read order."
            );
        }
        tracing::info!(
            streams = logical_loops,
            score_threads,
            "Starting parallel pipeline"
        );
        LineByLine::new(config, aln)?.process_parallel()
    } else {
        let mut aln: SmallVec<[Box<dyn AlignmentStream<BamRecord>>; 2]> = smallvec![];
        for i in 0..logical_loops {
            let idx = if config.alignment.len() == 1 { 0 } else { i };
            tracing::debug!(stream = i, path = %config.alignment[idx], "Opening stream");
            aln.push(Box::new(AlnStream::<BamRecord>::new(&mut config, idx)?));
            ensure!(
                aln[i].next_qname() == aln[0].next_qname(),
                "Input alignments must have the same read order."
            );
        }
        tracing::info!(streams = logical_loops, "Starting sequential pipeline");
        LineByLine::new(config, aln)?.process_sequential()
    }
}

// ---------------------------------------------------------------------------
// Hash-lookup — position-sorted BAMs, in-memory region filtering
// ---------------------------------------------------------------------------

fn run_hashlookup(mut config: Config) -> Result<()> {
    ensure!(
        config.alignment.len() == 2,
        "--matching-algorithm hashlookup requires exactly 2 alignment streams"
    );

    let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];
    for i in 0..2 {
        tracing::debug!(stream = i, path = %config.alignment[i], "Opening stream");
        aln.push(Box::new(AlnStream::<RecordBuf>::new(&mut config, i)?));
    }

    // Build contig name → 0-based integer ID from stream-0 header.
    let name_to_id = header_name_to_id(aln[0].header());

    let bed = load_ambiguous_regions(&config.ambiguous_regions, &name_to_id)?;
    let vcf = load_diagnostic_variants(&config.diagnostic_variants, &name_to_id)?;

    tracing::info!(streams = 2, "Starting HashLookup pipeline");
    HashLookup::new(config, aln, bed, vcf)?.process()
}

// ---------------------------------------------------------------------------
// Collated — individually name-sorted streams, tabix-indexed region queries
// ---------------------------------------------------------------------------

fn run_collated(mut config: Config) -> Result<()> {
    use region::tabix_query::{TabixBed, TabixVcf};

    ensure!(
        config.alignment.len() == 2,
        "--matching-algorithm collated requires exactly 2 alignment streams"
    );

    let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];
    for i in 0..2 {
        tracing::debug!(stream = i, path = %config.alignment[i], "Opening stream");
        aln.push(Box::new(AlnStream::<RecordBuf>::new(&mut config, i)?));
    }

    let bed: [Option<TabixBed>; 2] = [
        config
            .ambiguous_regions.first()
            .filter(|s| !s.is_empty())
            .map(|s| TabixBed::open(Path::new(s)))
            .transpose()?,
        config
            .ambiguous_regions
            .get(1)
            .filter(|s| !s.is_empty())
            .map(|s| TabixBed::open(Path::new(s)))
            .transpose()?,
    ];

    let vcf: [Option<TabixVcf>; 2] = [
        config
            .diagnostic_variants.first()
            .filter(|s| !s.is_empty())
            .map(|s| TabixVcf::open(Path::new(s)))
            .transpose()?,
        config
            .diagnostic_variants
            .get(1)
            .filter(|s| !s.is_empty())
            .map(|s| TabixVcf::open(Path::new(s)))
            .transpose()?,
    ];

    tracing::info!(streams = 2, "Starting Collated pipeline");
    CollatedMatcher::new(config, aln, bed, vcf)?.process()
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// Map reference sequence name → 0-based integer ID from a SAM header.
/// The integer IDs match what noodles BAM records carry in their ref_seq_id field.
fn header_name_to_id(header: &Header) -> HashMap<String, usize> {
    header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_string(), i))
        .collect()
}

fn load_ambiguous_regions(
    specs: &[String],
    name_to_id: &HashMap<String, usize>,
) -> Result<[Option<AmbiguousRegions>; 2]> {
    Ok([
        specs.first()
            .filter(|s| !s.is_empty())
            .map(|s| AmbiguousRegions::from_bed(Path::new(s), name_to_id))
            .transpose()?,
        specs
            .get(1)
            .filter(|s| !s.is_empty())
            .map(|s| AmbiguousRegions::from_bed(Path::new(s), name_to_id))
            .transpose()?,
    ])
}

fn load_diagnostic_variants(
    specs: &[String],
    name_to_id: &HashMap<String, usize>,
) -> Result<[Option<DiagnosticVariants>; 2]> {
    Ok([
        specs.first()
            .filter(|s| !s.is_empty())
            .map(|s| DiagnosticVariants::from_vcf(Path::new(s), name_to_id))
            .transpose()?,
        specs
            .get(1)
            .filter(|s| !s.is_empty())
            .map(|s| DiagnosticVariants::from_vcf(Path::new(s), name_to_id))
            .transpose()?,
    ])
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    pub(crate) use alignment::tests::*;
    pub(crate) use aln_stream::tests::*;
}
