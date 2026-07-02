extern crate anyhow;
extern crate clap;
extern crate core;
extern crate noodles;
extern crate smallvec;

mod alignment;
mod aln_stream;
mod bam;
mod config;
mod error;
mod filter_algorithm;
mod penalty;
mod progress;
mod region;
mod stats;
mod variant;

use crate::aln_stream::{CramStream, SamStdinStream};
use crate::bam::AlnFormat;
use aln_stream::{AlignmentStream, AlnStream};
use clap::Parser;
use config::{Config, MatchingAlgorithm};
pub(crate) use error::Error;
pub(crate) use filter_algorithm::line_by_line::print_routing_counters;
use filter_algorithm::{
    collated::CollatedMatcher, hash_lookup::HashLookup, line_by_line::LineByLine,
};
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::Header;
use region::{AmbiguousRegions, DiagnosticVariants};
use smallvec::{smallvec, SmallVec};
use std::collections::HashMap;
use std::path::Path;
use tracing_subscriber::{fmt, EnvFilter};

fn get_log_level(verbose_count: u8) -> &'static str {
    match verbose_count {
        0 => "warn",
        1 => "info",
        _ => "debug",
    }
}

fn run(mut config: Config) -> Result<(), Error> {
    config.validate_and_init()?;
    match config.matching_algorithm {
        MatchingAlgorithm::Namesorted => run_namesorted(config),
        MatchingAlgorithm::Hashlookup => run_hashlookup(config),
        MatchingAlgorithm::Collated => run_collated(config),
    }
}

fn main() -> Result<(), Error> {
    let config = Config::parse();

    fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| EnvFilter::new(get_log_level(config.verbose))),
        )
        .with_target(false)
        .with_writer(std::io::stderr)
        .init();

    tracing::info!(version = env!("CARGO_PKG_VERSION"), "xenofilters starting");

    run(config)
}

// ---------------------------------------------------------------------------
// Name-sorted — streaming merge, sequential or parallel
// ---------------------------------------------------------------------------
fn run_namesorted(mut config: Config) -> Result<(), Error> {
    let stats_path = config.stats_output.clone();
    let stream_labels = config.stream_labels.clone();
    let score_threads = config.score_threads;

    // Stream construction: dispatch on input_format and path "-"
    let logical_loops = if config.alignment.len() == 1 {
        2
    } else {
        config.alignment.len()
    };

    let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];
    for i in 0..logical_loops {
        let idx = if config.alignment.len() == 1 { 0 } else { i };
        let path = config.alignment[idx].as_str();
        let stream: Box<dyn AlignmentStream<RecordBuf> + Unpin> = match config.input_format {
            AlnFormat::Cram => Box::new(CramStream::new(
                Path::new(path),
                config
                    .reference
                    .as_deref()
                    .ok_or(Error::MissingReference(path.to_string()))?,
            )?),
            AlnFormat::Sam if path == "-" => Box::new(SamStdinStream::new(&mut config, i)?),
            AlnFormat::Sam => {
                // File-based SAM: wrap in AlnStream with SAM reader.
                // TODO: add SamFileStream; for now error.
                return Err(Error::SamFileNotYetSupported);
            }
            AlnFormat::Bam => Box::new(AlnStream::<RecordBuf>::new(&mut config, idx)?),
        };
        tracing::debug!(stream = i, path = %config.alignment[idx], "Opening stream");
        aln.push(stream);
        if aln[i].next_qname() != aln[0].next_qname() {
            return Err(Error::InputAlignmentsMustHaveSameReadOrder);
        }
    }
    if score_threads > 1 {
        let mut lbl = LineByLine::new(config, aln)?;
        let result = lbl.process_parallel();
        if let Some(p) = stats_path {
            stats::write_stats(
                &p,
                &lbl.routing_counters,
                logical_loops,
                &stream_labels,
                "sample",
            )?;
        }
        result
    } else {
        tracing::info!(streams = logical_loops, "Starting sequential pipeline");
        LineByLine::new(config, aln)?.process_sequential()
    }
}

// ---------------------------------------------------------------------------
// Hash-lookup — position-sorted BAMs, in-memory region filtering
// ---------------------------------------------------------------------------

fn run_hashlookup(mut config: Config) -> Result<(), Error> {
    if config.alignment.len() != 2 {
        return Err(Error::AlgoRequiresTwoStreams);
    }

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

fn run_collated(mut config: Config) -> Result<(), Error> {
    use region::tabix_query::{TabixBed, TabixVcf};

    if config.alignment.len() != 2 {
        return Err(Error::AlgoRequiresTwoStreams);
    }

    let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];
    for i in 0..2 {
        tracing::debug!(stream = i, path = %config.alignment[i], "Opening stream");
        aln.push(Box::new(AlnStream::<RecordBuf>::new(&mut config, i)?));
    }

    let bed: [Option<TabixBed>; 2] = [
        config
            .ambiguous_regions
            .first()
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
            .diagnostic_variants
            .first()
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
) -> Result<[Option<AmbiguousRegions>; 2], Error> {
    Ok([
        specs
            .first()
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
) -> Result<[Option<DiagnosticVariants>; 2], Error> {
    Ok([
        specs
            .first()
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
    pub(crate) mod common;
    mod exploratory;
    use super::*;
    pub(crate) use alignment::tests::*;
    pub(crate) use aln_stream::tests::*;

    // Kills mutations in `header_name_to_id` (HashMap::new(), HashMap::from_iter, etc.)
    #[test]
    fn test_header_name_to_id() {
        // Construct a realistic SAM header to ensure iteration and indexing are correct
        let header_str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200";
        let header: Header = header_str.parse().expect("Failed to parse SAM header");

        let map = header_name_to_id(&header);

        assert_eq!(map.len(), 2);
        assert_eq!(map.get("chr1"), Some(&0));
        assert_eq!(map.get("chr2"), Some(&1));
    }

    #[test]
    fn test_load_ambiguous_regions_ignores_empty_strings() {
        let name_to_id = HashMap::new();

        // Two empty strings should return [None, None] without triggering a file-read error
        let specs = vec!["".to_string(), "".to_string()];
        let result = load_ambiguous_regions(&specs, &name_to_id).unwrap();
        assert!(result[0].is_none());
        assert!(result[1].is_none());

        // A single empty string should return [None, None]
        let specs_single = vec!["".to_string()];
        let result_single = load_ambiguous_regions(&specs_single, &name_to_id).unwrap();
        assert!(result_single[0].is_none());
        assert!(result_single[1].is_none());
    }

    #[test]
    fn test_load_diagnostic_variants_ignores_empty_strings() {
        let name_to_id = HashMap::new();

        let specs = vec!["".to_string(), "".to_string()];
        let result = load_diagnostic_variants(&specs, &name_to_id).unwrap();
        assert!(result[0].is_none());
        assert!(result[1].is_none());
    }

    #[test]
    fn test_namesorted_sequential_single_alignment() {
        let config = Config {
            matching_algorithm: MatchingAlgorithm::Namesorted,
            score_threads: 1, // Forces the sequential path
            alignment: vec!["tests/fixtures/dummy1.bam".into()],
            // .. populate remaining necessary fields
            ..Default::default()
        };

        // We just care that the logic branches correctly.
        let _ = run_namesorted(config);
    }

    #[test]
    fn test_namesorted_parallel_dual_alignment() {
        let config = Config {
            matching_algorithm: MatchingAlgorithm::Namesorted,
            score_threads: 2, // Forces the parallel path
            alignment: vec![
                "tests/fixtures/dummy1.bam".into(),
                "tests/fixtures/dummy2.bam".into(),
            ],
            ..Default::default()
        };

        let _ = run_namesorted(config);
    }

    #[test]
    fn test_get_log_level() {
        assert_eq!(get_log_level(0), "warn");
        assert_eq!(get_log_level(1), "info");
        assert_eq!(get_log_level(2), "debug");
        assert_eq!(get_log_level(5), "debug");
    }
}
