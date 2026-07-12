// Public API surface for benches and integration tests.
// Items that must remain crate-private to the binary stay in main.rs.

pub mod alignment;
pub mod aln_stream;
pub mod bam;
pub mod config;
pub mod error;
pub mod filter_algorithm;
pub mod penalty;
pub mod progress;
pub mod region;
pub mod stats;
pub mod variant;

use crate::region::ScoredRegions;
use aln_stream::{AlignmentStream, AlnStream};
use config::{Config, MatchingAlgorithm};
pub use error::Error;
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

// Test helpers exposed only under cfg(test) or cfg(bench).
#[cfg(any(test, feature = "bench-internals"))]
pub mod tests;

pub fn run(mut config: Config) -> Result<(), Error> {
    // existing main() body moved here verbatim
    init_tracing(&config);
    config.validate_and_init()?;

    match config.matching_algorithm {
        MatchingAlgorithm::Namesorted => run_namesorted(config),
        MatchingAlgorithm::Hashlookup => run_hashlookup(config),
        MatchingAlgorithm::Collated => run_collated(config),
    }
}

fn init_tracing(config: &Config) {
    // JSON logging when RUST_LOG_FORMAT=json (e.g. for log aggregators).
    // Plain fmt is the default.
    let use_json = std::env::var("RUST_LOG_FORMAT")
        .map(|v| v.eq_ignore_ascii_case("json"))
        .unwrap_or(false);

    let default_level = match config.verbose {
        0 => "warn",
        1 => "info",
        _ => "debug",
    };
    let filter =
        EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(default_level));

    if use_json {
        tracing_subscriber::fmt()
            .json()
            .with_env_filter(filter)
            .with_writer(std::io::stderr)
            .init();
    } else {
        tracing_subscriber::fmt()
            .with_env_filter(filter)
            .with_target(false)
            .with_writer(std::io::stderr)
            .init();
    }

    fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| EnvFilter::new(get_log_level(config.verbose))),
        )
        .with_target(false)
        .with_writer(std::io::stderr)
        .init();

    tracing::info!(version = env!("CARGO_PKG_VERSION"), "xenofilters starting");
}

fn get_log_level(verbose_count: u8) -> &'static str {
    match verbose_count {
        0 => "warn",
        1 => "info",
        _ => "debug",
    }
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
        let positive_regions = config
            .positive_regions
            .get(i)
            .filter(|s| !s.is_empty())
            .map(|s| ScoredRegions::from_bed(Path::new(s), &header_name_to_id(aln[idx].header())))
            .transpose()?
            .map(|sr| (sr, config.region_score_fn));
        tracing::debug!(stream = i, path = %config.alignment[idx], "Opening stream");
        aln.push(Box::new(AlnStream::<RecordBuf>::new(
            &mut config,
            idx,
            positive_regions,
        )?));
        if aln[i].next_qname() != aln[0].next_qname() {
            return Err(Error::InputAlignmentsMustHaveSameReadOrder);
        }
    }
    if score_threads > 1 {
        let mut lbl = LineByLine::new(&config, aln)?;
        let result = lbl.process_parallel(&config);
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
        LineByLine::new(&config, aln)?.process_sequential(&config)
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
        let positive_regions = config
            .positive_regions
            .get(i)
            .filter(|s| !s.is_empty())
            .map(|s| ScoredRegions::from_bed(Path::new(s), &header_name_to_id(aln[i].header())))
            .transpose()?
            .map(|sr| (sr, config.region_score_fn));
        aln.push(Box::new(AlnStream::<RecordBuf>::new(
            &mut config,
            i,
            positive_regions,
        )?));
    }

    // Build contig name → 0-based integer ID from stream-0 header.
    let name_to_id = header_name_to_id(aln[0].header());

    let bed = load_ambiguous_regions(&config.ambiguous_regions, &name_to_id)?;
    let vcf = load_diagnostic_variants(&config.diagnostic_variants, &name_to_id)?;

    tracing::info!(streams = 2, "Starting HashLookup pipeline");
    HashLookup::new(&config, aln, bed, vcf)?.process(&config)
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
        let positive_regions = config
            .positive_regions
            .get(i)
            .filter(|s| !s.is_empty())
            .map(|s| ScoredRegions::from_bed(Path::new(s), &header_name_to_id(aln[i].header())))
            .transpose()?
            .map(|sr| (sr, config.region_score_fn));
        aln.push(Box::new(AlnStream::<RecordBuf>::new(
            &mut config,
            i,
            positive_regions,
        )?));
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
    CollatedMatcher::new(&config, aln, bed, vcf)?.process(&config)
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
