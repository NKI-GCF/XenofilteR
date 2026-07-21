// Public API surface for benches and integration tests.
// Items that must remain crate-private to the binary stay in main.rs.

pub mod alignment;
pub mod aln_stream;
pub mod bam;
pub mod config;
pub mod error;
pub mod file_spec;
pub mod filter_algorithm;
pub mod penalty;
pub mod progress;
pub mod region;
pub mod reporting;
pub mod stats;
pub mod variant;

use clap::CommandFactory;
use clap_complete::generate;
use config::{
    run_config::RunConfig, AlgorithmCommand, Cli, CollatedArgs, HashlookupArgs, MatchingAlgorithm,
    NamesortedArgs,
};
pub use error::Error;
use file_spec::path_for_stream;
use filter_algorithm::{
    hash_lookup::HashLookup, line_by_line::LineByLine, strain::StrainArgs,
    viral_integration::ViralIntegrationArgs,
};
use noodles::sam::alignment::record_buf::RecordBuf;
use region::tabix_query::{TabixBed, TabixVcf};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use tracing_subscriber::{fmt, EnvFilter};

// Test helpers exposed only under cfg(test) or cfg(bench).
#[cfg(any(test, feature = "bench-internals"))]
pub mod tests;

pub fn run(mut cli: Cli) -> Result<(), Error> {
    // Completion must run before logging init (stdout must be clean).
    if let AlgorithmCommand::Completion(ref comp) = cli.command {
        let shell = comp.shell;
        let mut cmd = Cli::command();
        let bin_name = cmd.get_name().to_string();
        match shell.to_clap_shell() {
            Some(sh) => generate(sh, &mut cmd, bin_name, &mut std::io::stdout()),
            None => {
                // Nu shell via clap_complete_nushell
                #[cfg(feature = "nu-completion")]
                {
                    use clap_complete_nushell::Nushell;
                    generate(Nushell, &mut cmd, bin_name, &mut std::io::stdout());
                }
                #[cfg(not(feature = "nu-completion"))]
                {
                    eprintln!(
                        "Nu shell completion not compiled in. \
                               Rebuild with --features nu-completion."
                    );
                    std::process::exit(1);
                }
            }
        }
        return Ok(());
    }
    let common = cli.command.common_mut();
    init_tracing(common.io.verbose);

    match cli.command {
        AlgorithmCommand::Namesorted(args) => run_namesorted(args),
        AlgorithmCommand::Hashlookup(args) => run_hashlookup(args),
        AlgorithmCommand::Collated(args) => run_collated(args),
        AlgorithmCommand::Strain(args) => run_strain(args),
        AlgorithmCommand::ViralIntegration(args) => run_viral_integration(args),
        AlgorithmCommand::Completion(_) => unreachable!(),
    }
}

fn init_tracing(verbose: u8) {
    // JSON logging when RUST_LOG_FORMAT=json (e.g. for log aggregators).
    // Plain fmt is the default.
    let use_json = std::env::var("RUST_LOG_FORMAT")
        .map(|v| v.eq_ignore_ascii_case("json"))
        .unwrap_or(false);

    let default_level = match verbose {
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
                .unwrap_or_else(|_| EnvFilter::new(get_log_level(verbose))),
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
fn run_namesorted(args: NamesortedArgs) -> Result<(), Error> {
    let stats_path = args.common.io.stats_output.clone();
    let chimeric = args.chimeric.clone();
    let score_threads = args.parallel.score_threads;
    let run = args.into_run_config()?;
    let stream_labels = chimeric.stream_labels.clone();
    let chimeric_pairs = chimeric.parse_pairs(stream_labels.len())?;
    run_line_by_line(
        stats_path,
        stream_labels,
        chimeric_pairs,
        score_threads,
        run,
    )
}

fn run_line_by_line(
    stats_path: Option<PathBuf>,
    stream_labels: Vec<String>,
    chimeric_pairs: Vec<[usize; 2]>,
    mut score_threads: usize,
    mut run: RunConfig,
) -> Result<(), Error> {
    let aln = run.open_streams_unified(MatchingAlgorithm::Namesorted, run.threads)?;

    // Clamp score_threads to aln.len() * 2 as a sanity bound; beyond that
    // the IO thread becomes the bottleneck and extra workers are idle
    let n = aln.len();
    score_threads = score_threads.max(n * 2);

    let mut lbl = LineByLine::new(&run, aln, stream_labels.clone(), chimeric_pairs)?;
    if score_threads > 1 {
        tracing::info!(
            score_threads,
            "Fragment scoring will use a parallel worker pool. \
             Output order is nondeterministic."
        );
        let result = lbl.process_parallel(&run, score_threads);
        if let Some(p) = stats_path {
            crate::stats::write_stats(&p, &lbl.routing_counters, n, &stream_labels, "sample")?;
        }
        result
    } else {
        lbl.process_sequential(&run)
    }
}

// ---------------------------------------------------------------------------
// Hash-lookup — position-sorted BAMs, in-memory region filtering
// ---------------------------------------------------------------------------

fn run_hashlookup(args: HashlookupArgs) -> Result<(), Error> {
    let bed: [Option<TabixBed>; 2] = [
        path_for_stream(&args.ambiguous_regions, 0)
            .map(|s| TabixBed::open(Path::new(s)))
            .transpose()?,
        path_for_stream(&args.ambiguous_regions, 1)
            .map(|s| TabixBed::open(Path::new(s)))
            .transpose()?,
    ];
    let vcf: [Option<TabixVcf>; 2] = [
        path_for_stream(&args.distinct_variants, 0)
            .map(|s| TabixVcf::open(Path::new(s)))
            .transpose()?,
        path_for_stream(&args.distinct_variants, 1)
            .map(|s| TabixVcf::open(Path::new(s)))
            .transpose()?,
    ];
    let mut run = args.into_run_config()?;
    let aln = run.open_streams_unified(MatchingAlgorithm::Hashlookup, run.threads)?;
    let _name_to_id = aln
        .iter()
        .map(|s| crate::variant::name_to_id::header_name_to_id(s.header()))
        .collect::<Vec<HashMap<_, _>>>();

    HashLookup::<RecordBuf>::new(&run, aln, bed, vcf, [None, None])?.process(&run)
}

// ---------------------------------------------------------------------------
// Collated — individually name-sorted streams
// ---------------------------------------------------------------------------

fn run_collated(args: CollatedArgs) -> Result<(), Error> {
    let bed: [Option<TabixBed>; 2] = [
        path_for_stream(&args.ambiguous_regions, 0)
            .map(|s| TabixBed::open(Path::new(s)))
            .transpose()?,
        path_for_stream(&args.ambiguous_regions, 1)
            .map(|s| TabixBed::open(Path::new(s)))
            .transpose()?,
    ];
    let vcf: [Option<TabixVcf>; 2] = [
        path_for_stream(&args.distinct_variants, 0)
            .map(|s| TabixVcf::open(Path::new(s)))
            .transpose()?,
        path_for_stream(&args.distinct_variants, 1)
            .map(|s| TabixVcf::open(Path::new(s)))
            .transpose()?,
    ];
    let mut run = args.into_run_config()?;
    let aln = run.open_streams_unified(MatchingAlgorithm::Collated, run.threads)?;

    crate::filter_algorithm::collated::CollatedMatcher::<RecordBuf>::new(&run, aln, bed, vcf)?
        .process(&run)
}

// ---------------------------------------------------------------------------
// Strain — specialized collated matcher for strain detection
// ---------------------------------------------------------------------------

fn run_strain(args: StrainArgs) -> Result<(), Error> {
    let stats_path = args.io.stats_output.clone();
    let score_threads = args.parallel.score_threads;
    let stream_labels = vec![];
    let chimeric_pairs = vec![];
    let run = args.into_run_config()?;
    run_line_by_line(
        stats_path,
        stream_labels,
        chimeric_pairs,
        score_threads,
        run,
    )
}

// ---------------------------------------------------------------------------
// Viral integration — specialized collated matcher for virus-host chimeras
// ---------------------------------------------------------------------------

fn run_viral_integration(args: ViralIntegrationArgs) -> Result<(), Error> {
    let stats_path = args.io.stats_output.clone();
    let score_threads = args.score_threads;
    let chimeric_pairs = vec![[0, 1]];
    let stream_labels = args.stream_labels.clone();
    let run = args.into_run_config()?;
    run_line_by_line(
        stats_path,
        stream_labels,
        chimeric_pairs,
        score_threads,
        run,
    )
}
