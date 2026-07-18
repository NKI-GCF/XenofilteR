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

use crate::file_spec::{path_for_stream, FileSpec};
use crate::filter_algorithm::{strain::StrainArgs, viral_integration::ViralIntegrationArgs};
use crate::region::tabix_load::{open_tabix_bed, open_tabix_scored, open_tabix_vcf};
use clap::CommandFactory;
use clap_complete::generate;
use config::{
    args::parse_chimeric_pairs, AlgorithmCommand, Cli, CollatedArgs, HashlookupArgs, NamesortedArgs,
};
pub use error::Error;
use filter_algorithm::{
    collated::CollatedMatcher,
    hash_lookup::HashLookup,
    line_by_line::{LineByLine, MAX_STREAMS},
};
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::Header;
use region::{
    load::{
        load_ambiguous_regions_memory, load_diagnostic_variants_memory,
        load_positive_regions_memory,
    },
    AmbiguousRegions, DiagnosticVariants,
};
use std::collections::HashMap;
use std::path::Path;
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

    // Validate common args, detect pass 2, resolve thresholds.
    common.validate_and_init()?;

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
    let stats_path = args.common.output.stats_output.clone();
    let stream_labels = args.chimeric.stream_labels.clone();
    let score_threads = args.parallel.score_threads;
    let mut runconfig = args.to_runconfig();
    let aln = runconfig.open_streams_unified(args.parallel.threads)?;
    let n = aln.len();

    if score_threads > 1 {
        let mut lbl = LineByLine::new(&args, aln)?;
        let result = lbl.process_parallel(&args);
        if let Some(p) = stats_path {
            crate::stats::write_stats(&p, &lbl.routing_counters, n, &stream_labels, "sample")?;
        }
        result
    } else {
        LineByLine::new(&args, aln)?.process_sequential(&args)
    }
}

// ---------------------------------------------------------------------------
// Hash-lookup — position-sorted BAMs, in-memory region filtering
// ---------------------------------------------------------------------------

fn run_hashlookup(mut args: HashlookupArgs) -> Result<(), Error> {
    let bed = load_ambiguous_regions(&args.ambiguous_regions, &name_to_id)?;
    let vcf = load_diagnostic_variants(&args.diagnostic_variants, &name_to_id)?;
    let mut runconfig = args.to_runconfig()?;
    let aln = runconfig.open_streams_raw_bam()?;
    let name_to_id = aln
        .iter()
        .map(|s| crate::variant::name_to_id::header_name_to_id(s.header()))
        .collect::<Vec<HashMap<_, _>>>();

    HashLookup::<RecordBuf>::new_from_hashlookup(&args, aln, bed, vcf, [None, None])?.process(&args)
}

// ---------------------------------------------------------------------------
// Collated — individually name-sorted streams, tabix-indexed region queries
// ---------------------------------------------------------------------------

fn run_collated(mut args: CollatedArgs) -> Result<(), Error> {
    use crate::region::tabix_query::{TabixBed, TabixVcf};
    use std::path::Path;

    let mut runconfig = args.to_runconfig()?;

    let mut aln: smallvec::SmallVec<[Box<dyn crate::aln_stream::AlignmentStream<RecordBuf>>; 2]> =
        smallvec::smallvec![];
    for i in 0..2 {
        aln.push(Box::new(crate::aln_stream::AlnStream::<RecordBuf>::new(
            &mut run_cfg,
            i,
            None,
        )?));
    }

    let bed: [Option<TabixBed>; 2] = [
        path_for_stream(&args.ambiguous_regions, 0)
            .map(|s| TabixBed::open(Path::new(s)))
            .transpose()?,
        path_for_stream(&args.ambiguous_regions, 1)
            .map(|s| TabixBed::open(Path::new(s)))
            .transpose()?,
    ];
    let vcf: [Option<TabixVcf>; 2] = [
        path_for_stream(&args.diagnostic_variants, 0)
            .map(|s| TabixVcf::open(Path::new(s)))
            .transpose()?,
        path_for_stream(&args.diagnostic_variants, 1)
            .map(|s| TabixVcf::open(Path::new(s)))
            .transpose()?,
    ];

    crate::filter_algorithm::collated::CollatedMatcher::<RecordBuf>::new_from_collated(
        &args, aln, bed, vcf,
    )?
    .process(&args.common)
}

// ---------------------------------------------------------------------------
// Strain — specialized collated matcher for strain detection
// ---------------------------------------------------------------------------

fn run_strain(mut args: StrainArgs) -> Result<(), Error> {
    let mut run = args.into_run_config()?;
    let threads = run.threads;
    let aln = open_streams_unified(&mut run, threads)?;
    let n = aln.len();
    let score_threads = run.score_threads;
    let stats_path = run.output.stats_output.clone();
    let stream_labels = run.stream_labels.clone();
    if score_threads > 1 {
        let mut lbl = LineByLine::new(&run, aln)?;
        let result = lbl.process_parallel(&run);
        if let Some(p) = stats_path {
            crate::stats::write_stats(&p, &lbl.routing_counters, n, &stream_labels, "sample")?;
        }
        result
    } else {
        LineByLine::new(&run, aln)?.process_sequential(&run)
    }
}

// ---------------------------------------------------------------------------
// Viral integration — specialized collated matcher for virus-host chimeras
// ---------------------------------------------------------------------------

fn run_viral_integration(mut args: ViralIntegrationArgs) -> Result<(), Error> {
    let mut run = args.into_run_config()?;
    let threads = run.threads;
    let aln = open_streams_unified(&mut run, threads)?;
    let n = aln.len();
    let score_threads = run.score_threads;
    let stats_path = run.output.stats_output.clone();
    let stream_labels = run.stream_labels.clone();
    if score_threads > 1 {
        let mut lbl = LineByLine::new(&run, aln)?;
        let result = lbl.process_parallel(&run);
        if let Some(p) = stats_path {
            crate::stats::write_stats(&p, &lbl.routing_counters, n, &stream_labels, "sample")?;
        }
        result
    } else {
        LineByLine::new(&run, aln)?.process_sequential(&run)
    }
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

fn load_specs<T, F>(
    specs: &[FileSpec],
    name_to_id: &[HashMap<String, usize>],
    mut load_fn: F,
) -> Result<[Option<T>; 2], Error>
where
    F: FnMut(&Path, &HashMap<String, usize>) -> Result<T, Error>,
{
    Ok([
        path_for_stream(specs, 0)
            .map(|s| load_fn(Path::new(s), &name_to_id[0]))
            .transpose()?,
        path_for_stream(specs, 1)
            .map(|s| load_fn(Path::new(s), &name_to_id[1]))
            .transpose()?,
    ])
}

fn load_ambiguous_regions(
    specs: &[FileSpec],
    name_to_id: &[HashMap<String, usize>],
) -> Result<[Option<AmbiguousRegions>; 2], Error> {
    load_specs(specs, name_to_id, |p, name_to_id| {
        AmbiguousRegions::from_bed(p, name_to_id)
    })
}

fn load_diagnostic_variants(
    specs: &[FileSpec],
    name_to_id: &Vec<HashMap<String, usize>>,
) -> Result<[Option<DiagnosticVariants>; 2], Error> {
    load_specs(specs, name_to_id, |p, name_to_id| {
        DiagnosticVariants::from_vcf(p, name_to_id)
    })
}
