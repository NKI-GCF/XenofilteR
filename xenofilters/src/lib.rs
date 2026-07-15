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
pub mod reporting;
pub mod stats;
pub mod variant;

use aln_stream::{open::open_streams_unified, AlignmentStream, AlnStream};
use config::{args::parse_chimeric_pairs, AlgorithmCommand, Cli, CollatedArgs, HashlookupArgs, NamesortedArgs};
use clap::CommandFactory;
use clap_complete::generate;
pub use error::Error;
use filter_algorithm::{
    collated::CollatedMatcher,
    hash_lookup::HashLookup,
    line_by_line::{LineByLine, MAX_STREAMS},
};
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::Header;
use region::{AmbiguousRegions, DiagnosticVariants, ScoredRegions, load::{


}};
use smallvec::{smallvec, SmallVec};
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
    init_tracing(common.verbose);

    // Validate common args, detect pass 2, resolve thresholds.
    common.validate_and_init()?;

    match cli.command {
        AlgorithmCommand::Namesorted(args) => run_namesorted(args),
        AlgorithmCommand::Hashlookup(args) => run_hashlookup(args),
        AlgorithmCommand::Collated(args) => run_collated(args),
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
fn run_namesorted(mut args: NamesortedArgs) -> Result<(), Error> {
    let n = args.common.alignment.len();
    if n == 0 {
        return Err(Error::InsufficientAlignmentStreams { count: n });
    }
    if n > MAX_STREAMS {
        return Err(Error::MaxStreamsExceeded {
            count: n,
            max: MAX_STREAMS,
        });
    }

    // Parse chimeric pairs
    args.common.parsed_chimeric_pairs = parse_chimeric_pairs(&args.chimeric_pairs, n)?;

    let aln = open_streams_unified(&mut args.common, args.threads)?;
    let counters = if args.score_threads > 1 {
        LineByLine::new_from_namesorted(&args, aln)?.process_parallel()?
    } else {
        LineByLine::new_from_namesorted(&args, aln)?.process_sequential()?
    };
    print_routing_counters(&counters, "namesorted", &args.common);
    write_stats_if_configured(&counters, &args.common, n, "namesorted")
}

// ---------------------------------------------------------------------------
// Hash-lookup — position-sorted BAMs, in-memory region filtering
// ---------------------------------------------------------------------------

fn run_hashlookup(mut args: HashlookupArgs) -> Result<(), Error> {
    let n = args.common.alignment.len();
    if n < 2 {
        return Err(Error::InsufficientAlignmentStreams { count: n });
    }

    let header_name_to_id = open_first_header_name_to_id(&args.common)?;
    let bed = load_ambiguous_regions_memory(&args.ambiguous_regions, &header_name_to_id)?;
    let vcf = load_diagnostic_variants_memory(&args.diagnostic_variants, &header_name_to_id)?;
    let pos = load_positive_regions_memory(
        &args.positive_regions,
        &header_name_to_id,
        args.region_score_fn,
    )?;
    let aln = open_streams_raw_bam(&mut args.common)?;
    let counters = HashLookup::new_from_hashlookup(&args, aln, bed, vcf, pos)?.process()?;
    print_routing_counters(&counters, "hashlookup", &args.common);
    write_stats_if_configured(&counters, &args.common, n, "hashlookup")
}

// ---------------------------------------------------------------------------
// Collated — individually name-sorted streams, tabix-indexed region queries
// ---------------------------------------------------------------------------

fn run_collated(mut args: CollatedArgs) -> Result<(), Error> {
    let n = args.common.alignment.len();
    if n < 2 {
        return Err(Error::InsufficientAlignmentStreams { count: n });
    }

    let bed = open_tabix_bed(&args.ambiguous_regions)?;
    let vcf = open_tabix_vcf(&args.diagnostic_variants)?;
    let pos = open_tabix_scored(&args.positive_regions, args.region_score_fn)?;
    let aln = open_streams_unified(&mut args.common, 1)?;
    let counters = CollatedMatcher::new_from_collated(&args, aln, bed, vcf, pos)?.process()?;
    print_routing_counters(&counters, "collated", &args.common);
    write_stats_if_configured(&counters, &args.common, n, "collated")
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

fn load_specs<T, F>(specs: &[String], mut load_fn: F) -> Result<[Option<T>; 2], Error>
where
    F: FnMut(&Path) -> Result<T, Error>,
{
    Ok([
        specs
            .first()
            .filter(|s| !s.is_empty())
            .map(|s| load_fn(Path::new(s)))
            .transpose()?,
        specs
            .get(1)
            .filter(|s| !s.is_empty())
            .map(|s| load_fn(Path::new(s)))
            .transpose()?,
    ])
}

fn load_ambiguous_regions(
    specs: &[String],
    name_to_id: &HashMap<String, usize>,
) -> Result<[Option<AmbiguousRegions>; 2], Error> {
    load_specs(specs, |p| AmbiguousRegions::from_bed(p, name_to_id))
}

fn load_diagnostic_variants(
    specs: &[String],
    name_to_id: &HashMap<String, usize>,
) -> Result<[Option<DiagnosticVariants>; 2], Error> {
    load_specs(specs, |p| DiagnosticVariants::from_vcf(p, name_to_id))
}
