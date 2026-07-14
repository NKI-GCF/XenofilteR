// src/config/args.rs

use crate::penalty::ErrorModel;
use crate::region::ScoreFn;
use clap::Args;
use std::path::PathBuf;
use crate::Error;
use crate::config::StripReadSuffix;

/// Shared by every subcommand. No arity-specific fields here.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct IoArgs {
    /// Reference FASTA for CRAM decoding.
    #[arg(long, help_heading = "Input")]
    pub(crate) reference: Option<PathBuf>,

    /// Strip /1 /2 read-name suffix. auto | true | false | variable
    #[arg(short = 'R', long, default_value = "auto", help_heading = "Input")]
    pub(crate) strip_read_suffix: StripReadSuffix,

    /// Add XF:C / XR:C decision-confidence aux tags to output records.
    #[arg(short = 'A', long, default_value = "false", help_heading = "Output")]
    pub(crate) add_decision_tag: bool,

    /// Suppress @PG header line.
    #[arg(short = 'P', long, default_value = "false", help_heading = "Output")]
    pub(crate) no_program_line: bool,

    /// Write JSON summary statistics (MultiQC-compatible).
    #[arg(long, env = "XENOFILTERS_STATS_OUTPUT", help_heading = "Output")]
    pub(crate) stats_output: Option<PathBuf>,

    /// Suppress progress output to stderr.
    #[arg(long, default_value = "false", help_heading = "Misc")]
    pub(crate) quiet: bool,

    /// Verbosity: -v = INFO, -vv = DEBUG.
    #[arg(short, long, action = clap::ArgAction::Count, help_heading = "Misc")]
    pub(crate) verbose: u8,
}

/// Scoring flags, shared everywhere.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct ScoringArgs {
    #[arg(long, default_value = "illumina", help_heading = "Scoring")]
    pub(crate) error_model: ErrorModel,

    #[arg(short = 'm', long, default_value = "4.0", help_heading = "Scoring")]
    pub(crate) mismatch_penalty: f64,

    #[arg(short = 'g', long, default_value = "6.0", help_heading = "Scoring")]
    pub(crate) gap_open: f64,

    #[arg(short = 'e', long, default_value = "1.0", help_heading = "Scoring")]
    pub(crate) gap_extend: f64,

    #[arg(long, default_value_t = u32::MAX, value_name = "PHRED|auto",
          help_heading = "Scoring")]
    pub(crate) ambiguous_threshold: u32,

    #[arg(long, default_value = "0.05", help_heading = "Scoring")]
    pub(crate) warn_ambig_fraction: f64,

    #[arg(long, default_value = "false", help_heading = "Scoring")]
    pub(crate) bisulfite: bool,
}

/// Variant-rescue flags. Arity varies (see AlignmentArgs*); the flags
/// themselves are identical everywhere.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct VariantArgs {
    #[arg(short = 's', long, num_args = 0..=32, value_name = "[IDX:]FILE",
          help_heading = "Variants")]
    pub(crate) sample_variants: Vec<String>,

    #[arg(short = 'p', long, num_args = 0..=32, value_name = "[IDX:]FILE",
          help_heading = "Variants")]
    pub(crate) population_variants: Vec<String>,
}

/// Parallelism — only meaningful for namesorted (hashlookup/collated force 1).
#[derive(Args, Debug, Clone)]
pub(crate) struct ParallelArgs {
    #[arg(
        short = 't',
        long,
        default_value = "4",
        env = "XENOFILTERS_THREADS",
        help_heading = "Parallelism"
    )]
    pub(crate) threads: usize,

    #[arg(
        short = 'S',
        long,
        default_value = "1",
        env = "XENOFILTERS_SCORE_THREADS",
        help_heading = "Parallelism"
    )]
    pub(crate) score_threads: usize,
}

/// Chimeric-pair detection — only meaningful for namesorted (paired-end,
/// multi-stream). Shared verbatim between `namesorted` and `viral-integration`.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct ChimericArgs {
    /// Chimeric stream-index pairs (format A:B). Reads spanning species
    /// boundaries get XC:Z:<other_label> and are written to both outputs.
    #[arg(long, num_args = 0.., value_name = "A:B", help_heading = "Chimeric")]
    pub(crate) chimeric_pairs: Vec<String>,

    /// Labels per stream, used in XC:Z tags and stats JSON.
    #[arg(long, num_args = 0..=32, help_heading = "Chimeric")]
    pub(crate) stream_labels: Vec<String>,
}

/// In-memory region flags (hashlookup).
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct RegionArgsMemory {
    #[arg(long, num_args = 0..=2, value_name = "FILE", help_heading = "Regions")]
    pub(crate) ambiguous_regions: Vec<String>,

    #[arg(long, num_args = 0..=2, value_name = "FILE", help_heading = "Regions")]
    pub(crate) diagnostic_variants: Vec<String>,

    #[arg(long, num_args = 0..=2, value_name = "FILE", help_heading = "Regions")]
    pub(crate) positive_regions: Vec<String>,

    #[arg(long, default_value = "linear:1.0", help_heading = "Regions")]
    pub(crate) region_score_fn: ScoreFn,
}

/// Tabix-indexed region flags (collated). Same flag *names*, different
/// runtime loader (TabixBed vs AmbiguousRegions::from_bed). Kept as a
/// separate struct only because the doc strings differ (tabix requirement).
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct RegionArgsTabix {
    #[arg(long, num_args = 0..=2, value_name = "FILE.bed.gz",
          help_heading = "Regions")]
    pub(crate) ambiguous_regions: Vec<String>,

    #[arg(long, num_args = 0..=2, value_name = "FILE.vcf.gz",
          help_heading = "Regions")]
    pub(crate) diagnostic_variants: Vec<String>,

    #[arg(long, num_args = 0..=2, value_name = "FILE.bed.gz",
          help_heading = "Regions")]
    pub(crate) positive_regions: Vec<String>,

    #[arg(long, default_value = "linear:1.0", help_heading = "Regions")]
    pub(crate) region_score_fn: ScoreFn,
}

/// Output paths: arity is 1..=32 for general multi-stream use.
/// Kept separate from IoArgs because strain/viral variants want stricter caps.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct OutputArgsMulti {
    #[arg(short = 'o', long, num_args = 1..=32, help_heading = "Output")]
    pub(crate) output: Vec<String>,

    #[arg(short = 'a', long, num_args = 0..=32, help_heading = "Output")]
    pub(crate) ambiguous_output: Vec<String>,
}

/// Output paths for exactly-2-logical-stream subcommands (strain, hashlookup,
/// collated, viral-integration). Same flags, tighter arity — better --help
/// and earlier validation than a shared 1..=32 bound.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct OutputArgsPair {
    #[arg(short = 'o', long, num_args = 1..=2, help_heading = "Output")]
    pub(crate) output: Vec<String>,

    #[arg(short = 'a', long, num_args = 0..=2, help_heading = "Output")]
    pub(crate) ambiguous_output: Vec<String>,
}

// src/config/args.rs

impl IoArgs {
    pub(crate) fn validate(&self) -> Result<(), Error> {
        if let Some(ref p) = self.reference {
            if !p.exists() {
                return Err(Error::ReferenceNotFound {
                    path: p.display().to_string(),
                });
            }
        }
        Ok(())
    }
}


impl ScoringArgs {
    pub(crate) fn validate(&mut self) -> Result<(), Error> {
        if self.mismatch_penalty <= 0.0
            || self.gap_open <= 0.0
            || self.gap_extend < 0.0
        {
            return Err(Error::InvalidPenalties);
        }
        if !(0.0..=1.0).contains(&self.warn_ambig_fraction) {
            return Err(Error::InvalidWarnAmbigFraction {
                value: self.warn_ambig_fraction,
            });
        }
        Ok(())
    }

    pub(crate) fn to_penalty(&self) -> crate::penalty::Penalty {
        crate::penalty::Penalty::build(
            self.gap_open,
            self.gap_extend,
            self.mismatch_penalty,
            20,
            self.error_model,
        )
    }
}

impl VariantArgs {
    pub(crate) fn parse_indexed(
        specs: &[String],
    ) -> Result<Vec<(usize, std::path::PathBuf)>, Error> {
        specs.iter().enumerate().map(|(i, s)| {
            match s.split_once(':') {
                Some((idx_str, path))
                    if idx_str.chars().all(|c| c.is_ascii_digit()) =>
                {
                    let idx = idx_str
                        .parse::<usize>()
                        .map_err(|_| Error::InvalidVariantStreamIndex {
                            spec: s.clone(),
                        })?;
                    Ok((idx, std::path::PathBuf::from(path)))
                }
                _ => Ok((i, std::path::PathBuf::from(s))),
            }
        }).collect()
    }

    pub(crate) fn has_index(&self, idx: usize) -> Result<bool, Error> {
        let s = Self::parse_indexed(&self.sample_variants)?
            .iter().any(|(i, _)| *i == idx);
        let p = Self::parse_indexed(&self.population_variants)?
            .iter().any(|(i, _)| *i == idx);
        Ok(s || p)
    }
}

impl OutputArgsMulti {
    pub(crate) fn validate(&self, n_streams: usize) -> Result<(), Error> {
        if self.output.len() > n_streams {
            return Err(Error::TooManyOutputPaths {
                count: self.output.len(),
                max: n_streams,
            });
        }
        if !self.ambiguous_output.is_empty()
            && self.ambiguous_output.len() > n_streams
        {
            return Err(Error::TooManyAmbiguousPaths {
                given: self.ambiguous_output.len(),
                streams: n_streams,
            });
        }
        Ok(())
    }
}

impl OutputArgsPair {
    pub(crate) fn validate(&self) -> Result<(), Error> {
        if self.output.len() > 2 {
            return Err(Error::TooManyOutputPathsPair);
        }
        if self.ambiguous_output.len() > 2 {
            return Err(Error::TooManyAmbiguousPathsPair);
        }
        Ok(())
    }
}

impl ChimericArgs {
    pub(crate) fn parse_pairs(
        &self,
        n_streams: usize,
    ) -> Result<Vec<[usize; 2]>, Error> {
        parse_chimeric_pairs(&self.chimeric_pairs, n_streams)
    }
}

impl Default for ParallelArgs {
    fn default() -> Self {
        Self {
            threads: 1,
            score_threads: 1,
        }
    }
}

pub(crate) fn parse_chimeric_pairs(
    specs: &[String],
    n_streams: usize,
) -> Result<Vec<[usize; 2]>, Error> {
    let mut pairs = Vec::with_capacity(specs.len());
    for raw in specs {
        let (a_str, b_str) =
            raw.split_once(':').ok_or_else(|| Error::InvalidChimericPairFormat {
                raw: raw.clone(),
            })?;
        let a = a_str
            .trim()
            .parse::<usize>()
            .map_err(|_| Error::InvalidChimericPairFormat { raw: raw.clone() })?;
        let b = b_str
            .trim()
            .parse::<usize>()
            .map_err(|_| Error::InvalidChimericPairFormat { raw: raw.clone() })?;
        if a == b {
            return Err(Error::ChimericPairSameIndex { raw: raw.clone() });
        }
        if a >= n_streams || b >= n_streams {
            return Err(Error::ChimericPairIndexOutOfRange {
                raw: raw.clone(),
                streams: n_streams,
            });
        }
        pairs.push([a.min(b), a.max(b)]);
    }
    pairs.sort_unstable();
    pairs.dedup();
    Ok(pairs)
}

pub(crate) fn resolve_threshold(phred: u32, is_pass2: bool) -> f64 {
    let p = match phred {
        u32::MAX => if is_pass2 { 0 } else { 10 },
        p        => p,
    };
    (p as f64) * std::f64::consts::LN_10 / 10.0
}
