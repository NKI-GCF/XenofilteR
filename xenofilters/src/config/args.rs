// src/config/args.rs

use crate::bam::AlnFormat;
use crate::file_spec::{FileSpec, path_for_stream};
use crate::filter_algorithm::line_by_line::MAX_STREAMS;
use crate::penalty::ErrorModel;
use crate::penalty::Penalty;
use crate::region::ScoreFn;
use crate::{Error, ensure};
use clap::Args;
use std::ops::RangeInclusive;
use std::path::PathBuf;

/// Shared by every subcommand. No arity-specific fields here.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct IoArgs {
    /// Input alignments to compare. Only one required for subcommand 'strain'.
    #[arg(required = true, num_args = 1..MAX_STREAMS, help_heading = "Input")]
    pub alignment: Vec<PathBuf>,

    /// Reference FASTA for CRAM decoding.
    #[arg(long, value_name = "[IDX:]FILE", help_heading = "Input")]
    pub(crate) reference: Vec<FileSpec>,

    #[arg(short = 'o', long, num_args = 1..=MAX_STREAMS, help_heading = "Output")]
    pub(crate) output: Vec<PathBuf>,

    #[arg(short = 'a', long, num_args = 0..=MAX_STREAMS, help_heading = "Output")]
    pub(crate) ambiguous_output: Vec<PathBuf>,

    /// Write JSON summary statistics (MultiQC-compatible).
    #[arg(long, env = "XENOFILTERS_STATS_OUTPUT", help_heading = "Output")]
    pub(crate) stats_output: Option<PathBuf>,

    /// Output format: sam | bam | cram. Default: bam.
    #[arg(short = 'O', long, default_value = "bam", help_heading = "Output")]
    pub(crate) stdout_format: AlnFormat,

    /// Add XF:C / XR:C decision-confidence aux tags to output records.
    #[arg(short = 'A', long, default_value = "false", help_heading = "Output")]
    pub(crate) add_decision_tag: bool,

    /// Suppress @PG header line.
    #[arg(short = 'P', long, default_value = "false", help_heading = "Output")]
    pub(crate) no_program_line: bool,

    /// Input format: sam | bam | cram. Default: bam.
    #[arg(long, default_value = "bam", help_heading = "Input")]
    pub(crate) input_format: AlnFormat,

    /// Keep discarded reads alongside ambiguous.
    #[arg(long, default_value_t = false, hide = true)]
    pub(crate) write_discarded: bool,

    /// Discard unmapped reads (default true). Hidden because it's a
    #[arg(short = 'U', long, default_value_t = true, hide = true)]
    pub(crate) discard_unmapped: bool,

    /// Skip secondary alignments (default false).
    #[arg(short = 'k', long, default_value_t = false, hide = true)]
    pub(crate) skip_secondary: bool,

    /// Suppress progress output to stderr.
    #[arg(long, default_value = "false")]
    pub(crate) quiet: bool,

    /// Verbosity: -v = INFO, -vv = DEBUG.
    #[arg(short, long, action = clap::ArgAction::Count)]
    pub(crate) verbose: u8,

    // Mutable state set during stream init -- not CLI args.
    #[arg(skip)]
    pub(crate) is_pass2: bool,
}

/// Scoring flags, shared everywhere.
#[derive(Args, Debug, Clone, Default)]
pub struct ScoringArgs {
    /// Error model for scoring. Default: illumina.
    #[arg(long, default_value = "illumina", help_heading = "Scoring")]
    pub error_model: ErrorModel,

    /// Mismatch penalty (PHRED). Default: from error model: 4.0 for Illumina.
    #[arg(short = 'm', long, default_value_t = f64::INFINITY, help_heading = "Scoring")]
    pub mismatch_penalty: f64,

    /// Gap open penalty (PHRED). Default: from error model: 6.0 for Illumina.
    #[arg(short = 'g', long, default_value_t = f64::INFINITY, help_heading = "Scoring")]
    pub gap_open: f64,

    /// Gap extend penalty (PHRED). Default: from error model: 1.0 for Illumina.
    #[arg(short = 'e', long, default_value_t = f64::INFINITY, help_heading = "Scoring")]
    pub gap_extend: f64,

    /// A supplemnetary read counts as a gap with bases. Default: 20.
    #[arg(short = 'J', long, default_value_t = 20, help_heading = "Scoring")]
    pub chimeric_junction_bases: u32,

    /// Threshold for ambiguous reads (PHRED). Default: auto (10 for pass1, 0 for pass2).
    #[arg(long, default_value_t = u32::MAX, value_name = "PHRED|auto",
          help_heading = "Scoring", default_value = "auto")]
    pub ambiguous_threshold: u32,

    /// Warn if ambiguous fraction exceeds this value (0.0-1.0). Default: 0.05.
    #[arg(long, default_value = "0.05", help_heading = "Scoring")]
    pub warn_ambig_fraction: f64,

    /// Bisulfite scoring mode. Default: false.
    #[arg(long, default_value = "false", help_heading = "Scoring")]
    pub bisulfite: bool,

    // TODO: if variants are below this, skip them at load time.
    /// Minimum p_variant required for a variant to contribute rescue
    /// evidence, in addition to the structural requirement p_variant > 0.5.
    /// Default 0.5 (no additional gate beyond the structural minimum).
    /// Raise this (e.g. 0.9) to require materially stronger evidence before
    /// a variant is allowed to influence a fragment's score -- analogous to
    /// a significance cutoff on top of the sign requirement.
    #[arg(long, default_value_t = 0.5, help_heading = "Variants")]
    pub min_rescue_p_variant: f64,
}

/// Variant-rescue flags. Arity varies (see AlignmentArgs*); the flags
/// themselves are identical everywhere.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct RelatedArgs {
    #[arg(short = 's', long, num_args = 0..=MAX_STREAMS, value_name = "[IDX:]FILE",
          help_heading = "Variants")]
    pub(crate) sample_variants: Vec<FileSpec>,

    #[arg(short = 'p', long, num_args = 0..=MAX_STREAMS, value_name = "[IDX:]FILE",
          help_heading = "Variants")]
    pub(crate) population_variants: Vec<FileSpec>,

    /// BED file(s) of regions giving reads a positive score bonus.
    #[arg(long, num_args = 0..=MAX_STREAMS, value_name = "[IDX:]FILE.bed.gz",
          help_heading = "Regions")]
    pub(crate) positive_regions: Vec<FileSpec>,

    /// Score function for --positive-regions BED score column.
    /// Format: fn[:weight]  fn  in  {linear, log, constant, overlap_fraction}
    #[arg(long, default_value = "linear:1.0", help_heading = "Regions")]
    pub(crate) region_score_fn: ScoreFn,

    // FIXME: only expand if necessary.
    #[arg(
        long,
        default_value_t = false,
        requires = "reference",
        help_heading = "Variants"
    )]
    pub(crate) expand_indels: bool,

    #[arg(
        long,
        default_value_t = 50,
        requires = "expand-indels",
        help_heading = "Variants"
    )]
    pub(crate) indel_expand_padding: usize,

    /// Minimum population allele frequency required for a variant to be
    /// loaded from --population-variants. Default 0.9. Any variant below
    /// this AF cannot structurally produce a positive rescue delta anyway
    /// (p_variant = AF must exceed 0.5), so a low default wastes memory and
    /// overlap-query time loading variants that can never fire. 0.9 aligns
    /// the default with this tool's actual supported use case -- near-fixed,
    /// diagnostic AF differences between species/strains -- rather than
    /// general population-frequency-weighted rescue, which the underlying
    /// scoring model does not support for AF < 0.5 regardless of this flag.
    #[arg(long, default_value_t = 0.9, help_heading = "Variants")]
    pub(crate) min_population_af: f64,

    /// Minimum GQ required for a variant to be loaded from --sample-variants.
    /// Default 20.
    #[arg(long, default_value_t = 20.0, help_heading = "Variants")]
    pub(crate) min_sample_gq: f64,
}

/// Parallelism -- only meaningful for namesorted (hashlookup/collated force 1).
#[derive(Args, Debug, Clone)]
pub(crate) struct ParallelArgs {
    /// Number of threads to use for alignment scoring. Default: 4.
    #[arg(
        short = 't',
        long,
        default_value = "4",
        env = "XENOFILTERS_THREADS",
        help_heading = "Parallelism"
    )]
    pub(crate) threads: usize,

    /// Number of threads to use for scoring. Default: 1.
    #[arg(
        short = 'S',
        long,
        default_value = "1",
        env = "XENOFILTERS_SCORE_THREADS",
        help_heading = "Parallelism"
    )]
    pub(crate) score_threads: usize,
}

/// Chimeric-pair detection -- only meaningful for namesorted (paired-end,
/// multi-stream). Shared verbatim between `namesorted` and `viral-integration`.
#[derive(Args, Debug, Clone)]
pub(crate) struct ChimericArgs {
    /// Chimeric stream-index pairs (format A:B). Reads spanning species
    /// boundaries get XC:Z:<other_label> and are written to both outputs.
    #[arg(long, num_args = 0.., value_name = "A:B", help_heading = "Chimeric")]
    pub(crate) chimeric_pairs: Vec<String>,

    /// Labels per stream, used in XC:Z tags and stats JSON.
    #[arg(long, num_args = 0..=MAX_STREAMS, help_heading = "Chimeric")]
    pub(crate) stream_labels: Vec<String>,

    /// Minimum mapped bases required on EACH arm of a candidate read-split,
    /// expressed as a FRACTION of read length (not an absolute base count).
    /// Default 0.1 (10%). The previous implementation used a fixed 15bp
    /// floor, which is trivially satisfied by long reads (ONT/HiFi) even in
    /// low-complexity regions with no real breakpoint signal; expressing
    /// this as a fraction keeps the same relative stringency across read
    /// length regimes. Needs empirical validation against a simulated
    /// integration-breakpoint dataset before trusting defaults in a
    /// publication.
    #[arg(long, default_value_t = 0.10, help_heading = "Chimeric")]
    pub(crate) chimeric_min_mapped_frac: f64,

    /// Maximum allowed overlap between the two arms' mapped read ranges, as
    /// a fraction of read length. Default 0.20 (20%).
    #[arg(long, default_value_t = 0.20, help_heading = "Chimeric")]
    pub(crate) chimeric_max_overlap_frac: f64,

    /// Minimum combined (union) coverage of the read by both arms, as a
    /// fraction of read length. Default 0.80 (80%).
    #[arg(long, default_value_t = 0.80, help_heading = "Chimeric")]
    pub(crate) chimeric_min_union_frac: f64,
}

impl Default for ChimericArgs {
    fn default() -> Self {
        Self {
            chimeric_pairs: vec![],
            stream_labels: vec![],
            chimeric_min_mapped_frac: 0.10,
            chimeric_max_overlap_frac: 0.20,
            chimeric_min_union_frac: 0.80,
        }
    }
}

/// Tabix-indexed region flags. Hashmap and collated algorithm only.
#[derive(Args, Debug, Clone, Default)]
pub(crate) struct SegregateArgs {
    /// BED file(s) of regions to segregate reads into a separate output stream.
    #[arg(long, num_args = 0..=2, value_name = "[IDX:]FILE.bed.gz",
          help_heading = "Regions")]
    pub(crate) ambiguous_regions: Vec<FileSpec>,

    /// VCF file(s) of distinct variants to segregate reads into a separate output stream.
    #[arg(long, num_args = 0..=MAX_STREAMS, value_name = "[IDX:]FILE.vcf.gz",
          help_heading = "Regions")]
    pub(crate) distinct_variants: Vec<FileSpec>,
}

impl IoArgs {
    pub(crate) fn validate(
        &self,
        max_stdin: usize,
        streams: RangeInclusive<usize>,
    ) -> Result<(), Error> {
        // stdin: at most one stream, namesorted only.
        let stdin_count = self
            .alignment
            .iter()
            .filter(|p| p.to_string_lossy() == "-")
            .count();
        ensure!(stdin_count <= max_stdin, Error::TooManyStdinStreams);

        // CRAM sanity: any .cram input requires --reference.
        let cram_lacks_ref =
            self.alignment.iter().enumerate().any(|(i, p)| {
                p.ends_with(".cram") && path_for_stream(&self.reference, i).is_none()
            });
        ensure!(!cram_lacks_ref, Error::CramRequiresReference);

        let n = self.alignment.len();
        let max = *streams.end();
        ensure!(
            streams.contains(&n),
            Error::InvalidStreamCount {
                n,
                min: *streams.start(),
                max
            }
        );
        let n_streams = if max == 1 { 2 } else { n };

        ensure!(
            self.output.len() <= n_streams,
            Error::TooManyOutputPaths {
                count: self.output.len(),
                max: n_streams,
            }
        );
        ensure!(
            self.ambiguous_output.is_empty() || self.ambiguous_output.len() <= n_streams,
            Error::TooManyAmbiguousPaths {
                given: self.ambiguous_output.len(),
                streams: n_streams,
            }
        );
        Ok(())
    }
}

impl ScoringArgs {
    pub fn validate(&mut self) -> Result<(), Error> {
        if self.mismatch_penalty <= 0.0 || self.gap_open <= 0.0 || self.gap_extend < 0.0 {
            return Err(Error::InvalidPenalties);
        }
        if !(0.0..=1.0).contains(&self.warn_ambig_fraction) {
            return Err(Error::InvalidWarnAmbigFraction {
                value: self.warn_ambig_fraction,
            });
        }
        // An affine gap model requires opening a gap to cost at least as much as extending one,
        // or the optimizer prefers many independent small gaps over one contiguous gap.
        if self.gap_open < self.gap_extend {
            return Err(Error::GapOpenBelowGapExtend {
                gap_open: self.gap_open,
                gap_extend: self.gap_extend,
            });
        }
        if !(0.0..=1.0).contains(&self.min_rescue_p_variant) {
            return Err(Error::InvalidMinRescuePVariant {
                value: self.min_rescue_p_variant,
            });
        }
        Ok(())
    }

    pub fn to_penalty(&self) -> Penalty {
        let gap_open = match self.gap_open {
            f64::INFINITY => self.error_model.default_gap_open(),
            _ => self.gap_open,
        };
        let gap_extend = match self.gap_extend {
            f64::INFINITY => self.error_model.default_gap_extend(),
            _ => self.gap_extend,
        };
        let mismatch_penalty = match self.mismatch_penalty {
            f64::INFINITY => self.error_model.default_mismatch_penalty(),
            _ => self.mismatch_penalty,
        };

        Penalty::build(
            gap_open,
            gap_extend,
            mismatch_penalty,
            self.chimeric_junction_bases,
            self.error_model,
            self.min_rescue_p_variant,
        )
    }
}

impl RelatedArgs {
    pub(crate) fn validate(&self) -> Result<(), Error> {
        if !(0.0..=1.0).contains(&self.min_population_af) {
            return Err(Error::InvalidMinPopulationAf {
                value: self.min_population_af,
            });
        }
        if self.min_sample_gq < 0.0 {
            return Err(Error::InvalidMinSampleGq {
                value: self.min_sample_gq,
            });
        }
        Ok(())
    }
    pub(crate) fn has_index(&self, idx: usize) -> bool {
        path_for_stream(&self.sample_variants, idx).is_some()
            || path_for_stream(&self.population_variants, idx).is_some()
    }
}

impl ChimericArgs {
    pub(crate) fn parse_pairs(&self, n_streams: usize) -> Result<Vec<[usize; 2]>, Error> {
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
        let (a_str, b_str) = raw
            .split_once(':')
            .ok_or_else(|| Error::InvalidChimericPairFormat { raw: raw.clone() })?;
        let a = a_str
            .trim()
            .parse::<usize>()
            .map_err(|_| Error::InvalidChimericPairFormat { raw: raw.clone() })?;
        let b = b_str
            .trim()
            .parse::<usize>()
            .map_err(|_| Error::InvalidChimericPairFormat { raw: raw.clone() })?;
        if a == b {
            return Err(Error::ChimericPairsIdenticalIndices { raw: raw.clone() });
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
        u32::MAX => {
            if is_pass2 {
                0
            } else {
                10
            }
        }
        p => p,
    };
    (p as f64) * std::f64::consts::LN_10 / 10.0
}

#[cfg(test)]
mod chimeric_pair_tests {
    use super::*;

    #[test]
    fn valid_pair_normalizes_order() {
        assert_eq!(
            parse_chimeric_pairs(&["1:0".into()], 2).unwrap(),
            vec![[0, 1]]
        );
    }

    #[test]
    fn missing_colon_is_format_error() {
        assert!(matches!(
            parse_chimeric_pairs(&["01".into()], 2),
            Err(Error::InvalidChimericPairFormat { .. })
        ));
    }

    #[test]
    fn non_numeric_index_is_format_error() {
        assert!(matches!(
            parse_chimeric_pairs(&["a:b".into()], 2),
            Err(Error::InvalidChimericPairFormat { .. })
        ));
    }

    #[test]
    fn identical_indices_rejected() {
        assert!(matches!(
            parse_chimeric_pairs(&["1:1".into()], 2),
            Err(Error::ChimericPairsIdenticalIndices { .. })
        ));
    }

    #[test]
    fn out_of_range_index_rejected() {
        assert!(matches!(
            parse_chimeric_pairs(&["0:5".into()], 2),
            Err(Error::ChimericPairIndexOutOfRange { .. })
        ));
    }

    #[test]
    fn duplicate_pairs_after_normalization_are_deduped() {
        // "0:1" and "1:0" normalize to the same pair.
        let got = parse_chimeric_pairs(&["0:1".into(), "1:0".into()], 2).unwrap();
        assert_eq!(got, vec![[0, 1]]);
    }

    #[test]
    fn multiple_pairs_are_sorted() {
        let got = parse_chimeric_pairs(&["2:3".into(), "0:1".into()], 4).unwrap();
        assert_eq!(got, vec![[0, 1], [2, 3]]);
    }

    #[test]
    fn whitespace_around_indices_is_trimmed() {
        assert_eq!(
            parse_chimeric_pairs(&[" 0 : 1 ".into()], 2).unwrap(),
            vec![[0, 1]]
        );
    }

    #[test]
    fn empty_specs_yields_empty_vec() {
        assert!(parse_chimeric_pairs(&[], 2).unwrap().is_empty());
    }

    #[test]
    fn min_population_af_out_of_range_rejected() {
        let args = RelatedArgs {
            min_population_af: 1.5,
            ..Default::default()
        };
        assert!(matches!(
            args.validate(),
            Err(Error::InvalidMinPopulationAf { .. })
        ));
    }

    #[test]
    fn min_rescue_p_variant_out_of_range_rejected() {
        let mut s = ScoringArgs {
            gap_open: 6.0,
            gap_extend: 1.0,
            mismatch_penalty: 4.0,
            min_rescue_p_variant: -0.1,
            ..Default::default()
        };
        assert!(matches!(
            s.validate(),
            Err(Error::InvalidMinRescuePVariant { .. })
        ));
    }

    #[test]
    fn min_sample_gq_negative_is_rejected() {
        let args = RelatedArgs {
            min_sample_gq: -5.0,
            ..Default::default()
        };
        assert!(matches!(
            args.validate(),
            Err(Error::InvalidMinSampleGq { .. })
        ));
    }
}
