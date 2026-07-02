use crate::Error;
use crate::filter_algorithm::line_by_line::MAX_STREAMS;
use crate::{
    bam::AlnFormat,
    penalty::{MAX_Q, Penalty},
};
use clap::{Parser, ValueEnum};
use std::path::PathBuf;

const ARG_MAX: usize = 4;

#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Default)]
pub(crate) enum StripReadSuffix {
    #[default]
    Auto,
    True,
    False,
    Variable,
}

/// Fragment-matching algorithm.
#[derive(Copy, Clone, Debug, ValueEnum, PartialEq, Default)]
pub(crate) enum MatchingAlgorithm {
    /// Streaming merge. Both BAMs must be in identical query-name order.
    /// Lowest memory. Fastest. Supports multithreading. Default.
    #[default]
    Namesorted,

    /// Hash-table matching. Each BAM must be collated (all records for a read
    /// name contiguous) but the two BAMs may present fragments in different
    /// orders. Memory proportional to name-order skew. Single-threaded only.
    /// Output order is not guaranteed.
    Collated,

    /// Hash-table matching. Works on arbitrary (non-sorted, non-collated) BAM
    /// input. High memory usage — proportional to the number of in-flight
    /// fragments. Single-threaded only. Preserves driving-stream (stream 0)
    /// output order.
    Hashlookup,
}

#[derive(Parser, Debug, Default, Clone)]
#[command(
    author, version,
    about = "Fast alignment-based read classifier for xenograft / PDX sequencing data",
    long_about = None,
)]
pub(crate) struct Config {
    /// Input alignments to compare. If the same readnames are consecutive and in the same order for
    /// all inputs, a low memory non-hashing strategy is adopted.
    #[arg(required = true, num_args = 1..MAX_STREAMS)]
    pub(crate) alignment: Vec<String>,

    // FIXME: make this merged_output, and not discarded_output, ambiguous_output, and output.
    // Then we can have a single output file for all alignments, and the user can choose to discard
    // or keep the discarded and ambiguous reads based on read groups.
    /// Assign fragments matching alignment to these respective files. Writes first alignment to stdout when omitted
    #[arg(short, long, num_args = 1..MAX_STREAMS)]
    pub(crate) output: Vec<PathBuf>,

    /// Output file for all alignments (winners, discarded, and ambiguous).
    /// If set, overrides --output, --discarded-output, and --ambiguous-output.
    #[arg(short, long)]
    pub(crate) merged_output: Option<PathBuf>,

    /// Discard fragments distancing more in alignment to these files. Default: do not discard
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) discarded_output: Vec<PathBuf>,

    /// Write ambiguous reads (equally good mappings) to these files. Default: do not write
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) ambiguous_output: Vec<PathBuf>,

    /// Input format. Must match actual file content.
    #[arg(long, default_value = "bam")]
    pub(crate) input_format: AlnFormat,

    /// Output format of stdout
    #[arg(short = 'O', long, default_value = "sam")]
    pub(crate) stdout_format: AlnFormat,

    /// Reference FASTA for CRAM decoding (required when --input-format cram).
    #[arg(long)]
    pub(crate) reference: Option<PathBuf>,

    /// Suppress per-fragment progress output to stderr.
    #[arg(long)]
    pub(crate) quiet: bool,

    /// Write JSON summary statistics to this path (MultiQC-compatible).
    #[arg(long)]
    pub(crate) stats_output: Option<PathBuf>,

    /// Number of bgzf (de)compression threads per reader/writer.
    #[arg(short = 't', long, default_value = "4")]
    pub(crate) threads: usize,

    /// Number of parallel scoring worker threads.
    ///
    /// Each worker owns its own DP scratch space and scores fragments
    /// independently.  The IO thread (reading + writing) is always single-
    /// threaded; only the log-likelihood and variant-aware scoring is
    /// parallelised.
    ///
    /// Set to 1 (the default) for deterministic output order.
    /// Set to 0 to use all available logical CPUs.
    #[arg(short = 'S', long, default_value = "1")]
    pub(crate) score_threads: usize,

    /// Increase log verbosity (-v = INFO, -vv = DEBUG). Overridden by RUST_LOG.
    #[arg(short, long, action = clap::ArgAction::Count)]
    pub(crate) verbose: u8,

    /// Read first alignment from stdin; enforced with only one input alignment
    #[arg(short, long, default_value = "false")]
    pub(crate) read_from_stdin: bool,

    /// Exclude read(pair)s, unmapped in both alignments, even from the filter output.
    #[arg(short = 'U', long, default_value = "false")]
    pub(crate) discard_unmapped: bool,

    /// Mismatch penalty (affects mismatches)
    #[arg(short, long, default_value = "4", value_parser = clap::value_parser!(f64))]
    pub(crate) mismatch_penalty: f64,

    /// Gap open penalty for deletions and insertions
    #[arg(short, long, default_value = "6", value_parser = clap::value_parser!(f64))]
    pub(crate) gap_open: f64,

    /// Gap extend penalty (affects indels)
    #[arg(short = 'e', long, default_value = "1", value_parser = clap::value_parser!(f64))]
    pub(crate) gap_extend: f64,

    /// Penalty for 5'- and 3'-end clipping
    #[arg(short = 'c', long, default_value = "5", value_parser = clap::value_parser!(f64))]
    pub(crate) clipping_penalty: f64,

    /// Strip fastq-style /1 and /2 from read names when comparing
    #[arg(short = 'R', long, default_value = "auto")]
    pub(crate) strip_read_suffix: StripReadSuffix,

    /// Sample-specific variants used for variant-aware scoring.
    /// For single alignments, prefix with index (e.g., '0:file.vcf').
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) sample_variants: Vec<String>,

    /// Population variants used for variant-aware scoring.
    /// For single alignments, prefix with index (e.g., '1:file.vcf').
    #[arg(short, long, num_args = 0..ARG_MAX)]
    pub(crate) population_variants: Vec<String>,

    /// Add an XF tag to the records.
    #[arg(short = 'A', long, default_value = "false")]
    pub(crate) add_decision_tag: bool,

    /// Don't add a PG line to the output BAM header.
    #[arg(short = 'P', long, default_value = "false")]
    pub(crate) no_program_line: bool,

    /// Threshold (in phred scale) for considering two alignments equally good and thus ambiguous.
    /// Set to 0 to disable.
    #[arg(short, long, default_value = "0", value_parser = clap::value_parser!(u32).range(..=0x8000))]
    pub(crate) ambiguous_threshold: u32,

    /// Skip secondary mappings even if the primary mapping is written
    #[arg(short, long, default_value = "false")]
    pub(crate) skip_secondary: bool,

    /// Explicitly indicate that reads are paired-end
    #[arg(short, long)]
    pub(crate) is_paired: Option<bool>,

    /// Required explicit flag to allow running with only a single alignment stream
    /// using strain-specific variant profiles.
    #[arg(long)]
    pub(crate) single_alignment_mode: bool,

    /// Fragment-matching algorithm.
    ///
    /// namesorted (default): streaming merge; requires identical query-name order across
    /// all input BAMs; lowest memory; fastest; supports -t/--threads.
    ///
    /// collated: each BAM must be collated (all records for a read name contiguous) but
    /// the two BAMs may present fragments in different orders; memory proportional to
    /// name-order skew; single-threaded only; output order not guaranteed.
    ///
    /// hashlookup: works on arbitrary (non-sorted, non-collated) BAM input; highest
    /// memory usage; single-threaded only; preserves driving-stream (stream 0) output order.
    #[arg(long, default_value = "namesorted")]
    pub(crate) matching_algorithm: MatchingAlgorithm,

    /// BED file of ambiguous genomic regions per stream (positional: stream 0, then 1).
    /// Reads overlapping these regions are forced through full log-likelihood scoring.
    /// Collated: must be bgzf-compressed and tabix-indexed (.bed.gz + .tbi).
    /// HashLookup: loaded fully into memory.
    #[arg(long, num_args = 0..=2)]
    pub(crate) ambiguous_regions: Vec<String>,

    /// VCF/BCF of species-diagnostic positions per stream (positional: stream 0, then 1).
    /// Reads overlapping these positions are forced through full scoring.
    /// Same indexing and compression rules as --ambiguous-regions.
    #[arg(long, num_args = 0..=2)]
    pub(crate) diagnostic_variants: Vec<String>,

    /// Base-length constant used in the supplementary-alignment chimeric-junction
    /// penalty:  penalty = gap_open + chimeric_junction_bases × gap_extend.
    ///
    /// This value replaces the per-record non-clipped base count used in earlier
    /// versions.  The default of 20 is chosen to make one supplementary alignment
    /// of typical chimeric span costlier than a 20-base insertion but cheaper
    /// than mapping to the wrong species entirely.
    #[arg(short = 'J', long, default_value = "20",
          value_parser = clap::value_parser!(u32).range(0..=10000))]
    pub(crate) chimeric_junction_bases: u32,

    /// Pairs of stream indices that may produce chimeric (cross-species) fragments.
    ///
    /// Format: "A:B" where A and B are 0-based stream indices.
    /// When a paired-end fragment has mates that split across a configured pair
    /// (some mates mapped only in stream A, complementary mates only in stream B),
    /// both streams' records are written to their assigned outputs with an `XC:Z:`
    /// SAM tag identifying the other stream — no stream is discarded.
    ///
    /// Example: `--chimeric-pairs 0:1` for human + HPV integration analysis.
    ///          `--chimeric-pairs 0:1 --chimeric-pairs 1:2` for human + HPV + mouse
    ///          where HPV can integrate into human and human+HPV tissue is xenografted
    ///          in mouse.  Pairs not listed compete normally in the tournament.
    #[arg(long, num_args = 0..)]
    pub(crate) chimeric_pairs: Vec<String>,

    /// Human-readable labels for each alignment stream (positional: stream 0, 1, …).
    ///
    /// Used as the value of the `XC:Z:` SAM aux tag written to chimeric reads.
    /// The tag on a read from stream N reads `XC:Z:<label of the other stream>`.
    /// Defaults to `stream_N` when not supplied.
    ///
    /// Example: `--stream-labels human hpv mouse`
    #[arg(long, num_args = 0..)]
    pub(crate) stream_labels: Vec<String>,

    /// Parsed and validated chimeric stream pairs.
    /// Stored in canonical order (lower index first) so lookups are O(pairs).
    #[arg(skip)]
    pub(crate) parsed_chimeric_pairs: Vec<[usize; 2]>,
}

impl Config {
    pub(super) fn validate_and_init(&mut self) -> Result<(), Error> {
        let aln_count = self.alignment.len();
        // -- Chimeric pair parsing --------------------------------------------
        let mut parsed_chimeric_pairs: Vec<[usize; 2]> = Vec::new();
        for raw in &self.chimeric_pairs {
            let (a_str, b_str) = raw
                .split_once(':')
                .ok_or(Error::ChimericPairsInvalidFormat { raw: raw.clone() })?;
            let a =
                a_str
                    .trim()
                    .parse::<usize>()
                    .map_err(|_| Error::ChimericPairsInvalidIndex {
                        index_str: a_str.to_string(),
                    })?;
            let b =
                b_str
                    .trim()
                    .parse::<usize>()
                    .map_err(|_| Error::ChimericPairsInvalidIndex {
                        index_str: b_str.to_string(),
                    })?;

            if a == b {
                return Err(Error::ChimericPairsIdenticalIndices { raw: raw.clone() });
            }
            if self.matching_algorithm != MatchingAlgorithm::Namesorted {
                return Err(Error::ChimericPairsRequiresNamesorted);
            }
            // Canonical order: lower index first.
            parsed_chimeric_pairs.push([a.min(b), a.max(b)]);
        }
        // Deduplicate.
        parsed_chimeric_pairs.sort_unstable();
        parsed_chimeric_pairs.dedup();
        self.parsed_chimeric_pairs = parsed_chimeric_pairs;

        if !self.parsed_chimeric_pairs.is_empty() {
            tracing::info!(
                pairs = ?self.parsed_chimeric_pairs,
                "Chimeric pair detection enabled"
            );
        }

        if self.input_format == AlnFormat::Cram && self.reference.is_none() {
            return Err(Error::CramRequiresReference);
        }
        if self.input_format == AlnFormat::Cram
            && self.matching_algorithm != MatchingAlgorithm::Namesorted
        {
            return Err(Error::CramNamesortedOnly);
        }
        // Stdin: path "-" implies SAM unless explicitly overridden.
        if self.alignment.iter().any(|p| p == "-") {
            if self.alignment.iter().filter(|p| p.as_str() == "-").count() > 1 {
                return Err(Error::MultipleStdinStreams);
            }
            if self.matching_algorithm != MatchingAlgorithm::Namesorted {
                return Err(Error::StdinRequiresNamesorted);
            }
        }

        // Reject multi-threaded modes for non-namesorted algorithms.
        if self.matching_algorithm == MatchingAlgorithm::Namesorted {
            if !self.ambiguous_regions.is_empty() || !self.diagnostic_variants.is_empty() {
                return Err(Error::NamesortedUnsupportedOptions);
            }
            // Resolve score_threads = 0 → available logical CPUs with a max of 16.
            if self.score_threads == 0 {
                self.score_threads = std::thread::available_parallelism()
                    .map(|n| n.get().min(16))
                    .unwrap_or(1);
                tracing::info!(
                    score_threads = self.score_threads,
                    "score_threads=0: using all available logical CPUs"
                );
            }
        } else {
            if self.score_threads != 1 {
                return Err(Error::MultiThreadedScoringRequiresNamesorted);
            }
        }

        // 1. Guard against single-stream niche execution accidents first
        if aln_count == 1 {
            if !self.single_alignment_mode {
                return Err(Error::SingleStreamMissingFlag);
            }
            if self.read_from_stdin {
                return Err(Error::SingleStreamStdinUnsupported);
            }
            if self.matching_algorithm != MatchingAlgorithm::Namesorted {
                return Err(Error::SingleStreamRequiresNamesorted);
            }
        } else {
            if self.single_alignment_mode {
                return Err(Error::SingleAlignmentModeExpectsOneStream { count: aln_count });
            }
            if aln_count < 2 {
                return Err(Error::InsufficientAlignmentStreams { count: aln_count });
            }
            if aln_count > 2 {
                if self.matching_algorithm != MatchingAlgorithm::Namesorted {
                    return Err(Error::MultiStreamRequiresNamesorted);
                }
                if aln_count > MAX_STREAMS {
                    return Err(Error::MaxStreamsExceeded {
                        count: aln_count,
                        max: MAX_STREAMS,
                    });
                }
            }
        }
        if self.merged_output.is_some()
            && (!self.output.is_empty()
                || !self.discarded_output.is_empty()
                || !self.ambiguous_output.is_empty())
            {
                return Err(Error::MergedOutputConflict);
            }
        if self.ambiguous_regions.len() > 2 {
            return Err(Error::TooManyAmbiguousRegionsFiles {
                count: self.ambiguous_regions.len(),
            });
        }
        if self.diagnostic_variants.len() > 2 {
            return Err(Error::TooManyDiagnosticVariantsFiles {
                count: self.diagnostic_variants.len(),
            });
        }
        // Determine effective dimensions (logical comparisons)
        let logical_len = if aln_count == 1 { 2 } else { aln_count };

        // 2. Parse, validate, and normalize variant arrays to size matching comparisons
        let mut normalized_samples = vec![PathBuf::new(); logical_len];
        let mut normalized_populations = vec![PathBuf::new(); logical_len];

        let mut stream_has_variants = vec![false; logical_len];

        for (i, arg) in self.sample_variants.iter().enumerate() {
            let (idx, path) = Self::parse_variant_string(arg, i)?;
            if idx >= logical_len {
                return Err(Error::SampleVariantIndexOutOfBounds {
                    idx,
                    max: logical_len,
                });
            }
            normalized_samples[idx] = path;
            stream_has_variants[idx] = true;
        }

        for (i, arg) in self.population_variants.iter().enumerate() {
            let (idx, path) = Self::parse_variant_string(arg, i)?;
            if idx >= logical_len {
                return Err(Error::PopulationVariantIndexOutOfBounds {
                    idx,
                    max: logical_len,
                });
            }
            normalized_populations[idx] = path;
            stream_has_variants[idx] = true;
        }

        if aln_count == 1
            && (!stream_has_variants[0] || !stream_has_variants[1]) {
                return Err(Error::SingleStreamMissingVariantProfiles);
            }

        self.sample_variants = normalized_samples
            .into_iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();
        self.population_variants = normalized_populations
            .into_iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect();

        // 3. Output bounds validation
        if self.output.len() > logical_len {
            return Err(Error::TooManyOutputPaths {
                count: self.output.len(),
                max: logical_len,
            });
        }
        if self.discarded_output.len() > logical_len {
            return Err(Error::TooManyDiscardedOutputPaths {
                count: self.discarded_output.len(),
                max: logical_len,
            });
        }

        if aln_count == 1 {
            if self.ambiguous_output.len() > 1 {
                return Err(Error::SingleStreamTooManyAmbiguousOutputs {
                    count: self.ambiguous_output.len(),
                });
            }
        } else {
            if self.ambiguous_output.len() > logical_len {
                return Err(Error::TooManyAmbiguousOutputPaths {
                    count: self.ambiguous_output.len(),
                    max: logical_len,
                });
            }
        }

        if self.gap_open > 0.0 {
            self.gap_open = -self.gap_open;
        }
        if self.gap_extend > 0.0 {
            self.gap_extend = -self.gap_extend;
        }

        if self.gap_open == 0.0 || self.mismatch_penalty <= 0.0 {
            return Err(Error::InvalidPenalties);
        }

        tracing::debug!(
            threads = self.threads,
            score_threads = self.score_threads,
            alignments = aln_count,
            gap_open = self.gap_open,
            gap_extend = self.gap_extend,
            mismatch = self.mismatch_penalty,
            "Configuration validated"
        );

        Ok(())
    }

    /// Parse `"<idx>:<path>"` or fall back to `(default_idx, path)`.
    fn parse_variant_string(arg: &str, default_idx: usize) -> Result<(usize, PathBuf), Error> {
        if let Some((idx_str, path_str)) = arg.split_once(':')
            && let Ok(idx) = idx_str.parse::<usize>()
        {
            return Ok((idx, PathBuf::from(path_str)));
        }
        Ok((default_idx, PathBuf::from(arg)))
    }

    /// Build a [`Penalty`] from the current penalty parameters.
    pub(super) fn to_penalties(&self) -> Penalty {
        let mut error_prob = [0.0_f64; MAX_Q];
        for (q, item) in error_prob.iter_mut().enumerate() {
            *item = 10f64.powf(-(q as f64) / 10.0);
        }

        let mut log_likelihood_match = [0.0_f64; MAX_Q];
        for (q, item) in log_likelihood_match.iter_mut().enumerate() {
            *item = (1.0 - error_prob[q]).log10();
        }

        let reference_penalty = 4.0;
        let scaling_factor = self.mismatch_penalty / reference_penalty;

        let mut log_likelihood_mismatch = [0.0f64; MAX_Q];
        for (q, item) in log_likelihood_mismatch.iter_mut().enumerate() {
            *item = -(q as f64) / 10.0 * scaling_factor;
        }

        Penalty {
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            log_likelihood_mismatch,
            log_likelihood_match,
            chimeric_junction_penalty: self.gap_open
                + (self.chimeric_junction_bases as f64) * self.gap_extend,
        }
    }
}

#[cfg(test)]
pub(crate) mod tests;
