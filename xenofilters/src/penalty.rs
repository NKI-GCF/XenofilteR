use std::f64::consts::LN_10;

pub(super) const MAX_Q: usize = 93;

/// Platform-specific quality-score calibration and gap cost presets.
///
/// Quality calibration factor: effective_q = reported_q × factor.
/// ONT quality scores are systematically inflated; a factor of 0.7
/// gives empirically better-calibrated effective error rates for R10.
/// HiFi qualities are accurate but gap costs should reflect the lower
/// indel rate of CCS chemistry.
///
/// All values can be individually overridden with the corresponding
/// CLI flags; the model only sets defaults.
#[derive(Clone, Copy, Debug, PartialEq, clap::ValueEnum, Default)]
pub enum ErrorModel {
    /// Standard Illumina short-read (SBS): qualities well-calibrated,
    /// indel rate ~0.1%.  Default gap_open=6, gap_extend=1.
    #[default]
    Illumina,
    /// PacBio HiFi (CCS): qualities highly accurate, SNR-based;
    /// indel rate ~0.01%.  Tighter gap penalties. gap_open=4, gap_extend=0.5.
    HiFi,
    /// Oxford Nanopore (R9.4 / R10.4): quality calibration ~0.7×,
    /// indel rate 1–10%.  Low gap penalties. gap_open=2, gap_extend=0.3.
    Ont,
}

impl ErrorModel {
    /// Fraction by which reported Phred scores are scaled before computing
    /// per-base log-likelihoods.  Values < 1.0 reduce confidence in
    /// reported qualities (conservative scoring).
    pub(crate) fn quality_calibration(&self) -> f64 {
        match self {
            ErrorModel::Illumina => 1.0,
            ErrorModel::HiFi => 1.0,
            ErrorModel::Ont => 0.7, // R10.4 empirical; adjust per basecaller
        }
    }

    /// Suggested gap_open penalty for this platform.
    pub(crate) fn default_gap_open(&self) -> f64 {
        match self {
            ErrorModel::Illumina => 6.0,
            ErrorModel::HiFi => 4.0,
            ErrorModel::Ont => 2.0,
        }
    }

    /// Suggested gap_extend penalty for this platform.
    pub(crate) fn default_gap_extend(&self) -> f64 {
        match self {
            ErrorModel::Illumina => 1.0,
            ErrorModel::HiFi => 0.5,
            ErrorModel::Ont => 0.3,
        }
    }

    /// Suggested mismatch_penalty for this platform.
    pub(crate) fn default_mismatch_penalty(&self) -> f64 {
        match self {
            ErrorModel::Illumina => 4.0,
            ErrorModel::HiFi => 4.0,
            ErrorModel::Ont => 2.0, // mismatches are less surprising in ONT
        }
    }
}

#[derive(Clone, Copy)]
pub struct Penalty {
    pub(crate) gap_open: f64,
    pub(crate) gap_extend: f64,
    pub(crate) log_likelihood_mismatch: [f64; MAX_Q],
    pub(crate) log_likelihood_match: [f64; MAX_Q],
    /// Flat penalty applied once per supplementary alignment:
    ///   `gap_open + chimeric_junction_bases × gap_extend`
    ///
    /// Using a constant base-length (rather than the record's actual non-clipped
    /// length) decouples the penalty from read-length variation across
    /// supplementary records and makes it independently tunable.
    /// Both terms are negative after sign normalisation in `Config::to_penalties`.
    pub(crate) chimeric_junction_penalty: f64,
    /// Zero-cost penalty for bisulfite conversions (C→T or G→A).
    /// Used only when --bisulfite is active.
    pub(crate) log_likelihood_bisulfite: f64,
}

impl Penalty {
    pub(crate) fn build(
        gap_open: f64,
        gap_extend: f64,
        mismatch_penalty: f64,
        chimeric_junction_bases: u32,
        error_model: ErrorModel,
    ) -> Self {
        let calib = error_model.quality_calibration();
        let mut ll_match = [0.0f64; MAX_Q];
        let mut ll_mismatch = [0.0f64; MAX_Q];

        let reference_penalty = 4.0;
        let scaling_factor = mismatch_penalty / reference_penalty;

        for q in 0..MAX_Q {
            // Apply platform-specific quality calibration
            let eff_q = (q as f64 * calib).min(MAX_Q as f64 - 1.0);

            // 1. Match Log-Likelihood (Base-10)
            let p_err = 10f64.powf(-eff_q / 10.0);
            ll_match[q] = (1.0 - p_err).log10();

            // 2. Mismatch Log-Likelihood (Base-10 back-compat formula)
            ll_mismatch[q] = (-eff_q / 10.0) * scaling_factor;
        }

        let chimeric_junction_penalty = gap_open + (chimeric_junction_bases as f64) * gap_extend;

        Penalty {
            gap_open, // Do not force negative here; let Config/algorithms handle signs
            gap_extend,
            log_likelihood_match: ll_match,
            log_likelihood_mismatch: ll_mismatch,
            chimeric_junction_penalty,
            log_likelihood_bisulfite: 0.0,
        }
    }
}
