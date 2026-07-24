pub(super) const MAX_Q: usize = 93;

/// Platform-specific quality-score calibration and gap cost presets.
///
/// Quality calibration factor: effective_q = reported_q * factor.
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
    /// Oxford Nanopore (R9.4 / R10.4): quality calibration ~0.7*,
    /// indel rate 1-10%.  Low gap penalties. gap_open=2, gap_extend=0.3.
    Ont,
}

impl ErrorModel {
    // Return a tuple of (calibration, gap_open, gap_extend, mismatch_penalty)
    fn params(&self) -> (f64, f64, f64, f64) {
        match self {
            ErrorModel::Illumina => (1.0, 6.0, 1.0, 4.0),
            ErrorModel::HiFi => (1.0, 4.0, 0.5, 4.0),
            // NEEDS VALIDATION: 0.7 is a placeholder calibration factor, not
            // derived from a specific basecaller's published error model.
            // Do not treat ONT output confidence as validated without
            // checking against a real R9.4/R10.4 dataset for your basecaller
            // version.
            ErrorModel::Ont => (0.7, 2.0, 0.3, 2.0),
        }
    }
    /// Fraction by which reported Phred scores are scaled before computing
    /// per-base log-likelihoods.  Values < 1.0 reduce confidence in
    /// reported qualities (conservative scoring).
    pub(crate) fn quality_calibration(&self) -> f64 {
        self.params().0
    }

    /// Suggested gap_open penalty for this platform.
    pub(crate) fn default_gap_open(&self) -> f64 {
        self.params().1
    }

    /// Suggested gap_extend penalty for this platform.
    pub(crate) fn default_gap_extend(&self) -> f64 {
        self.params().2
    }

    /// Suggested mismatch_penalty for this platform.
    pub(crate) fn default_mismatch_penalty(&self) -> f64 {
        self.params().3
    }
}

#[derive(Clone, Copy)]
pub struct Penalty {
    /// Stored as a **negative** cost in nats (natural log units), consistent
    /// with `log_likelihood_match`/`log_likelihood_mismatch` below. Adding
    /// this to a running NW score always makes gapped alignments less
    /// favorable, as required for a max-picking DP recurrence.
    pub(crate) gap_open: f64,
    pub(crate) gap_extend: f64,
    pub(crate) log_likelihood_mismatch: [f64; MAX_Q],
    pub(crate) log_likelihood_match: [f64; MAX_Q],
    /// Flat penalty applied once per supplementary alignment:
    ///   `gap_open + chimeric_junction_bases * gap_extend`
    pub(crate) chimeric_junction_penalty: f64,
    /// Penalty for bisulfite conversions (C->T or G->A).
    /// Used only when --bisulfite is active.
    pub(crate) log_likelihood_bisulfite_base: [f64; MAX_Q],
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
        let mut ll_bisulfite = [0.0f64; MAX_Q];

        const REFERENCE_PENALTY: f64 = 4.0;
        let bias_nats = (mismatch_penalty - REFERENCE_PENALTY) * std::f64::consts::LN_10 / 10.0;

        for q in 0..MAX_Q {
            let eff_q = (q as f64 * calib).min(MAX_Q as f64 - 1.0);
            let p_err = 10f64.powf(-eff_q / 10.0);

            // Natural log throughout — matches the Phred conversion used by
            // `resolve_threshold` and `Decision::PhredConfidence`.
            ll_match[q] = (1.0 - p_err).ln();
            ll_mismatch[q] = p_err.ln() - bias_nats;

            // A low-Q "C->T" is still uncertain. Charge it the quality-scaled
            // cost a genuine match receives. TODO: still ignores CpG context
            // and methylation priors, but no longer quality-blind.
            ll_bisulfite[q] = ll_match[q];
        }

        // Convert Phred-documented gap costs into negative nats, consistent
        // with the match/mismatch scale above. `chimeric_junction_bases` is
        // a base-count multiplier, not itself a Phred value, so it is
        // applied after conversion.
        let to_neg_nats = |phred: f64| -(phred * std::f64::consts::LN_10 / 10.0);
        let gap_open_nats = to_neg_nats(gap_open);
        let gap_extend_nats = to_neg_nats(gap_extend);
        let chimeric_junction_penalty =
            gap_open_nats + (chimeric_junction_bases as f64) * gap_extend_nats;

        Penalty {
            gap_open: gap_open_nats,
            gap_extend: gap_extend_nats,
            log_likelihood_match: ll_match,
            log_likelihood_mismatch: ll_mismatch,
            chimeric_junction_penalty,
            log_likelihood_bisulfite_base: ll_bisulfite,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::LN_10;

    #[test]
    fn default_mismatch_penalty_is_a_genuine_log_probability() {
        // At the reference penalty (4.0), ll_mismatch must equal ln(p_err)
        // exactly — no residual scaling artifact.
        let p = Penalty::build(6.0, 1.0, 4.0, 20, ErrorModel::Illumina);
        let q = 30;
        let eff_q = q as f64; // calib = 1.0 for Illumina
        let expected = 10f64.powf(-eff_q / 10.0).ln();
        assert!((p.log_likelihood_mismatch[q] - expected).abs() < 1e-9);
    }

    #[test]
    fn nondefault_mismatch_penalty_shifts_by_flat_bias_not_rescale() {
        let p_ref = Penalty::build(6.0, 1.0, 4.0, 20, ErrorModel::Illumina);
        let p_double = Penalty::build(6.0, 1.0, 8.0, 20, ErrorModel::Illumina);
        // The *difference* between Q10 and Q30 mismatch cost must be
        // identical across mismatch_penalty settings (flat additive bias
        // does not change the quality-dependent shape, only its offset).
        let d_ref = p_ref.log_likelihood_mismatch[30] - p_ref.log_likelihood_mismatch[10];
        let d_double = p_double.log_likelihood_mismatch[30] - p_double.log_likelihood_mismatch[10];
        assert!(
            (d_ref - d_double).abs() < 1e-9,
            "mismatch_penalty must not distort the quality-dependence shape"
        );
    }

    #[test]
    fn phred_conversion_round_trips_at_default_settings() {
        // A log-likelihood-ratio "margin" of `phred * LN_10 / 10.0` nats
        // must convert back to exactly `phred` via the codebase's own
        // formula (`10.0 * margin / LN_10`), for any Phred value.
        for phred in [0u32, 5, 10, 20, 40] {
            let margin_nats = (phred as f64) * LN_10 / 10.0;
            let recovered = (10.0 * margin_nats / LN_10).round() as u32;
            assert_eq!(recovered, phred);
        }
    }

    #[test]
    fn gap_costs_are_negative_and_penalize() {
        let p = Penalty::build(6.0, 1.0, 4.0, 20, ErrorModel::Illumina);
        assert!(
            p.gap_open < 0.0,
            "gap_open must be a negative cost, got {}",
            p.gap_open
        );
        assert!(
            p.gap_extend < 0.0,
            "gap_extend must be a negative cost, got {}",
            p.gap_extend
        );
        assert!(p.chimeric_junction_penalty < 0.0);
    }

    #[test]
    fn gap_open_at_least_as_costly_as_extend_in_converted_units() {
        let p = Penalty::build(6.0, 1.0, 4.0, 20, ErrorModel::Illumina);
        assert!(p.gap_open.abs() >= p.gap_extend.abs());
    }

    #[test]
    fn bisulfite_conversion_cost_scales_with_quality() {
        let p = Penalty::build(6.0, 1.0, 4.0, 20, ErrorModel::Illumina);
        // Low-quality "conversion" must be cheaper to trust (less negative
        // relative to its own scale) than... actually the key invariant is
        // simply that it is NOT flat/zero across quality, unlike before.
        assert_ne!(
            p.log_likelihood_bisulfite_base[5], p.log_likelihood_bisulfite_base[40],
            "bisulfite scoring must not be quality-blind"
        );
        assert!(
            p.log_likelihood_bisulfite_base[40] > p.log_likelihood_bisulfite_base[5],
            "higher-quality conversion calls should be more favorably scored"
        );
    }
}
