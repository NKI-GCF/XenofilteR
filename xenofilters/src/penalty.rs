pub(super) const MAX_Q: usize = 93;

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
}
