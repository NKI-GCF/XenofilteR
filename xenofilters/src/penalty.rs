pub(super) const MAX_Q: usize = 93;

pub(super) struct Penalty {
    pub(crate) gap_open: f64,
    pub(crate) gap_extend: f64,
    pub(crate) log_likelihood_mismatch: [f64; MAX_Q],
    pub(crate) log_likelihood_match: [f64; MAX_Q],
}
