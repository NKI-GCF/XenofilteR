// src/reporting.rs -- new shared module, no longer on Config
use crate::filter_algorithm::line_by_line::COUNTER_STRIDE;
use crate::Error;
use std::path::Path;

pub(crate) struct RunReport<'a> {
    pub(crate) counters: &'a [u64],
    pub(crate) backend: &'static str,
    pub(crate) warn_ambig_fraction: f64,
    pub(crate) is_pass2: bool,
    pub(crate) threshold_phred: u32,
}

impl<'a> RunReport<'a> {
    pub(crate) fn print_routing_counters(&self) {
        let stream_count = self.counters.len() / COUNTER_STRIDE;
        for nr in 0..stream_count {
            let b = nr * COUNTER_STRIDE;
            tracing::info!(
                stream = nr,
                backend = self.backend,
                discard = self.counters[b],
                out = self.counters[b + 1],
                ambiguous = self.counters[b + 2],
                chimeric = self.counters[b + 3],
                "Stream summary"
            );
        }
        let total: u64 = self.counters.iter().sum();
        if total == 0 {
            return;
        }
        let ambiguous: u64 = self
            .counters
            .chunks_exact(COUNTER_STRIDE)
            .map(|c| c[2])
            .sum();
        let frac = ambiguous as f64 / total as f64;
        if frac > self.warn_ambig_fraction {
            let threshold_phred = match self.threshold_phred {
                u32::MAX => {
                    if self.is_pass2 {
                        0
                    } else {
                        10
                    }
                }
                p => p,
            };
            tracing::warn!(
                ambiguous_pct = format!("{:.1}", frac * 100.0),
                threshold_phred,
                "Ambiguous fraction exceeds warning level. \
                 Consider pass 2 with BQSR or lowering --ambiguous-threshold."
            );
        }
    }

    pub(crate) fn write_stats_if_configured(
        &self,
        path: Option<&Path>,
        labels: &[String],
        sample: &str,
    ) -> Result<(), Error> {
        let Some(path) = path else { return Ok(()) };
        let stream_count = self.counters.len() / COUNTER_STRIDE;
        crate::stats::write_stats(path, self.counters, stream_count, labels, sample)
    }
}
