//! One-line stderr progress updated every 100 k fragments.

use std::io::Write;

const REPORT_INTERVAL: u64 = 100_000;

pub(crate) struct ProgressReporter {
    fragments:  u64,
    start:      std::time::Instant,
}

impl ProgressReporter {
    pub(crate) fn new() -> Self {
        Self { fragments: 0, start: std::time::Instant::now() }
    }

    /// Call once per completed fragment. Returns true on every report tick.
    pub(crate) fn tick(&mut self) -> bool {
        self.fragments += 1;
        if self.fragments % REPORT_INTERVAL != 0 {
            return false;
        }
        let elapsed = self.start.elapsed().as_secs_f64();
        let rate    = self.fragments as f64 / elapsed.max(1e-9) / 1_000.0;
        eprint!(
            "\r[xenofilters] {:>12} fragments  {:.1} k/s  {:.1}s elapsed   ",
            self.fragments, rate, elapsed
        );
        let _ = std::io::stderr().flush();
        true
    }

    pub(crate) fn finish(&self) {
        eprintln!();
    }

    pub(crate) fn count(&self) -> u64 { self.fragments }
}
