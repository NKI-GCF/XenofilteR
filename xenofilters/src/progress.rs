//! One-line stderr progress updated every 100 k fragments.

use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};

const REPORT_INTERVAL: u64 = 100_000;

const STYLE: &str = "[xenofilters] {human_pos} frags  {wide_bar:.cyan/blue}  \
     {per_sec}  {elapsed_precise} elapsed  {eta_precise} ETA";

pub(crate) struct ProgressReporter {
    bar: ProgressBar,
}

impl ProgressReporter {
    pub(crate) fn new() -> Self {
        let bar = ProgressBar::with_draw_target(
            None,
            ProgressDrawTarget::stderr_with_hz(4), // 4 redraws/sec
        );
        bar.set_style(
            ProgressStyle::with_template(STYLE)
                .expect("progress template is valid")
                .progress_chars("#."),
        );
        Self { bar }
    }

    /// Call once per completed fragment. Returns true on every draw tick.
    pub(crate) fn tick(&mut self) -> bool {
        self.bar.inc(1);
        self.bar.is_finished() // always false; keeps caller API stable
    }

    pub(crate) fn finish(&self) {
        self.bar.finish_and_clear();
    }
}
