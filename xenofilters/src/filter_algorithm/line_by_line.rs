pub(super) mod chimeric;
pub mod core;
pub(super) mod io;
pub(super) mod parallel;

#[cfg(not(test))]
pub(super) mod ordering;
#[cfg(test)]
pub(crate) mod ordering;

#[cfg(test)]
mod tests;

// Re-export the main type and the most important items
pub(crate) use chimeric::{ChimericDecision, detect_chimeric_event};
pub(crate) use core::{COUNTER_STRIDE, LineByLine, MAX_STREAMS, READ_CT, Scratch};
pub(crate) use ordering::apply_decision_tag;
