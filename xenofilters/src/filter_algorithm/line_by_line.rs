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
pub(crate) use chimeric::{detect_chimeric_event, ChimericDecision};
pub(crate) use core::{LineByLine, Scratch, MAX_STREAMS, READ_CT};
