pub(super) mod core;
pub(super) mod io;
pub(super) mod score;
pub(super) mod chimeric;

#[cfg(not(test))]
pub(super) mod ordering;
#[cfg(test)]
pub(crate) mod ordering;

#[cfg(test)]
mod tests;

// Re-export the main type and the most important items
pub(crate) use core::{LineByLine, Scratch, MAX_STREAMS, READ_CT};
pub(crate) use chimeric::{ChimericDecision, detect_chimeric_event};
