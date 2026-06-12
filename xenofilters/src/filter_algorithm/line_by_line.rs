pub(super) mod core;
pub(super) mod io;
pub(super) mod score;

#[cfg(not(test))]
pub(super) mod ordering;
#[cfg(test)]
pub(crate) mod ordering;

#[cfg(test)]
mod tests;


// Re-export the main type and the most important items
pub(crate) use core::LineByLine;
pub(crate) use score::NeedlemanWunsch;
