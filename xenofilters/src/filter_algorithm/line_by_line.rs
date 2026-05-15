pub(super) mod core;
pub(super) mod io;
pub(super) mod ordering;

#[cfg(test)]
mod tests;

// Re-export the main type and the most important items
pub(crate) use core::LineByLine;
