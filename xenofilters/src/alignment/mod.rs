mod errors;
mod mdopiterator;
mod ops;
mod prepared;
mod stitched_alignment;

pub use errors::{AlignmentError, MdOpIteratorError, PrepareError};
pub use mdopiterator::{MdOp, MdOpIterator};
pub use ops::{UnifiedOp, UnifiedOpIterator};
pub use prepared::{PreparedAlignmentPair, PreparedAlignmentPairIter};
pub use stitched_alignment::stitched_fragment;

#[cfg(test)]
pub mod tests {
    use super::*;
    pub use ops::tests::*;
    pub use prepared::tests::*;
}
