mod errors;
mod iterator;
mod ops;
mod prepared;
mod stitched_alignment;

pub use errors::{AlignmentError, PrepareError};
pub use iterator::AlignmentCompareIterator;
pub use ops::{AlnCmpOp,  MdOpIterator, UnifiedOp, UnifiedOpIterator, lift_alignment_ops};
pub use prepared::{PreparedAlignmentPair, PreparedAlignmentPairIter};
