mod errors;
mod iterator;
mod ops;
mod prepared;

pub use errors::{AlignmentError, PrepareError};
pub use iterator::{AlignmentCompareIterator, AlignmentIterator};
pub use ops::{AlignmentOp, AlnCmpOp, MdOp, MdOpIterator};
pub use prepared::{PreparedAlignmentPair, PreparedAlignmentPairIter};
