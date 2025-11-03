mod errors;
mod iterator;
mod ops;
mod prepared;
mod cigar_operation;
mod stitched_alignment;

pub use errors::{AlignmentError, PrepareError};
pub use iterator::{AlignmentCompareIterator, AlignmentIterator};
pub use ops::{AlignmentOp, AlnCmpOp, MdOp, MdOpIterator, UnifiedOp, lift_alignment_ops};
pub use prepared::{PreparedAlignmentPair, PreparedAlignmentPairIter};
pub use cigar_operation::{CigarOp, CigarOpIterator};
pub use stitched_alignment::{StitchedAlignment, stitch_alignment_segments};

