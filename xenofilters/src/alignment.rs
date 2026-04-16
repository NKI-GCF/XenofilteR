mod errors;
mod mdopiterator;
mod ops;
mod prepared;
mod fragment;
use smallvec::SmallVec;

pub(crate) use errors::{AlignmentError, MdOpIteratorError, PrepareError};
pub(crate) use mdopiterator::{MdOp, MdOpIterator};
pub(crate) use ops::{UnifiedOp, UnifiedOpIterator};
pub(crate) use prepared::{PreparedAlignmentPairIter, stringify_record};
pub(crate) use fragment::{build_fragment, Segment};
pub(crate) type SegmentVec<'a> = SmallVec<[Segment<'a>; 2]>;

#[cfg(test)]
pub(crate) use prepared::PreparedAlignmentPair;

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    pub(crate) use ops::tests::*;
    pub(crate) use prepared::tests::*;
    pub(crate) use fragment::tests::*;
}
