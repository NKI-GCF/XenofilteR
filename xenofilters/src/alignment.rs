mod errors;
mod mdopiterator;
mod ops;
mod prepared;
mod fragment;

pub(crate) use errors::{AlignmentError, MdOpIteratorError, PrepareError};
pub(crate) use mdopiterator::{MdOp, MdOpIterator};
pub(crate) use ops::{UnifiedOp, UnifiedOpIterator};
pub(crate) use prepared::{PreparedAlignmentPairIter, stringify_record};
pub(crate) use fragment::build_fragment;

#[cfg(test)]
pub(crate) use prepared::PreparedAlignmentPair;

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    pub(crate) use ops::tests::*;
    pub(crate) use prepared::tests::*;
    pub(crate) use fragment::tests::*;
}
