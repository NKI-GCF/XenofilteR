use thiserror::Error;
use noodles::sam::alignment::record::cigar::Op;
use crate::alignment::ops::MdOp;

#[derive(Debug, Error)]
pub(crate) enum AlignmentError {
    #[error("Op inconsistency: cigar: ({0:?}) and md: ({1:?})")]
    MdCigMis(Option<Op>, Option<MdOp>),

    #[error("Operation not implemented")]
    UnImplemented,

    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),

    #[error(transparent)]
    MdError(#[from] std::io::Error),
}

#[derive(Debug, Error)]
pub(crate) enum PrepareError {
    #[error(transparent)]
    Alignment(#[from] AlignmentError),

    #[error("Wrong MD tag type found")]
    BadMdTag,

    #[error(transparent)]
    MdError(#[from] std::io::Error),

    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),
}
