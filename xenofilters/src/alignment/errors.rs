use thiserror::Error;
use noodles::sam::alignment::record::cigar::Op;

#[derive(Debug, Error)]
pub(crate) enum AlignmentError {
    #[error("Op inconsistency: cigar: ({0:?}) and md: ({1:?})")]
    MdCigMis(Option<Op>, Option<u8>),

    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),

    #[error(transparent)]
    MdError(#[from] std::io::Error),
}
