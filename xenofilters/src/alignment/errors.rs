use thiserror::Error;

#[derive(Debug, Error)]
pub(crate) enum AlignmentError {
    #[error(
        "MD/CIGAR inconsistency: A CIGAR deletion ('D') was not matched by a corresponding MD deletion ('^')."
    )]
    MissingMdDeletion,

    #[error("MD/CIGAR inconsistency: Excess mismatches in MD tag after processing CIGAR.")]
    MdCigarMismatch,

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

    #[error(transparent)]
    MdOpIterator(#[from] MdOpIteratorError),
}

#[derive(Debug, Error)]
pub(crate) enum MdOpIteratorError {
    #[error("Wrong MD tag type found")]
    BadMdTag,

    #[error("MD parsing error: invalid byte '{0}'")]
    MdParse(u8),

    #[error(transparent)]
    MdError(#[from] std::io::Error),
}
