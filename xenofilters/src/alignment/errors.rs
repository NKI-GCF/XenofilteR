use thiserror::Error;

#[derive(Debug, Error)]
pub enum AlignmentError {
    #[error(
        "MD/CIGAR inconsistency: A CIGAR deletion ('D') was not matched by a corresponding MD deletion ('^')."
    )]
    MissingMdDeletion,

    #[error("MD/CIGAR inconsistency: Excess mismatches in MD tag after processing CIGAR.")]
    MdCigarMismatch,

    #[error("Operation not implemented")]
    UnImplemented,

    #[error(transparent)]
    Aux(#[from] rust_htslib::errors::Error),

    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),
}

#[derive(Debug, Error)]
pub enum PrepareError {
    #[error(transparent)]
    Alignment(#[from] AlignmentError),

    #[error(transparent)]
    MdOpIterator(#[from] MdOpIteratorError),
}

#[derive(Debug, Error)]
pub enum MdOpIteratorError {
    #[error("Wrong MD tag type found")]
    BadMdTag,

    #[error("MD parsing error: invalid character '{0}'")]
    MdParse(char),

    #[error(transparent)]
    Aux(#[from] rust_htslib::errors::Error),
}
