use thiserror::Error;

#[derive(Debug, Error, PartialEq)]
pub enum AlignmentError {
    #[error(
        "MD/CIGAR inconsistency: Encountered an MD deletion ('^') during a CIGAR match ('M'/'='/'X')."
    )]
    MismatchedDeletion,

    #[error(
        "MD/CIGAR inconsistency: The MD tag ended before the CIGAR string was fully processed."
    )]
    UnexpectedMdEnd,

    #[error(
        "MD/CIGAR inconsistency: A CIGAR deletion ('D') was not matched by a corresponding MD deletion ('^')."
    )]
    MissingMdDeletion,

    #[error("Quality score index out of bounds")]
    QualIndexOutOfBounds,

    #[error("MD/CIGAR inconsistency: Excess mismatches in MD tag after processing CIGAR.")]
    MdCigarMismatch,

    #[error("Unexpected translocation operation in CIGAR string")]
    UnexpectedTranslocate,
}

#[derive(Debug, Error)]
pub enum PrepareError {
    #[error("No MD tag found")]
    NoMdTag,

    #[error("Aux error {0}")]
    AuxError(String),

    #[error(transparent)]
    AlignmentError(#[from] AlignmentError),
}
