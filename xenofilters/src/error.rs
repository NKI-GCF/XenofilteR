use thiserror::Error;
use std::result::Result as StdResult;
use noodles::sam::alignment::record::cigar::Op;
use std::path::PathBuf;

/// A type alias for `Result` that uses `Error` as the error type.
pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub(crate) enum Error {
    // --- Existing Core Errors ---
    #[error("Op inconsistency: cigar: ({0:?}) and md: ({1:?})")]
    MdCigMis(Option<Op>, Option<u8>),

    #[error(transparent)]
    MdError(#[from] std::io::Error),

    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),

    // --- BAM & Stream Records Management ---
    #[error("{bam_str} has no records")]
    BamHasNoRecords { bam_str: String },

    #[error("Record has no name")]
    RecordHasNoName,

    #[error("Record has no read name")]
    RecordHasNoReadName,

    #[error("Cannot un-next more than one record")]
    CannotUnNext,

    #[error("No BAM reader available for seek")]
    NoBamReaderForSeek,

    #[error("No record at virtual offset {virtual_offset}")]
    NoRecordAtVirtualOffset { virtual_offset: u64 },

    #[error("MockStream: no record at virtual offset {virtual_offset}")]
    MockStreamNoRecordAtVirtualOffset { virtual_offset: u64 },

    #[error("BAM read error: {0}")]
    BamReadError(String),

    #[error("RecordBuf conversion: {0}")]
    RecordBufConversion(String),

    #[error("fetch_by_virtual_offset not supported for this stream type")]
    FetchByVirtualOffsetNotSupported,

    #[error("No stream {nr}")]
    NoStream { nr: usize },

    // --- CLI & Configuration Argument Handling ---
    #[error("--chimeric-pairs: '{index}' is not a valid stream index")]
    InvalidChimericPairsIndex { index: String },

    #[error("--chimeric-pairs: expected format 'A:B' (e.g. '0:1'), got '{raw}'")]
    InvalidChimericPairsFormat { raw: String },

    #[error("Gap open/mismatch penalties must be positive")]
    InvalidPenalties,

    // --- IO & Filesystem Path Verification ---
    #[error("Path '{path}' is not valid UTF-8", path = path.display())]
    InvalidPathUtf8 { path: PathBuf },

    #[error("Cannot create output file '{path}': {source}", path = path.display())]
    CreateOutputFileFailed { path: PathBuf, source: std::io::Error },

    #[error("Cannot open index {path}: {source}", path = path.display())]
    CannotOpenIndex { path: PathBuf, source: std::io::Error },

    #[error("Cannot read tabix index {path}: {source}", path = path.display())]
    CannotReadTabixIndex { path: PathBuf, source: std::io::Error },

    #[error("No tabix index found for {path} (tried .tbi and .<ext>.tbi)", path = path.display())]
    TabixIndexNotFound { path: PathBuf },

    #[error("Cannot open VCF {path}: {source}", path = path.display())]
    CannotOpenVcf { path: PathBuf, source: std::io::Error },

    #[error("Cannot open VCF/BCF {path}: {source}", path = path.display())]
    CannotOpenVcfBcf { path: PathBuf, source: std::io::Error },

    #[error("Failed to open VCF/BCF {path}: {source}", path = path.display())]
    FailedToOpenVcfBcf { path: PathBuf, source: std::io::Error },

    #[error("Cannot open BED file {path}: {source}", path = path.display())]
    CannotOpenBedFile { path: PathBuf, source: std::io::Error },

    // --- Format Parsing & Genomics Metadata Validation ---
    //#[error("Invalid region {region_str}: {source}")]
    //InvalidRegion { region_str: String, source: String },

    #[error("Tabix BED query failed: {0}")]
    TabixBedQueryFailed(String),

    #[error("Tabix VCF query failed: {0}")]
    TabixVcfQueryFailed(String),

    #[error("VCF header read error: {0}")]
    VcfHeaderReadError(String),

    #[error("Multiple ALT alleles not supported for population variants")]
    MultipleAltAllelesNotSupported,

    #[error("Missing AF tag or AF tag is not a float")]
    MissingOrInvalidAfTag,

    #[error("Multiple samples not supported")]
    MultipleSamplesNotSupported,

    #[error("Missing GQ tag or not an integer")]
    MissingOrInvalidGqTag,

    #[error("Missing GT tag or not an integer")]
    MissingOrInvalidGtTag,

    //#[error("BED parse error in {path}: {source}", path = path.display())]
    //BedParseError { path: PathBuf, source: String },

    #[error("BED start error: {0}")]
    BedStartError(String),

    #[error("BED record missing end field")]
    BedRecordMissingEnd,

    #[error("BED end error: {0}")]
    BedEndError(String),

    #[error("No sample data in record")]
    NoSampleData,

    #[error("Quality score index {nt_i} out of bounds for segment {seg_i}")]
    QualityScoreOutOfBounds { nt_i: usize, seg_i: usize },

    #[error("BCF contig index {chrom_idx} not in header")]
    BcfContigMissing { chrom_idx: usize },

    #[error("Invalid position {0}")]
    InvalidPosition(i64),

    #[error("Unknown CIGAR op {0}")]
    UnknownCigarOp(char),

    #[error("Invalid CIGAR character: {c}")]
    InvalidCigarChar { c: char },

    #[error("MD not UTF-8: {0}")]
    MdNotUtf8(#[from] std::string::FromUtf8Error),

    #[error("missing MD tag")]
    MissingMdTag,

    #[error("unexpected MD tag value type")]
    UnexpectedMdTagValueType,

    #[error("RG aux tag has unexpected type {0}; expected String")]
    UnexpectedRgTagType(String),

    // --- Structural Alignment Checks & Flags ---
    #[error("No ref seq id")]
    NoRefSeqId,

    #[error("No reference sequence ID")]
    NoReferenceSequenceId,

    #[error("No alignment start")]
    NoAlignmentStart,

    #[error("No alignment for index {aln_idx}")]
    NoAlignmentForIndex { aln_idx: usize },

    #[error("No alignment for index {nr}")]
    NoAlignmentForNr { nr: usize },

    #[error("No flags for record {idx}")]
    NoFlagsForRecord { idx: usize },

    #[error("No flags for record index {idx}")]
    NoFlagsForRecordIndex { idx: usize },

    #[error("No flags for record index {idx} in alignment {aln_idx}")]
    NoFlagsForRecordInAlignment { idx: usize, aln_idx: usize },

    #[error("Mapped record has no reference sequence ID")]
    MappedRecordNoReferenceSequenceId,

    #[error("Mapped record has no alignment start")]
    MappedRecordNoAlignmentStart,

    #[error("MdCigFlags missing for {idx}")]
    MdCigFlagsMissing { idx: usize },

    #[error("MdCigFlags missing for record index {idx}")]
    MdCigFlagsMissingForIndex { idx: usize },

    #[error("MdCigFlags missing for index {idx}")]
    MdCigFlagsMissingIndex { idx: usize },

    #[error("MdCigFlags consumed for {idx}")]
    MdCigFlagsConsumed { idx: usize },

    #[error("MdCigFlags already consumed for index {idx}")]
    MdCigFlagsAlreadyConsumed { idx: usize },

    // --- Scoring Workers & Concurrency ---
    #[error("Missing driving records for full scoring")]
    MissingDrivingRecords,

    #[error("Missing lookup records for full scoring")]
    MissingLookupRecords,

    //#[error("Score error stream {aln_idx}: {source}")]
    //ScoreStreamError { aln_idx: usize, source: String },

    #[error("Scorer worker exited unexpectedly")]
    ScorerWorkerExited,

    #[error("Work channel closed unexpectedly")]
    WorkChannelClosed,

    #[error("All scorer workers exited unexpectedly")]
    AllScorerWorkersExited,

    #[error("Scorer workers exited unexpectedly")]
    ScorerWorkersExited,

    #[error("Error scoring fragment for alignment {aln_idx}: {message}\n{state}")]
    FragmentScoringError { aln_idx: usize, message: String, state: String },

    #[error("Scoring error: {message}\n{state}")]
    ScoringError { message: String, state: String },

    #[error("BUG: unmapped record should already have been excluded")]
    UnmappedRecordInMdCigFlags
}
