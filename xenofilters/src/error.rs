use noodles::sam::alignment::record::cigar::Op;
use std::path::PathBuf;
use thiserror::Error;

/// A type alias for `Result` that uses `Error` as the error type.
pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    MdError(#[from] std::io::Error),

    #[error(transparent)]
    SerdeJson(#[from] serde_json::Error),

    // ---------------------------------------------------------------------------
    // BAM / SAM / BGZF I/O
    // ---------------------------------------------------------------------------
    #[error("{bam_str} has no records")]
    BamHasNoRecords { bam_str: String },

    #[error("Input read names do not have /1 or /2 suffixes, but strip_read_suffix is true.")]
    ReadNamesMissingSuffixes,

    #[error("Input read names have /1 or /2 suffixes, but strip_read_suffix is false.")]
    ReadNamesHaveSuffixes,

    #[error("All input BAMs must be either paired-end or single-end.")]
    MixedPairedAndSingleEnd,

    #[error("Record has no name")]
    RecordHasNoName,

    #[error("Record has no read name")]
    RecordHasNoReadName,

    #[error("Cannot un-next more than one record")]
    CannotUnNext,

    #[error("Input alignments must have the same read order.")]
    InputAlignmentsMustHaveSameReadOrder,

    #[error("No BAM reader available for seek")]
    NoBamReaderForSeek,

    #[error("Missing reference sequence for file {0}")]
    MissingReference(String),

    #[error("No record at virtual offset {virtual_offset}")]
    NoRecordAtVirtualOffset { virtual_offset: u64 },

    #[cfg(test)]
    #[error("MockStream: no record at virtual offset {virtual_offset}")]
    MockStreamNoRecordAtVirtualOffset { virtual_offset: u64 },

    #[error("--matching-algorithm hashlookup|collated requires exactly 2 alignment streams")]
    AlgoRequiresTwoStreams,

    #[error("BAM read error: {0}")]
    BamReadError(String),

    #[error("SAM file input (non-stdin) not yet implemented; use BAM or CRAM")]
    SamFileNotYetSupported,

    #[error("RecordBuf conversion: {0}")]
    RecordBufConversion(String),

    #[error("fetch_by_virtual_offset not supported for this stream type")]
    FetchByVirtualOffsetNotSupported,

    #[error(
        "Coordinate-sorted input detected; \
                     use --matching-algorithm hashlookup or collated."
    )]
    CoordinateSortedInputDetected,

    #[error("No stream {nr}")]
    NoStream { nr: usize },

    // ---------------------------------------------------------------------------
    // CLI / Configuration
    // ---------------------------------------------------------------------------
    #[error("--chimeric-pairs: '{index}' is not a valid stream index")]
    InvalidChimericPairsIndex { index: String },

    #[error("--chimeric-pairs: expected format 'A:B' (e.g. '0:1'), got '{raw}'")]
    InvalidChimericPairsFormat { raw: String },

    #[error("CRAM input requires --reference <fasta>")]
    CramRequiresReference,

    #[error("CRAM input is only supported with --matching-algorithm namesorted")]
    CramNamesortedOnly,

    #[error("stdin requires --input-format sam")]
    StdinRequiresSamFormat,

    #[error("Gap open/mismatch penalties must be positive")]
    InvalidPenalties,

    // --- IO & Filesystem Path Verification ---
    #[error("Path '{path}' is not valid UTF-8", path = path.display())]
    InvalidPathUtf8 { path: PathBuf },

    #[error("Cannot create output file '{path}': {source}", path = path.display())]
    CreateOutputFileFailed {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("Cannot open index {path}: {source}", path = path.display())]
    CannotOpenIndex {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("Cannot read tabix index {path}: {source}", path = path.display())]
    CannotReadTabixIndex {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("No tabix index found for {path} (tried .tbi and .<ext>.tbi)", path = path.display())]
    TabixIndexNotFound { path: PathBuf },

    #[error("Cannot open VCF {path}: {source}", path = path.display())]
    CannotOpenVcf {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("Cannot open VCF/BCF {path}: {source}", path = path.display())]
    CannotOpenVcfBcf {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("Failed to open VCF/BCF {path}: {source}", path = path.display())]
    FailedToOpenVcfBcf {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("Cannot open BED file {path}: {source}", path = path.display())]
    CannotOpenBedFile {
        path: PathBuf,
        source: std::io::Error,
    },

    // ---------------------------------------------------------------------------
    // Region / BED / Tabix
    // ---------------------------------------------------------------------------
    #[error("Invalid region {0}")]
    InvalidRegion(String),

    #[error("Tabix BED query failed: {0}")]
    TabixBedQueryFailed(String),

    #[error("Tabix VCF query failed: {0}")]
    TabixVcfQueryFailed(String),

    #[error("VCF header read error: {0}")]
    VcfHeaderReadError(String),

    // ---------------------------------------------------------------------------
    // Variant / VCF / BCF
    // ---------------------------------------------------------------------------
    #[error("Missing AF tag or AF tag is not a float")]
    MissingOrInvalidAfTag,

    #[error("Multiple samples not supported")]
    MultipleSamplesNotSupported,

    #[error("Missing GQ tag or not an integer")]
    MissingOrInvalidGqTag,

    #[error("Missing GT tag or not an integer")]
    MissingOrInvalidGtTag,

    #[error("BED parse error in {path}: {source}", path = path.display())]
    BedParseError {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("Invalid chromosome name {0}")]
    BedInvalidChromName(String),

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
    InvalidPosition(usize),

    // ---------------------------------------------------------------------------
    // CIGAR / MD alignment parsing
    // ---------------------------------------------------------------------------
    #[error("Op inconsistency: cigar: ({0:?}) and md: ({1:?})")]
    MdCigMis(Option<Op>, Option<u8>),

    #[error("Unknown CIGAR op {0}")]
    UnknownCigarOp(u32),

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

    #[error("alignment {i} still has reads")]
    AlignmentStillHasReads { i: usize },

    // ---------------------------------------------------------------------------
    // Concurrency / Worker channels
    // ---------------------------------------------------------------------------
    #[error("Missing driving records for full scoring")]
    MissingDrivingRecords,

    #[error("Missing lookup records for full scoring")]
    MissingLookupRecords,

    #[error("alignment {j} still has reads after parallel processing")]
    AlignmentStillHasReadsAfterParallelProcessing { j: usize },

    #[error("Score error stream {0}: {1}")]
    ScoreStreamError(usize, String),

    #[error("Scorer worker exited unexpectedly")]
    ScorerWorkerExited,

    #[error("Work channel closed unexpectedly")]
    WorkChannelClosed,

    #[error("All scorer workers exited unexpectedly")]
    AllScorerWorkersExited,

    #[error("NW error: {0}")]
    NeedlemanWunschError(String),

    #[error("Scorer workers exited unexpectedly")]
    ScorerWorkersExited,

    #[error("Error scoring fragment for alignment {aln_idx}: {message}\n{state}")]
    FragmentScoringError {
        aln_idx: usize,
        message: String,
        state: String,
    },

    #[error("Scoring error: {message}\n{state}")]
    ScoringError { message: String, state: String },

    #[error("BUG: unmapped record should already have been excluded")]
    UnmappedRecordInMdCigFlags,

    #[error("--chimeric-pairs: expected format 'A:B' (e.g. '0:1'), got '{raw}'")]
    ChimericPairsInvalidFormat { raw: String },

    #[error("--chimeric-pairs: '{index_str}' is not a valid stream index")]
    ChimericPairsInvalidIndex { index_str: String },

    #[error("--chimeric-pairs: stream index must differ (got '{raw}')")]
    ChimericPairsIdenticalIndices { raw: String },

    #[error("--chimeric-pairs is only supported with --matching-algorithm namesorted")]
    ChimericPairsRequiresNamesorted,

    #[error(
        "--ambiguous-regions / --diagnostic-variants have no effect with --matching-algorithm namesorted"
    )]
    NamesortedUnsupportedOptions,

    #[error("Multi-threaded scoring is only supported with --matching-algorithm namesorted.")]
    MultiThreadedScoringRequiresNamesorted,

    #[error(
        "Single alignment stream detected. If this is intentional for within-species disambiguation, \
        please pass --single-alignment-mode."
    )]
    SingleStreamMissingFlag,

    #[error(
        "Cannot use single alignment mode with stdin because the stream must be duplicated via file system access."
    )]
    SingleStreamStdinUnsupported,

    #[error("--single-alignment-mode is only supported with --matching-algorithm namesorted.")]
    SingleStreamRequiresNamesorted,

    #[error("at most one stream may be read from stdin")]
    MultipleStdinStreams,

    #[error("stdin is only supported with --matching-algorithm namesorted")]
    StdinRequiresNamesorted,

    #[error(
        "--single-alignment-mode can only be used with exactly 1 alignment stream (got {count})."
    )]
    SingleAlignmentModeExpectsOneStream { count: usize },

    #[error(
        "At least two alignments required when not running in single alignment mode (got {count})."
    )]
    InsufficientAlignmentStreams { count: usize },

    #[error("BUG: expected Scoring on both sides")]
    BugExpectedScoringOnBothSides,

    #[error(
        "More than 2 alignment streams requires --matching-algorithm namesorted \
        (hashlookup and collated are limited to 2 streams by design)"
    )]
    MultiStreamRequiresNamesorted,

    #[error("namesorted supports at most {max} alignment streams (got {count})")]
    MaxStreamsExceeded { count: usize, max: usize },

    #[error(
        "Cannot use --merged-output in combination with --output, --discarded-output, or --ambiguous-output."
    )]
    MergedOutputConflict,

    #[error("--ambiguous-regions: at most 2 files (one per stream), got {count}")]
    TooManyAmbiguousRegionsFiles { count: usize },

    #[error("--diagnostic-variants: at most 2 files (one per stream), got {count}")]
    TooManyDiagnosticVariantsFiles { count: usize },

    #[error("Sample variant stream index {idx} out of bounds (max {max})")]
    SampleVariantIndexOutOfBounds { idx: usize, max: usize },

    #[error("Population variant stream index {idx} out of bounds (max {max})")]
    PopulationVariantIndexOutOfBounds { idx: usize, max: usize },

    #[error(
        "Single alignment mode requires both strain slots (index 0 and 1) to have a variant profile. \
        An option like '--sample-variant 0:a.vcf --population-variant 0:b.vcf' is invalid because strain 1 \
        then has no variants."
    )]
    SingleStreamMissingVariantProfiles,

    #[error(
        "More output paths ({count}) than logical alignment processing streams ({max}) specified"
    )]
    TooManyOutputPaths { count: usize, max: usize },

    #[error(
        "More discarded output paths ({count}) than logical alignment processing streams ({max}) specified"
    )]
    TooManyDiscardedOutputPaths { count: usize, max: usize },

    #[error(
        "Only one ambiguous output file is allowed when operating on a single alignment stream (got {count})."
    )]
    SingleStreamTooManyAmbiguousOutputs { count: usize },

    #[error("More ambiguous output paths ({count}) than input ({max}) specified")]
    TooManyAmbiguousOutputPaths { count: usize, max: usize },
}
