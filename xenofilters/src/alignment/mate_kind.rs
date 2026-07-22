//! Shared per-mate alignment-quality classification used by all three
//! backends to detect mates whose NW scoring contribution is provably
//! identical across streams (same physical read, same quality string),
//! and can therefore be excluded from per-base scoring entirely.

use noodles::sam::alignment::record::Flags;

/// Per-mate alignment quality classification, used to detect when a mate's
/// contribution is provably identical across both streams and can be
/// excluded from NW scoring entirely.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MateKind {
    /// Mate unmapped: contributes 0 to NW scoring in this stream.
    Unmapped,
    /// Mate is a perfect match (no SA, single CIGAR op, all-digit MD):
    /// contributes `Sigma log_lik_match[q_i]`, identical across streams because
    /// the quality string is the same physical read regardless of reference.
    Perfect,
    /// Mapped but imperfect, malformed, or disqualified by an overlapping
    /// ambiguous region / diagnostic variant: must be scored normally.
    Other,
}

/// Encode the segment role of a record: `0x40` = read 1, `0x80` = read 2,
/// `0x00` = single-end.
pub(crate) fn segment_id(flags: &Flags) -> u8 {
    match (flags.is_first_segment(), flags.is_last_segment()) {
        (true, _) => 0x40,
        (false, true) => 0x80,
        (false, false) => 0x00,
    }
}

/// Map a segment id to a fixed mate slot: 0 = forward/single-end, 1 = reverse.
pub(crate) fn mate_slot(seg_id: u8) -> usize {
    if seg_id == 0x80 {
        1
    } else {
        0
    }
}

/// Implementors can classify their own per-mate alignment quality without
/// performing full per-base NW scoring.
///
/// - **HashLookup / Collated** (2 streams): a mate cancels when both
///   streams agree (`Unmapped`/`Unmapped` or `Perfect`/`Perfect`).
/// - **LineByLine** (N <= `MAX_STREAMS` streams): a mate cancels only when
///   *every* stream still competing agrees on the same non-`Other` kind.
pub(crate) trait MateClassifiable {
    /// Index 0 = forward/single-end mate, index 1 = reverse mate.
    /// `None` means that mate slot is absent for this fragment in this stream.
    fn mate_kinds(&self) -> [Option<MateKind>; 2];
}
