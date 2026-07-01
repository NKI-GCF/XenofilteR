//!
//! Design principle: always accumulate `RecordKind`s during pass 1, with the
//! cheapest representation that fits the record's role:
//!
//! - `Secondary`        — virtual offset + flags only; never scored, tags along.
//! - `UnmappedPrimary`   — virtual offset + flags + qualities; never NW-scored
//!                          (no CIGAR/MD exists), qualities reserved for future
//!                          quality-weighted ambiguous tie-breaking.
//! - `Mapped`            — boxed full record (CIGAR/MD/qualities/etc.) for
//!                          primary or supplementary mapped alignments; the
//!                          only variant ever passed to NW scoring.
//!
//! Boxing the expensive variant keeps `RecordKind` small (pointer + tag) so
//! `SmallVec<[RecordKind; 2]>` stays cheap even when most fragments are
//! perfect/unmapped and never reach the boxed branch.
//!
//! Early-assignability is evaluated at fragment *completion* time (see
//! `StreamAccumulator::classify`). When neither fast path applies, a
//! per-mate `MateKind` summary is computed alongside the records so that
//! `resolve_fragment` can skip NW scoring for individual mates that are
//! identical-contribution in both streams (Unmapped-Unmapped or
//! Perfect-Perfect), even when the fragment as a whole requires scoring.

use crate::region::{AmbiguousRegions, DiagnosticVariants};
use noodles::sam::alignment::record::Flags;
use smallvec::SmallVec;

// ---------------------------------------------------------------------------
// MappedRecord — full data, only ever built for mapped (primary or
// supplementary) records.
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub(crate) struct MappedRecord {
    pub(crate) flags: Flags,
    pub(crate) ref_id: usize,
    /// 1-based alignment start (BAM convention).
    pub(crate) pos: usize,
    /// Reference span derived from CIGAR.
    pub(crate) ref_len: usize,
    /// Raw BAM-encoded CIGAR bytes (4 bytes per op, little-endian u32).
    pub(crate) cigar_bytes: Vec<u8>,
    pub(crate) md: Vec<u8>,
    pub(crate) qualities: Vec<u8>,
    pub(crate) virtual_offset: u64,
    /// Supplementary alignment count from the SA:Z tag (semicolons counted).
    pub(crate) supp_count: usize,
}

impl MappedRecord {
    fn is_perfect(&self) -> bool {
        if self.supp_count > 0 {
            return false;
        }
        let cigar_ok = self.cigar_bytes.len() == 4 && (self.cigar_bytes[0] & 0x0F) == 0;
        let md_ok = !self.md.is_empty() && self.md.iter().all(|&b| b.is_ascii_digit());
        cigar_ok && md_ok
    }
}

// ---------------------------------------------------------------------------
// RecordKind — minimal-cost representation per record role.
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub(crate) enum RecordKind {
    /// Secondary alignment: never scored, tags along on the fragment's
    /// final decision. Only a virtual offset is needed for pass-2 retrieval.
    Secondary { flags: Flags, virtual_offset: u64 },
    /// Unmapped primary: never NW-scored (no CIGAR/MD exists for an unmapped
    /// read). Qualities are retained — not currently consumed by scoring, but
    /// reserved for a possible future quality-weighted ambiguous tie-break.
    UnmappedPrimary {
        flags: Flags,
        virtual_offset: u64,
        #[allow(dead_code)]
        qualities: Vec<u8>,
    },
    /// Mapped primary or supplementary record. The only variant ever fed
    /// into NW scoring or perfection checks.
    Mapped(Box<MappedRecord>),
}

impl RecordKind {
    pub(crate) fn flags(&self) -> Flags {
        match self {
            RecordKind::Secondary { flags, .. } => *flags,
            RecordKind::UnmappedPrimary { flags, .. } => *flags,
            RecordKind::Mapped(m) => m.flags,
        }
    }
    pub(crate) fn virtual_offset(&self) -> u64 {
        match self {
            RecordKind::Secondary { virtual_offset, .. } => *virtual_offset,
            RecordKind::UnmappedPrimary { virtual_offset, .. } => *virtual_offset,
            RecordKind::Mapped(m) => m.virtual_offset,
        }
    }
    fn is_primary(&self) -> bool {
        let f = self.flags();
        !f.is_secondary() && !f.is_supplementary()
    }
    pub(crate) fn is_supplementary(&self) -> bool {
        self.flags().is_supplementary()
    }
}

/// Encode the segment role of a primary record as a byte:
/// `0x40` = first segment (read 1), `0x80` = last segment (read 2),
/// `0x00` = single-end or unknown. Shared mate-slot encoding with the
/// chimeric-detection module.
pub(crate) fn segment_id(flags: &Flags) -> u8 {
    match (flags.is_first_segment(), flags.is_last_segment()) {
        (true, _) => 0x40,
        (false, true) => 0x80,
        (false, false) => 0x00,
    }
}

/// Map a segment id to a fixed mate slot (0 = forward/single-end, 1 = reverse).
pub(crate) fn mate_slot(seg_id: u8) -> usize {
    if seg_id == 0x80 {
        1
    } else {
        0
    }
}

// ---------------------------------------------------------------------------
// MateKind — per-mate classification used for partial-scoring cancellation
// ---------------------------------------------------------------------------

/// Per-mate alignment quality classification, used to detect when a mate's
/// contribution is provably identical across both streams and can be
/// excluded from NW scoring entirely.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MateKind {
    /// Mate is unmapped — never NW-scored in either stream.
    Unmapped,
    /// Mate is a perfect match — scores `log_lik_match[q_i]` at every
    /// position, identical in both streams (same quality string).
    Perfect,
    /// Mate is mapped but imperfect — must be scored normally.
    Other,
}

// ---------------------------------------------------------------------------
// StreamAccumulator — accumulates records for one stream before classification
// ---------------------------------------------------------------------------

#[derive(Default)]
struct StreamAccumulator {
    records: SmallVec<[RecordKind; 2]>,
    primary_count: usize,
}

impl StreamAccumulator {
    fn push(&mut self, rec: RecordKind) {
        if rec.is_primary() {
            self.primary_count += 1;
        }
        self.records.push(rec);
    }

    fn classify(
        self,
        bed: Option<&AmbiguousRegions>,
        vcf: Option<&DiagnosticVariants>,
    ) -> StreamKind {
        if self.primary_count == 0 {
            return StreamKind::Scoring {
                records: Box::new(self.records),
                mate_kinds: [None, None],
            };
        }

        let primaries: SmallVec<[&RecordKind; 2]> =
            self.records.iter().filter(|r| r.is_primary()).collect();

        // Fast-path 1: all primaries unmapped.
        let all_unmapped = primaries
            .iter()
            .all(|r| matches!(r, RecordKind::UnmappedPrimary { .. }));
        if all_unmapped {
            let offsets = self.records.iter().map(|r| r.virtual_offset()).collect();
            return StreamKind::Early {
                kind: EarlyKind::AllUnmapped,
                virtual_offsets: offsets,
            };
        }

        // Fast-path 2: all primaries perfect, no region overlap.
        let all_perfect = primaries.iter().all(|r| match r {
            RecordKind::Mapped(m) => {
                m.is_perfect()
                    && !bed.is_some_and(|b| b.overlaps(m.ref_id, m.pos, m.pos + m.ref_len))
                    && !vcf.is_some_and(|v| v.overlaps(m.ref_id, m.pos, m.pos + m.ref_len))
            }
            _ => false,
        });
        if all_perfect {
            let offsets = self.records.iter().map(|r| r.virtual_offset()).collect();
            return StreamKind::Early {
                kind: EarlyKind::AllPerfect,
                virtual_offsets: offsets,
            };
        }

        // Neither fast path: build per-mate classification for partial
        // scoring cancellation in `resolve_fragment`.
        let mut mate_kinds: [Option<MateKind>; 2] = [None, None];
        for r in &primaries {
            let (seg_id, kind) = match r {
                RecordKind::UnmappedPrimary { flags, .. } => {
                    (segment_id(flags), MateKind::Unmapped)
                }
                RecordKind::Mapped(m) => {
                    let perfect_here = m.is_perfect()
                        && !bed.is_some_and(|b| b.overlaps(m.ref_id, m.pos, m.pos + m.ref_len))
                        && !vcf.is_some_and(|v| v.overlaps(m.ref_id, m.pos, m.pos + m.ref_len));
                    (
                        segment_id(&m.flags),
                        if perfect_here {
                            MateKind::Perfect
                        } else {
                            MateKind::Other
                        },
                    )
                }
                RecordKind::Secondary { .. } => unreachable!("filtered to primaries"),
            };
            mate_kinds[mate_slot(seg_id)] = Some(kind);
        }
        drop(primaries);

        StreamKind::Scoring {
            records: Box::new(self.records),
            mate_kinds,
        }
    }
}

// ---------------------------------------------------------------------------
// StreamKind — post-classification state
// ---------------------------------------------------------------------------

/// Why a stream was fast-pathed without per-base NW scoring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EarlyKind {
    AllUnmapped,
    AllPerfect,
}

#[derive(Default)]
pub(crate) enum StreamKind {
    Early {
        virtual_offsets: SmallVec<[u64; 2]>,
        kind: EarlyKind,
    },
    Scoring {
        records: Box<SmallVec<[RecordKind; 2]>>,
        /// Per-mate classification: index 0 = forward/single-end, index 1 = reverse.
        /// `None` when that mate slot is absent in this stream.
        mate_kinds: [Option<MateKind>; 2],
    },
    #[default]
    Empty,
}

impl StreamKind {
    pub(crate) fn is_empty(&self) -> bool {
        matches!(self, StreamKind::Empty)
    }
    pub(crate) fn early_kind(&self) -> Option<EarlyKind> {
        match self {
            StreamKind::Early { kind, .. } => Some(*kind),
            _ => None,
        }
    }
    pub(crate) fn virtual_offsets(&self) -> SmallVec<[u64; 2]> {
        match self {
            StreamKind::Early {
                virtual_offsets, ..
            } => virtual_offsets.clone(),
            StreamKind::Scoring { records, .. } => {
                records.iter().map(|r| r.virtual_offset()).collect()
            }
            StreamKind::Empty => SmallVec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// PendingFragment — unchanged structurally; type updated to RecordKind
// ---------------------------------------------------------------------------

pub(crate) struct PendingFragment {
    driving_buf: StreamAccumulator,
    lookup_buf: StreamAccumulator,
    pub(crate) driving: StreamKind,
    pub(crate) lookup: StreamKind,
    pub(crate) supplementary_offsets: [SmallVec<[u64; 1]>; 2],
    pub(crate) seq_nr: u64,
    pub(crate) is_paired: Option<bool>,
}

impl PendingFragment {
    pub(crate) fn new(seq_nr: u64) -> Self {
        Self {
            driving_buf: StreamAccumulator::default(),
            lookup_buf: StreamAccumulator::default(),
            driving: StreamKind::Empty,
            lookup: StreamKind::Empty,
            supplementary_offsets: [SmallVec::new(), SmallVec::new()],
            seq_nr,
            is_paired: None,
        }
    }

    pub(crate) fn push(
        &mut self,
        rec: RecordKind,
        nr: usize,
        bed: Option<&AmbiguousRegions>,
        vcf: Option<&DiagnosticVariants>,
    ) -> bool {
        if self.is_paired.is_none() {
            self.is_paired = Some(rec.flags().is_segmented());
        }
        if rec.is_supplementary() {
            self.supplementary_offsets[nr].push(rec.virtual_offset());
            return self.check_complete(bed, vcf);
        }
        if nr == 0 {
            self.driving_buf.push(rec);
        } else {
            self.lookup_buf.push(rec);
        }
        self.check_complete(bed, vcf)
    }

    fn expected_primaries(&self) -> usize {
        self.is_paired.map_or(1, |p| p as usize + 1)
    }

    fn check_complete(
        &mut self,
        bed: Option<&AmbiguousRegions>,
        vcf: Option<&DiagnosticVariants>,
    ) -> bool {
        let exp = self.expected_primaries();
        if self.driving.is_empty() && self.driving_buf.primary_count >= exp {
            let buf = std::mem::take(&mut self.driving_buf);
            self.driving = buf.classify(bed, vcf);
        }
        if self.lookup.is_empty() && self.lookup_buf.primary_count >= exp {
            let buf = std::mem::take(&mut self.lookup_buf);
            self.lookup = buf.classify(bed, vcf);
        }
        !self.driving.is_empty() && !self.lookup.is_empty()
    }
}

// ---------------------------------------------------------------------------
// FragmentTable
// ---------------------------------------------------------------------------

pub(crate) type FragmentTable = std::collections::HashMap<Box<[u8]>, PendingFragment>;

pub(crate) fn insert(
    table: &mut FragmentTable,
    rec: RecordKind,
    canonical_name: Box<[u8]>,
    nr: usize,
    seq_counter: &mut u64,
    bed: Option<&AmbiguousRegions>,
    vcf: Option<&DiagnosticVariants>,
) -> (Box<[u8]>, bool) {
    let entry = table.entry(canonical_name.clone()).or_insert_with(|| {
        let sn = *seq_counter;
        *seq_counter += 1;
        PendingFragment::new(sn)
    });
    let complete = entry.push(rec, nr, bed, vcf);
    (canonical_name, complete)
}
