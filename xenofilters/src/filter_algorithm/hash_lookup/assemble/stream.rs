use crate::Error;
use crate::alignment::{MateKind, mate_slot, segment_id};
use crate::region::tabix_query::{TabixBed, TabixVcf};
use noodles::sam::alignment::record::Flags;
use smallvec::SmallVec;

// ---------------------------------------------------------------------------
// MappedRecord -- full data, only ever built for mapped (primary or
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
    pub(crate) sequence: Vec<u8>,
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
// RecordKind -- minimal-cost representation per record role.
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub(crate) enum RecordKind {
    /// Secondary alignment: never scored, tags along on the fragment's
    /// final decision. Only a virtual offset is needed for pass-2 retrieval.
    Secondary { flags: Flags, virtual_offset: u64 },
    /// Unmapped primary: never NW-scored (no CIGAR/MD exists for an unmapped
    /// read). Qualities are retained -- not currently consumed by scoring, but
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

macro_rules! extract_record_field {
    ($self:expr, $field:ident) => {
        match $self {
            RecordKind::Secondary { $field, .. } => *$field,
            RecordKind::UnmappedPrimary { $field, .. } => *$field,
            RecordKind::Mapped(m) => m.$field,
        }
    };
}

impl RecordKind {
    pub(crate) fn flags(&self) -> Flags {
        extract_record_field!(self, flags)
    }
    pub(crate) fn virtual_offset(&self) -> u64 {
        extract_record_field!(self, virtual_offset)
    }
    fn is_primary(&self) -> bool {
        let f = self.flags();
        !f.is_secondary() && !f.is_supplementary()
    }
    pub(crate) fn is_supplementary(&self) -> bool {
        self.flags().is_supplementary()
    }
}

// ---------------------------------------------------------------------------
// StreamAccumulator -- accumulates records for one stream before classification
// ---------------------------------------------------------------------------

#[derive(Default)]
pub(super) struct StreamAccumulator {
    records: SmallVec<[RecordKind; 2]>,
    pub(super) primary_count: usize,
}

impl StreamAccumulator {
    pub(super) fn push(&mut self, rec: RecordKind) {
        if rec.is_primary() {
            self.primary_count += 1;
        }
        self.records.push(rec);
    }

    pub(super) fn classify(
        self,
        bed: Option<&TabixBed>,
        vcf: Option<&TabixVcf>,
    ) -> Result<StreamKind, Error> {
        if self.primary_count == 0 {
            return Ok(StreamKind::Scoring {
                records: Box::new(self.records),
                mate_kinds: [None, None],
            });
        }

        let primaries: SmallVec<[&RecordKind; 2]> =
            self.records.iter().filter(|r| r.is_primary()).collect();

        // Fast-path 1: all primaries unmapped.
        let all_unmapped = primaries
            .iter()
            .all(|r| matches!(r, RecordKind::UnmappedPrimary { .. }));
        if all_unmapped {
            let offsets = self.records.iter().map(|r| r.virtual_offset()).collect();
            return Ok(StreamKind::Early {
                kind: EarlyKind::AllUnmapped,
                virtual_offsets: offsets,
            });
        }

        // Fast-path 2: all primaries perfect, no region overlap.
        let mut all_perfect = true;
        for r in &primaries {
            match r {
                RecordKind::Mapped(m) => {
                    if !m.is_perfect()
                        || bed
                            .map(|b| b.overlaps(m.ref_id, m.pos, m.pos + m.ref_len))
                            .transpose()?
                            != Some(true)
                        || vcf
                            .map(|v| v.overlaps(m.ref_id, m.pos, m.pos + m.ref_len))
                            .transpose()?
                            != Some(true)
                    {
                        all_perfect = false;
                        break;
                    }
                }
                _ => {
                    all_perfect = false;
                    break;
                }
            }
        }
        if all_perfect {
            let offsets = self.records.iter().map(|r| r.virtual_offset()).collect();
            return Ok(StreamKind::Early {
                kind: EarlyKind::AllPerfect,
                virtual_offsets: offsets,
            });
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
                        && bed
                            .map(|b| b.overlaps(m.ref_id, m.pos, m.pos + m.ref_len))
                            .transpose()?
                            != Some(true)
                        && vcf
                            .map(|v| v.overlaps(m.ref_id, m.pos, m.pos + m.ref_len))
                            .transpose()?
                            != Some(true);
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

        Ok(StreamKind::Scoring {
            records: Box::new(self.records),
            mate_kinds,
        })
    }
}

// ---------------------------------------------------------------------------
// StreamKind -- post-classification state
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

#[cfg(test)]
mod is_perfect_tests {
    use super::*;
    use noodles::sam::alignment::record::Flags;

    fn encode_op(len: u32, code: u32) -> [u8; 4] {
        ((len << 4) | code).to_le_bytes()
    }

    fn base(cigar_bytes: Vec<u8>, md: &[u8], supp_count: usize) -> MappedRecord {
        MappedRecord {
            flags: Flags::empty(),
            ref_id: 0,
            pos: 1,
            ref_len: 10,
            cigar_bytes,
            md: md.to_vec(),
            qualities: vec![30; 10],
            sequence: vec![b'A'; 10],
            virtual_offset: 0,
            supp_count,
        }
    }

    #[test]
    fn single_match_op_clean_md_is_perfect() {
        let m = base(encode_op(10, 0 /* Match */).to_vec(), b"10", 0);
        assert!(m.is_perfect());
    }

    #[test]
    fn pending_supplementary_forces_imperfect() {
        let m = base(encode_op(10, 0).to_vec(), b"10", 1);
        assert!(
            !m.is_perfect(),
            "supp_count > 0 must disqualify perfect fast-path"
        );
    }

    #[test]
    fn multi_op_cigar_forces_imperfect() {
        // 8 bytes = 2 ops (e.g. 5M5S) -- must not be misread as one op.
        let mut bytes = encode_op(5, 0).to_vec();
        bytes.extend_from_slice(&encode_op(5, 4 /* SoftClip */));
        let m = base(bytes, b"5", 0);
        assert!(!m.is_perfect());
    }

    #[test]
    fn single_non_match_op_forces_imperfect() {
        // e.g. a lone Insertion op (code 1) -- same length encoding as Match
        // but wrong op kind; must not be treated as perfect.
        let m = base(encode_op(10, 1).to_vec(), b"10", 0);
        assert!(!m.is_perfect());
    }

    #[test]
    fn empty_md_forces_imperfect() {
        let m = base(encode_op(10, 0).to_vec(), b"", 0);
        assert!(!m.is_perfect());
    }

    #[test]
    fn md_with_mismatch_letter_forces_imperfect() {
        let m = base(encode_op(10, 0).to_vec(), b"5A4", 0);
        assert!(!m.is_perfect());
    }
}
