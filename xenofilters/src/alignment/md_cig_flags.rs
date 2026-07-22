use crate::alignment::tie_break_bool;
use crate::Error;
use noodles::sam::alignment::record::{data::field::Tag, data::field::Value, Cigar, Flags, Record};
use std::cmp::Ordering;

pub struct MdCigFlags<'r> {
    flags: &'r Flags,
    md: &'r [u8],
    cig: Box<dyn Cigar + 'r>,
    /// Count of supplementary alignments declared in the SA:Z tag.
    ///
    /// SA:Z: format: "rname,pos,strand,CIGAR,mapQ,NM;" per supplementary,
    /// so semicolon count == supplementary count. A non-zero value means
    /// supplementary records will follow in the BAM stream for this read.
    /// We only need this on the PRIMARY record: once all primaries have been
    /// seen and none carries a non-empty SA tag, no supplementaries can follow.
    supp_count: usize,
    /// Read sequence (ASCII bases). `None` when bisulfite mode is off or
    /// the sequence is unavailable (hard-clipped supplementary).
    seq: Option<Vec<u8>>,
    /// Whether this record is on the reverse strand. Used for strand-specific
    /// bisulfite conversion detection (C->T forward, G->A reverse).
    is_reverse: bool,
}

impl<'r> MdCigFlags<'r> {
    /// Build an `MdCigRef` from a stored `MdCigFlags` and its matching record.
    pub fn try_from_record<R: Record>(
        record: &'r R,
        flags: &'r Flags,
        bisulfite: bool,
    ) -> Result<Self, Error> {
        if flags.is_unmapped() {
            return Err(Error::UnmappedRecordInMdCigFlags);
        }

        let md = match record
            .data()
            .get(&Tag::MISMATCHED_POSITIONS)
            .transpose()?
            .ok_or(Error::MissingMdTag)?
        {
            Value::String(bstr) => {
                // SAFETY: 'r is the lifetime of `record`, and `bstr` is derived from that borrow.
                let slice: &[u8] = bstr.as_ref();
                let md: &'r [u8] =
                    unsafe { std::slice::from_raw_parts(slice.as_ptr(), slice.len()) };
                md
            }
            _ => return Err(Error::UnexpectedMdTagValueType),
        };

        // SA:Z: "rname,pos,strand,CIGAR,mapQ,NM;" -- one ';' per supplementary.
        // Only meaningful on a primary record; supplementary records' own SA tags
        // list siblings, not additional supplementaries, so double-counting is
        // harmless (the `flags.is_supplementary()` guard in cmp_perfect already
        // blocks fast-pathing for supplementary records themselves).
        let supp_count = match record
            .data()
            .get(&Tag::OTHER_ALIGNMENTS) // SA
            .transpose()?
        {
            Some(Value::String(s)) => {
                let b: &[u8] = s.as_ref();
                b.iter().filter(|&&c| c == b';').count()
            }
            _ => 0,
        };

        // Sequence for bisulfite conversion detection.
        // Only extracted when bisulfite mode is on; zero cost otherwise.
        let seq: Option<Vec<u8>> = if bisulfite {
            let s = record.sequence();
            // Hard-clipped reads have an empty sequence in BAM; treat as None.
            if s.is_empty() {
                None
            } else {
                // Clone the sequence into a Vec for storage.
                if s.is_empty() {
                    None
                } else {
                    Some(s.iter().collect())
                }
            }
        } else {
            None
        };

        let cig: Box<dyn Cigar + 'r> = record.cigar();
        Ok(MdCigFlags {
            flags,
            md,
            cig,
            supp_count,
            seq,
            is_reverse: flags.is_reverse_complemented(),
        })
    }

    /// Returns `true` only when ALL of the following hold:
    ///  - `SA:Z` tag absent (no supplementary alignments pending in the stream).
    ///  - Exactly one CIGAR operation (read aligns end-to-end with no gaps or clips).
    ///  - MD tag is all digits (zero mismatches against the reference).
    ///
    /// Any supplementary record that would follow in the BAM stream carries a
    /// chimeric-junction penalty; a fragment with pending supplementaries therefore
    /// cannot be fast-pathed as perfect.
    ///
    /// In bisulfite mode, the perfect fast path is suppressed entirely (returns false)
    /// to ensure that all conversions are scored correctly.
    pub(crate) fn is_perfect(&self) -> bool {
        // Bisulfite mode: suppress perfect fast path for safety
        if self.seq.is_some() {
            return false;
        }
        self.supp_count == 0 && self.cig.len() == 1 && self.md.iter().all(|&b| b.is_ascii_digit())
    }

    pub(crate) fn supp_count(&self) -> usize {
        self.supp_count
    }
    pub(crate) fn is_supplementary(&self) -> bool {
        self.flags.is_supplementary()
    }
    pub(super) fn is_reverse_complemented(&self) -> bool {
        self.flags.is_reverse_complemented()
    }
    pub(super) fn is_last_segment(&self) -> bool {
        self.flags.is_last_segment()
    }
    pub(crate) fn get_md(&self) -> &[u8] {
        self.md
    }
    pub(crate) fn get_cigar(&self) -> &(dyn Cigar + 'r) {
        &self.cig
    }
    pub(crate) fn seq(&self) -> Option<&[u8]> {
        self.seq.as_deref()
    }
}

impl<'r> PartialOrd for MdCigFlags<'r> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        tie_break_bool(self.is_perfect(), other.is_perfect())
    }
}

impl<'r> PartialEq for MdCigFlags<'r> {
    fn eq(&self, other: &Self) -> bool {
        self.partial_cmp(other) == Some(Ordering::Equal)
    }
}

#[cfg(test)]
mod tests;
