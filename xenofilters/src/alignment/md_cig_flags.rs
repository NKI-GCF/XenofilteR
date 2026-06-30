use crate::Error;
use noodles::sam::alignment::record::{Cigar, Flags, Record, data::field::Tag, data::field::Value};
use std::cmp::Ordering;

pub(crate) struct MdCigFlags<'r> {
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
}

impl<'r> MdCigFlags<'r> {
    /// Build an `MdCigRef` from a stored `MdCigFlags` and its matching record.
    pub(crate) fn try_from_record<R: Record>(
        record: &'r R,
        flags: &'r Flags,
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

        // SA:Z: "rname,pos,strand,CIGAR,mapQ,NM;" — one ';' per supplementary.
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

        let cig: Box<dyn Cigar + 'r> = record.cigar();
        Ok(MdCigFlags {
            flags,
            md,
            cig,
            supp_count,
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
    pub(crate) fn is_perfect(&self) -> bool {
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
}

impl<'r> PartialOrd for MdCigFlags<'r> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self.is_perfect(), other.is_perfect()) {
            (true, true) => Some(Ordering::Equal), // both unmapped => tie-break with next pair
            (true, false) => Some(Ordering::Less), // self worse
            (false, true) => Some(Ordering::Greater), // other worse
            (false, false) => None, // Slow path, first cig/md after init, then per base evaluation.
        }
    }
}

impl<'r> PartialEq for MdCigFlags<'r> {
    fn eq(&self, other: &Self) -> bool {
        self.partial_cmp(other) == Some(Ordering::Equal)
    }
}

#[cfg(test)]
mod tests;
