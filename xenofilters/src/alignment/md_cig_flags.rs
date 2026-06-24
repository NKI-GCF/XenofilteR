use crate::alignment::{BaseOp, ScoreOpIter};
use anyhow::{anyhow, ensure, Result};
use noodles::sam::alignment::record::{data::field::Tag, data::field::Value, Cigar, Flags, Record};
use std::cmp::Ordering;

pub(crate) struct MdCigFlags<'r> {
    flags: &'r Flags,
    md: &'r [u8],
    cig: Box<dyn Cigar + 'r>,
}

impl<'r> MdCigFlags<'r> {
    /// Build an `MdCigRef` from a stored `MdCigFlags` and its matching record.
    pub(crate) fn try_from_record<R: Record>(record: &'r R, flags: &'r Flags) -> Result<Self> {
        ensure!(
            !flags.is_unmapped(),
            "BUG: unmapped record should already have been excluded"
        );
        match record
            .data()
            .get(&Tag::MISMATCHED_POSITIONS)
            .transpose()?
            .ok_or_else(|| anyhow!("missing MD tag"))?
        {
            Value::String(bstr) => {
                // SAFETY: 'r is the lifetime of `record`, and `bstr` is derived from that borrow.
                let slice: &[u8] = bstr.as_ref();
                let md: &'r [u8] =
                    unsafe { std::slice::from_raw_parts(slice.as_ptr(), slice.len()) };

                let cig: Box<dyn Cigar + 'r> = record.cigar();
                Ok(MdCigFlags { flags, md, cig })
            }
            _ => Err(anyhow!("unexpected MD tag value type")),
        }
    }

    pub(crate) fn is_perfect(&self) -> bool {
        // Single cigar operation and MD string is all digits (no mismatches).
        self.cig.len() == 1 && self.md.iter().all(|&b| b.is_ascii_digit())
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
