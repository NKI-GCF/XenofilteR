use noodles::sam::alignment::Record;
use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;
use crate::aln_stream::AlignmentStream;
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::record::Flags;
use anyhow::{Result, anyhow};

pub(crate) struct MdCigFlags<'r> {
    flags: &'r Flags,
    md: &'r [u8],
    cig: Box<dyn Cigar + 'r>,
}

impl<'r> MdCigFlags<'r> {
    /// Build an `MdCigRef` from a stored `MdCigFlags` and its matching record.
    /// Returns `None` when the read is unmapped or secondary (no MD/CIG needed).
    pub(crate) fn try_from_record<R: Record>(
        flags: &'r Flags,
        record: &'r R,
    ) -> Result<Self> {
        match record.data().get(&Tag::MISMATCHED_POSITIONS).transpose()?
            .ok_or_else(|| anyhow!("missing MD tag"))? {
            Value::String(bstr) => {
                // SAFETY: 'r is the lifetime of `record`, and `bstr` is derived from that borrow.
                let slice: &[u8] = bstr.as_ref();
                let md: &'r [u8] = unsafe { std::slice::from_raw_parts(slice.as_ptr(), slice.len()) };

                let cig: Box<dyn Cigar + 'r> = record.cigar();
                Ok(MdCigFlags { flags, md, cig })
            }
            _ => Err(anyhow!("unexpected MD tag value type")),
        }
    }

    fn is_perfect(&self) -> bool {
        // Single cigar operation and MD string is all digits (no mismatches).
        self.cig.len() == 1 && self.md.iter().all(|&b| b.is_ascii_digit())
    }
    pub(super) fn is_reverse_complemented(&self) -> bool {
        self.flags.is_reverse_complemented()
    }
    pub(super) fn is_last_segment(&self) -> bool {
        self.flags.is_last_segment()
    }
    pub(super) fn get_md(&self) -> &[u8] {
        self.md
    }
    pub(super) fn get_cigar(&self) -> &(dyn Cigar + 'r) {
        &self.cig
    }
}

impl<'r> PartialOrd for MdCigFlags<'r> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self.is_perfect(), other.is_perfect()) {
            (true, true) => Some(Ordering::Equal), // both unmapped => tie-break with next pair
            (true, false) => Some(Ordering::Less), // self worse
            (false, true) => Some(Ordering::Greater),    // other worse
            (false, false) => None, // Slow path, first cig/md after init, then per base evaluation.
        }
    }
}

impl<'r> PartialEq for MdCigFlags<'r> {
    fn eq(&self, other: &Self) -> bool {
        self.partial_cmp(other) == Some(Ordering::Equal)
    }
}

#[derive(PartialEq, Debug)]
pub(crate) struct FragmentState<R> {
    flags: SmallVec<[Flags; 8]>,
    records: SmallVec<[R; 8]>,
    species_nr: usize,
}

impl<R: Record> FragmentState<R> {
    pub(crate) fn from_record(r: R, species_nr: usize) -> Result<Self> {
        Ok(FragmentState {
            flags: smallvec![r.flags()?],
            records: smallvec![r],
            species_nr,
        })
    }
    pub(crate) fn add_record(&mut self, r: R) -> Result<()> {
        self.flags.push(r.flags()?);
        self.records.push(r);
        Ok(())
    }
    pub(crate) fn md_cig_refs(&self) -> Result<SmallVec<[MdCigFlags<'_>; 8]>> {
        self.flags
            .iter()
            .zip(self.records.iter())
            .map(|(flags, rec)| MdCigFlags::try_from_record(flags, rec))
            .collect()
    }
    pub(crate) fn get_records(&self) -> &[R] {
        &self.records
    }
    pub(crate) fn drain_records(&mut self) -> impl Iterator<Item = R> + '_ {
        self.records.drain(..)
    }
    #[must_use]
    pub(crate) fn first_qname(&self) -> &[u8] {
        self.records.first().and_then(|r| r.name()).map_or(b"", |n| n.as_ref())
    }
    #[must_use]
    pub(crate) fn get_nr(&self) -> usize {
        self.species_nr
    }
    pub(crate) fn flags(&self, i: usize) -> Option<&Flags> {
        self.flags.get(i)
    }
    // Rationale: if all unmapped, the first is unmapped and its mate also, when paired-end
    fn is_all_unmapped(&self) -> bool {
        let f = &self.flags[0];
        f.is_unmapped() && (!f.is_segmented() || f.is_mate_unmapped())
    }

    pub(crate) fn order_mates(&self, aln: &SmallVec<[Box<dyn AlignmentStream>; 2]>) -> SmallVec<[usize; 2]> {
        let len = self.records.len();
        let mut indices: SmallVec<[(u8, usize, usize, usize); 2]> = SmallVec::with_capacity(len);
        for i in 0..len {
            let r = &self.records[i];
            let start = match r.alignment_start() {
                Some(Ok(pos)) => pos.get(),
                _ => panic!("Mapped record has no alignment start"),
            };
            let tid = match r.reference_sequence_id(aln[i].header()) {
                Some(Ok(tid)) => tid,
                _ => panic!("Mapped record has no reference sequence ID"),
            };
            let flags = &self.flags[i];
            let pos = if flags.is_reverse_complemented() {
                start + r.cigar().len()
            } else {
                start
            };
           let ord = u8::from(flags.is_last_segment()) * 2 + u8::from(flags.is_secondary());
            indices.push((ord, tid, pos, i));
        }
        indices.sort();
        indices.iter().map(|t| t.3).collect()
    }
}

impl<R: Record + PartialEq> PartialOrd for FragmentState<R> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self.is_all_unmapped(), other.is_all_unmapped()) {
            (true,  true)  => Some(Ordering::Equal),
            (true,  false) => Some(Ordering::Less),
            (false, true)  => Some(Ordering::Greater),
            (false, false) => None,
        }
    }
}

#[cfg(test)]
mod tests;
