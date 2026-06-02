use noodles::sam::alignment::Record;
use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;
use crate::aln_stream::AlignmentStream;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::Flags;
use anyhow::{Result, anyhow};

#[derive(PartialEq, Debug, Default)]
pub(super) struct MdCigFlags {
    pub(super) flags: Flags,
    pub(super) md: SmallVec<[u8; 16]>,
    pub(super) cig: SmallVec<[Op; 4]>,
}

impl MdCigFlags {
    fn is_perfect(&self) -> bool {
        self.cig.len() == 1 && !self.md.iter().all(|&b| b.is_ascii_digit())
    }
    fn new(flags: Flags) -> Self {
        Self { flags, ..Default::default() }
    }
    fn complete<R: Record>(&mut self, r: &R) -> Result<()> {
        match r.data().get(&Tag::MISMATCHED_POSITIONS) {
            Some(Ok(Value::String(bstr))) => {
                let m: &[u8] = bstr.as_ref();
                self.md.extend_from_slice(m);
            },
            Some(Err(e)) => return Err(e.into()),
            Some(_) => return Err(anyhow!("Mapped record has non-string MD field")),
            None => return Err(anyhow!("Mapped record missing MD field")),
        }
        for c in r.cigar().iter() {
            self.cig.push(c?);
        }
        Ok(())
    }
}

impl PartialOrd for MdCigFlags {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self.is_perfect(), other.is_perfect()) {
            (true, true) => Some(Ordering::Equal), // perfect match for both => tie-break with next pair
            (true, false) => Some(Ordering::Greater), // self better
            (false, true) => Some(Ordering::Less),    // other better
            (false, false) => None, // Slow path: both imperfect => per-base evaluation
        }
    }
}

#[derive(PartialEq, Debug)]
pub(crate) struct FragmentState<R> {
    pub(super) ops: SmallVec<[MdCigFlags; 2]>,
    pub(super) records: SmallVec<[R; 2]>,
    species_nr: usize,
}

impl<R: Record> FragmentState<R> {
    pub(crate) fn from_record(r: R, species_nr: usize) -> Result<Self> {
        Ok(FragmentState {
            ops: smallvec![MdCigFlags::new(r.flags()?)],
            records: smallvec![r],
            species_nr,
        })
    }
    pub(crate) fn add_record(&mut self, r: R) -> Result<()> {
        self.ops.push(MdCigFlags::new(r.flags()?));
        self.records.push(r);
        Ok(())
    }
    pub(crate) fn init_md_cig(&mut self) -> Result<()> {
        if !self.ops.is_empty() {
            return Ok(());
        }
        for f in 0..self.ops.len() {
            let record = &self.records[f];
            let flags = record.flags()?;
            if !flags.is_secondary() && !flags.is_unmapped() {
                self.ops[f].complete(record)?;
            }
        }
        Ok(())
    }
    #[must_use]
    pub(crate) fn first_qname(&self) -> &[u8] {
        self.records.first().and_then(|r| r.name()).map_or(b"", |n| n.as_ref())
    }
    #[must_use]
    pub(crate) fn get_nr(&self) -> usize {
        self.species_nr
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
            let f = self.ops[i].flags;
            let pos = if f.is_reverse_complemented() {
                start + r.cigar().len()
            } else {
                start
            };
           let ord = u8::from(f.is_last_segment()) * 2 + u8::from(f.is_secondary());
            indices.push((ord, tid, pos, i));
        }
        indices.sort();
        indices.iter().map(|t| t.3).collect()
    }
    #[inline]
    fn quick_unmapped_cmp(&self, b: &Flags) -> Option<Ordering> {
        let a = self.ops[0].flags;
        if a.is_unmapped() && (!a.is_segmented() || a.is_mate_unmapped()) {
            if b.is_unmapped() && (!b.is_segmented() || b.is_mate_unmapped()) {
                Some(Ordering::Equal)
            } else {
                Some(Ordering::Less)
            }
        } else if b.is_unmapped() && (!b.is_segmented() || b.is_mate_unmapped()) {
            Some(Ordering::Greater)
        } else {
            None
        }
    }
}

impl<R: Record + PartialEq> PartialOrd for FragmentState<R> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.ops.is_empty() || other.ops.is_empty() {
            // Fast path: unmapped handling
            return self.quick_unmapped_cmp(&other.ops[0].flags);
        }
        self.ops[0].partial_cmp(&other.ops[0])
    }
}

#[cfg(test)]
mod tests;
