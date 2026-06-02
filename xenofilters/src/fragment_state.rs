use noodles::sam::alignment::Record;
use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;
use crate::aln_stream::AlignmentStream;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::Flags;
use anyhow::{Result, anyhow};

#[derive(PartialEq, Debug)]
struct MdCig {
    md: SmallVec<[u8; 16]>,
    cig: SmallVec<[Op; 2]>,
}

impl MdCig {
    fn is_perfect(&self) -> bool {
        self.cig.len() == 1 && !self.md.iter().all(|&b| b.is_ascii_digit())
    }
    fn new<R: Record>(r: &R) -> Result<Self> {
        let mut md = SmallVec::new();
        match r.data().get(&Tag::MISMATCHED_POSITIONS) {
            Some(Ok(Value::String(bstr))) => {
                let m: &[u8] = bstr.as_ref();
                md.extend_from_slice(m);
            },
            Some(Err(e)) => return Err(e.into()),
            Some(_) => return Err(anyhow!("Mapped record has non-string MD field")),
            None => return Err(anyhow!("Mapped record missing MD field")),
        }
        let mut cig_ops = SmallVec::new();
        for c in r.cigar().iter() {
            cig_ops.push(c.unwrap());
        }
        Ok(MdCig { md, cig: cig_ops })
    }
}

impl PartialOrd for MdCig {
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
    flags: SmallVec<[Flags; 2]>,
    ops: SmallVec<[MdCig; 2]>,
    pub(super) records: SmallVec<[R; 2]>,
    species_nr: usize,
}

impl<R: Record> FragmentState<R> {
    #[must_use]
    pub(crate) fn from_record(r: R, species_nr: usize) -> Result<Self> {
        Ok(FragmentState {
            flags: smallvec![r.flags()?],
            ops: smallvec![],
            records: smallvec![r],
            species_nr,
        })
    }
    pub(crate) fn add_record(&mut self, r: R) -> Result<()> {
        self.flags.push(r.flags()?);
        self.records.push(r);
        Ok(())
    }
    pub(crate) fn init_md_cig(&mut self) -> Result<()> {
        if self.ops.len() > 0 {
            return Ok(());
        }
        for f in 0..self.flags.len() {
            let flags = self.flags[f];
            if flags.is_secondary() || flags.is_unmapped() {
                continue;
            }
            self.ops.push(MdCig::new(&self.records[f])?);
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
            let f = self.flags[i];
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
        let a = self.flags[0];
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
        if self.ops.len() == 0 || other.ops.len() == 0 {
            // Fast path: unmapped handling
            return self.quick_unmapped_cmp(&other.flags[0]);
        }
        return self.ops[0].partial_cmp(&other.ops[0]);
    }
}

#[cfg(test)]
mod tests;
