use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::Flags;
use anyhow::Result;
use crate::alignment::MdCigFlags;
use crate::alignment::SimpleRec;

#[derive(PartialEq, Debug)]
pub(crate) struct FragmentState<R> {
    flags: SmallVec<[Flags; 8]>,
    records: SmallVec<[R; 8]>,
    species_nr: usize,
}

impl<R: SimpleRec> FragmentState<R> {
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
    pub(crate) fn is_all_perfect(&self) -> Result<bool> {
        for (flags, record) in self.flags.iter().zip(self.records.iter()) {
            // Secondary alignment = split hit = penalty → not perfect
            if flags.is_secondary() {
                return Ok(false);
            }
            let mcf = MdCigFlags::try_from_record(record, flags)?;
            if !mcf.is_perfect() {
                return Ok(false);
            }
        }
        Ok(true)
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

    pub(crate) fn order_mates(&self) -> SmallVec<[usize; 2]> {
        let len = self.records.len();
        let mut indices: SmallVec<[(u8, usize, usize, usize); 2]> = SmallVec::with_capacity(len);
        for i in 0..len {
            let r = &self.records[i];
            let start = match r.alignment_start() {
                Some(Ok(pos)) => pos.get(),
                _ => panic!("Mapped record has no alignment start"),
            };
            let tid = match r.ref_seq_id() {
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

impl<R: SimpleRec> PartialOrd for FragmentState<R> {
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
