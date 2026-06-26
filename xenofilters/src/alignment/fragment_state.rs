use crate::alignment::MdCigFlags;
use crate::alignment::SimpleRec;
use crate::filter_algorithm::line_by_line::READ_CT;
use anyhow::Result;
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::Flags;
use smallvec::{smallvec, SmallVec};
use std::cmp::Ordering;

#[derive(PartialEq, Debug)]
pub(crate) struct FragmentState<R> {
    flags: SmallVec<[Flags; 8]>,
    records: SmallVec<[R; 8]>,
    species_nr: usize,
}

type McfPair<'f> = (SmallVec<[MdCigFlags<'f>; 8]>, SmallVec<[MdCigFlags<'f>; 8]>);

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
    pub(crate) fn get_records(&self) -> &[R] {
        &self.records
    }
    pub(crate) fn drain_records(&mut self) -> impl Iterator<Item = R> + '_ {
        self.records.drain(..)
    }
    #[must_use]
    pub(crate) fn first_qname(&self) -> &[u8] {
        self.records
            .first()
            .and_then(|r| r.name())
            .map_or(b"", |n| n.as_ref())
    }
    #[must_use]
    pub(crate) fn get_nr(&self) -> usize {
        self.species_nr
    }
    pub(crate) fn flags(&self, i: usize) -> Option<&Flags> {
        self.flags.get(i)
    }
    // Sufficient check: for properly paired BAM the mate-unmapped flag mirrors
    // the mate's state, so reading only the first record's flags is valid.
    // Supplementary records cannot be unmapped by definition; they are never
    // present when all primaries are unmapped (SA:Z is absent on unmapped primaries).
    fn is_all_unmapped(&self) -> bool {
        let f = &self.flags[0];
        f.is_unmapped() && (!f.is_segmented() || f.is_mate_unmapped())
    }

    pub(crate) fn order_mates(&self) -> SmallVec<[usize; READ_CT]> {
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
    pub(crate) fn cmp_perfect<'f>(
        &'f self,
        other: &'f FragmentState<R>,
        ord: &mut Option<Ordering>,
    ) -> Result<McfPair<'f>> {
        let mut mcfs1: SmallVec<[MdCigFlags<'f>; 8]> =
            SmallVec::with_capacity(self.get_records().len());
        let mut perfect_self = true;
        for (flags, record) in self.flags.iter().zip(self.records.iter()) {
            let mcf = MdCigFlags::try_from_record(record, flags)?;
            if flags.is_supplementary() || !mcf.is_perfect() {
                perfect_self = false;
            }
            mcfs1.push(mcf);
        }

        let mut mcfs2: SmallVec<[MdCigFlags<'f>; 8]> =
            SmallVec::with_capacity(other.get_records().len());
        let mut perfect_other = true;
        for (flags, record) in other.flags.iter().zip(other.records.iter()) {
            let mcf = MdCigFlags::try_from_record(record, flags)?;
            if flags.is_supplementary() || !mcf.is_perfect() {
                perfect_other = false;
            }
            mcfs2.push(mcf);
        }

        *ord = match (perfect_self, perfect_other) {
            (true, true) => Some(Ordering::Equal),
            (false, true) => Some(Ordering::Less),
            (true, false) => Some(Ordering::Greater),
            (false, false) => None,
        };
        Ok((mcfs1, mcfs2))
    }
}

impl<R: SimpleRec> PartialOrd for FragmentState<R> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self.is_all_unmapped(), other.is_all_unmapped()) {
            (true, true) => Some(Ordering::Equal),
            (true, false) => Some(Ordering::Less),
            (false, true) => Some(Ordering::Greater),
            (false, false) => None,
        }
    }
}

#[cfg(test)]
mod tests;
