use crate::alignment::PreparedAlignmentPairIter;
use noodles::bam::record::Record;
use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;
use anyhow::Result;

#[derive(PartialEq, Debug)]
pub(crate) struct FragmentState {
    pub(crate) records: SmallVec<[Record; 2]>,
    pub(crate) species_nr: usize,
}

impl FragmentState {
    #[must_use]
    pub(crate) fn from_record(r: Record, species_nr: usize) -> Self {
        FragmentState {
            records: smallvec![r],
            species_nr,
        }
    }
    #[must_use]
    pub(crate) fn first_qname(&self) -> &[u8] {
        self.records.first().and_then(|r| r.name()).map_or(b"", |n| n.as_ref())
    }
    #[must_use]
    pub(crate) fn get_nr(&self) -> usize {
        self.species_nr
    }

    pub(crate) fn order_mates(&self) -> Result<SmallVec<[usize; 2]>> {
        let len = self.records.len();
        let mut indices: SmallVec<[(u8, usize, usize, usize); 2]> = SmallVec::with_capacity(len);
        for i in 0..len {
            let r = &self.records[i];
            let start = match r.alignment_start() {
                Some(Ok(pos)) => pos.get(),
                _ => panic!("Mapped record has no alignment start"),
            };
            let tid = match r.reference_sequence_id() {
                Some(Ok(tid)) => tid,
                _ => panic!("Mapped record has no reference sequence ID"),
            };
            let f = r.flags();
            let pos = if f.is_reverse_complemented() {
                start + r.cigar().len()
            } else {
                start
            };
           let ord = u8::from(f.is_last_segment()) * 2 + u8::from(f.is_secondary());
            indices.push((ord, tid, pos, i));
        }
        indices.sort();
        Ok(indices.iter().map(|t| t.3).collect())
    }
    #[inline]
    fn quick_unmapped_cmp(&self, a: &Record, b: &Record) -> Option<Ordering> {
        let a = a.flags();
        let b = b.flags();
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

impl PartialOrd for FragmentState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // Fast path: unmapped handling (unchanged)
        if let Some(ord) = self.quick_unmapped_cmp(&self.records[0], &other.records[0]) {
            return Some(ord);
        }

        let iter = PreparedAlignmentPairIter::new(&self.records, &other.records);

        for pair_result in iter {
            match pair_result {
                Ok(mut pair) => {
                    match pair.host_graft_are_perfect_match() {
                        (true, true) => continue, // perfect match for both => tie-break with next pair
                        (true, false) => return Some(Ordering::Greater), // host better
                        (false, true) => return Some(Ordering::Less),    // graft better
                        (false, false) => return None, // both imperfect => per-base evaluation
                    }
                },
                Err(e) => {
                    // FIXME
                    panic!("Alignment pair iteration failed: {:?}", e);
                }
            };

        }
        // All pairs matched equally perfect
        Some(Ordering::Equal)
    }
}

#[cfg(test)]
mod tests;
