use crate::alignment::PreparedAlignmentPairIter;
use rust_htslib::bam::record::Record;
use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;

#[derive(PartialEq, Debug)]
pub struct FragmentState {
    pub records: SmallVec<[Record; 2]>,
    pub species_nr: usize,
}

impl FragmentState {
    #[must_use]
    pub fn from_record(r: Record, species_nr: usize) -> Self {
        FragmentState {
            records: smallvec![r],
            species_nr,
        }
    }
    #[must_use]
    pub fn first_qname(&self) -> &[u8] {
        self.records.first().map_or(b"", |r| r.qname())
    }
    #[must_use]
    pub fn get_nr(&self) -> usize {
        self.species_nr
    }

    pub fn order_mates(&self) -> SmallVec<[usize; 2]> {
        let mut indices: SmallVec<[usize; 2]> = (0..self.records.len()).collect();
        indices.sort_by(|&a, &b| order_mates(&self.records[a]).cmp(&order_mates(&self.records[b])));
        indices
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
            let mut pair = match pair_result {
                Ok(p) => p,
                Err(_) => return None, // per-base evaluate
            };

            match pair.host_graft_are_perfect_match() {
                (true, true) => continue, // perfect match for both => tie-break with next pair
                (true, false) => return Some(Ordering::Greater), // host better
                (false, true) => return Some(Ordering::Less),    // graft better
                (false, false) => return None, // both imperfect => per-base evaluation
            }
        }
        // All pairs matched equally perfect
        Some(Ordering::Equal)
    }
}

impl FragmentState {
    #[inline]
    fn quick_unmapped_cmp(&self, a: &Record, b: &Record) -> Option<Ordering> {
        if a.is_unmapped() && (!a.is_paired() || a.is_mate_unmapped()) {
            if b.is_unmapped() && (!b.is_paired() || b.is_mate_unmapped()) {
                Some(Ordering::Equal)
            } else {
                Some(Ordering::Less)
            }
        } else if b.is_unmapped() && (!b.is_paired() || b.is_mate_unmapped()) {
            Some(Ordering::Greater)
        } else {
            None
        }
    }
}

fn order_mates(r: &Record) -> (u8, i32, i64) {
    let ord = u8::from(r.is_last_in_template()) * 2 + u8::from(r.is_secondary());
    let pos = if r.is_reverse() {
        r.cigar().end_pos()
    } else {
        r.pos()
    };
    (ord, r.tid(), pos)
}

#[cfg(test)]
mod tests;
