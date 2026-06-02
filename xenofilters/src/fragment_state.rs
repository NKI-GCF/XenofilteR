use noodles::sam::alignment::Record;
use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;
use crate::aln_stream::AlignmentStream;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::record::cigar::Op;

#[derive(PartialEq, Debug)]
pub(crate) struct FragmentState<R> {
    pub(crate) records: SmallVec<[R; 2]>,
    pub(crate) species_nr: usize,
}

impl<R: Record> FragmentState<R> {
    #[must_use]
    pub(crate) fn from_record(r: R, species_nr: usize) -> Self {
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
            let f = r.flags().expect("Bug: unchecked flags in order_mates");
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
    fn quick_unmapped_cmp(&self, a: &R, b: &R) -> Option<Ordering> {
        let a = a.flags().expect("Bug: unchecked host flags");
        let b = b.flags().expect("Bug: unchecked graft flags");
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
        // Fast path: unmapped handling
        if let Some(ord) = self.quick_unmapped_cmp(&self.records[0], &other.records[0]) {
            return Some(ord);
        }
        let mut cig_op: [SmallVec<[Op; 2]>; 2] = [SmallVec::new(), SmallVec::new()];
        let mut md: [SmallVec<[SmallVec<[u8; 4]>; 2]>; 2] = [SmallVec::new(), SmallVec::new()];
        for r in &self.records {
            let flags = r.flags().expect("Bug: unchecked Host flags");
            if flags.is_secondary() || flags.is_unmapped() {
                continue;
            }
            match r.data().get(&Tag::MISMATCHED_POSITIONS) {
                Some(Ok(Value::String(bstr))) => {
                    let m: &[u8] = bstr.as_ref();
                    md[0].push(m.into())
                },
                _ => panic!("Host record missing MD field"),
            }
            for c in r.cigar().iter() {
                cig_op[0].push(c.unwrap());
            }
        }
        for r in &other.records {
            let flags = r.flags().expect("Bug: unchecked Graft flags");
            if flags.is_secondary() || flags.is_unmapped() {
                continue;
            }
            match r.data().get(&Tag::MISMATCHED_POSITIONS) {
                Some(Ok(Value::String(bstr))) => {
                    let m: &[u8] = bstr.as_ref();
                    md[1].push(m.into())
                },
                _ => panic!("Graft record missing MD field"),
            }
            for c in r.cigar().iter() {
                cig_op[1].push(c.unwrap());
            }
        }
        let host_perfect = cig_op.len() == 1 && !md[0].iter().all(|m| m.iter().all(|&b| b.is_ascii_digit()));
        let graft_perfect = cig_op.len() == 1 && !md[1].iter().all(|m| m.iter().all(|&b| b.is_ascii_digit()));
        match (host_perfect, graft_perfect) {
            (true, true) => Some(Ordering::Equal), // perfect match for both => tie-break with next pair
            (true, false) => Some(Ordering::Greater), // host better
            (false, true) => Some(Ordering::Less),    // graft better
            (false, false) => None, // Slow path: both imperfect => per-base evaluation
        }
    }
}

#[cfg(test)]
mod tests;
