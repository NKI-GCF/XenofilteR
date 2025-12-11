use rust_htslib::bam::record::{Record, Cigar};
use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;
use crate::alignment::{PreparedAlignmentPairIter};

pub struct FragmentState {
    pub records: SmallVec<[Record; 2]>,
}
// TODO: if there are supplementary alignments, modify cigar to skip the clipped sections.
//       may require a penalty for the placement of a secondary, elsewhere.

impl FragmentState {
    #[must_use]
    pub fn from_record(r: Record) -> Self {
        FragmentState {
            records: smallvec![r],
        }
    }
    #[must_use]
    pub fn first_qname(&self) -> &[u8] {
        self.records.first().map_or(b"", |r| r.qname())
    }

    pub fn order_mates(&self) -> Vec<usize> {
        let mut indices: Vec<usize> = (0..self.records.len()).collect();
        indices.sort_by(|&a, &b| order_mates(&self.records[a]).cmp(&order_mates(&self.records[b])));
        indices
    }
}

impl PartialEq for FragmentState {
    fn eq(&self, other: &Self) -> bool {
        self.partial_cmp(other) == Some(Ordering::Equal)
    }
}

impl PartialOrd for FragmentState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // If the record is unmapped, it is always the first record.
        if self.records[0].is_unmapped() && self.records[0].is_mate_unmapped() {
            if other.records[0].is_unmapped() && other.records[0].is_mate_unmapped() {
                return Some(Ordering::Equal);
            }
            return Some(Ordering::Less);
        }

        if other.records[0].is_unmapped() && other.records[0].is_mate_unmapped() {
            return Some(Ordering::Greater);
        }

        if self.records[0].is_unmapped() && other.records[0].is_mate_unmapped() {
            return None;
        }
        if self.records[0].is_mate_unmapped() && other.records[0].is_unmapped() {
            return None;
        }

        let iter = PreparedAlignmentPairIter::new(&self.records, &other.records);
        let mut balance = None;
        for pair_result in iter {
            balance = match pair_result {
                Ok(mut pair) => {
                    let (first, second) = pair.are_perfect_match();
                    if first {
                        if second {
                            Some(Ordering::Equal)
                        } else {
                            match balance {
                                Some(Ordering::Less) => break,
                                None => Some(Ordering::Greater),
                                _ => return Some(Ordering::Greater),
                            }
                        }
                    } else if second {
                        match balance {
                            Some(Ordering::Greater) => break,
                            None => Some(Ordering::Less),
                            _ => return Some(Ordering::Less),
                        }
                    } else {
                        // both mapped but not perfect means no quick balance.
                        break;
                    }
                }
                Err(_) => return Some(Ordering::Equal),
            }
        }
        None
    }
}

fn order_mates(r: &Record) -> (u8, u32) {
    let ord = u8::from(r.is_last_in_template()) * 2 + u8::from(r.is_secondary());
    let mut read_start = 0;
    if r.is_reverse() {
        for op in r.cigar().iter().rev() {
            match *op {
                Cigar::SoftClip(len) | Cigar::HardClip(len) => read_start += len,
                _ => break,
            }
        }
    } else {
        for op in r.cigar().iter() {
            match *op {
                Cigar::SoftClip(len) | Cigar::HardClip(len) => read_start += len,
                _ => break,
            }
        }
    }
    (ord, read_start)
}
