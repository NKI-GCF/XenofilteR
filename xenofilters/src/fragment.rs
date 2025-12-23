use crate::alignment::PreparedAlignmentPairIter;
use rust_htslib::bam::record::{Cigar, Record};
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

    pub fn order_mates(&self) -> Vec<usize> {
        let mut indices: Vec<usize> = (0..self.records.len()).collect();
        indices.sort_by(|&a, &b| order_mates(&self.records[a]).cmp(&order_mates(&self.records[b])));
        indices
    }
}

impl PartialOrd for FragmentState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // If the record is unmapped, it is always the first record.
        if self.records[0].is_unmapped()
            && (!self.records[0].is_paired() || self.records[0].is_mate_unmapped())
        {
            if other.records[0].is_unmapped()
                && (!other.records[0].is_paired() || other.records[0].is_mate_unmapped())
            {
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
                    eprintln!(
                        "Comparing fragment states: first perfect: {}, second perfect: {}",
                        first, second
                    );
                    if first {
                        if second {
                            match balance {
                                Some(Ordering::Greater) => break,
                                Some(Ordering::Less) => break,
                                _ => Some(Ordering::Equal),
                            }
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
                        return None;
                    }
                }
                Err(_) => return Some(Ordering::Equal),
            }
        }
        balance
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

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Aux, CigarString};

    pub fn create_record(
        qname: &[u8],
        vec_cig: Vec<Cigar>,
        qual: &[u8],
        md: &str,
        is_rev: bool,
    ) -> rust_htslib::bam::Record {
        let mut record = Record::new();
        let cigar_string = CigarString(vec_cig);

        // Set basic fields
        record.set(qname, Some(&cigar_string), &vec![b'A'; qual.len()], qual);
        if is_rev {
            record.set_reverse();
        }

        // HTSlib requires the MD tag to be set manually for tests
        record.push_aux(b"MD", Aux::String(md)).unwrap();
        record
    }

    #[test]
    fn test_fragment_state_ordering() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 1);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Equal));
    }

    #[test]
    fn test_order_mates_function() {
        let qual = vec![37; 100];
        let mut rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let mut rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", true);
        rec1.set_first_in_template();
        rec2.set_last_in_template();
        let order1 = order_mates(&rec1);
        let order2 = order_mates(&rec2);
        assert!(order1 < order2); // Forward read should come before reverse read
    }

    #[test]
    fn test_fragment_state_first_qname() {
        let qual = vec![37; 100];
        let rec = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let state = FragmentState::from_record(rec, 0);
        assert_eq!(state.first_qname(), b"read1");
    }

    #[test]
    fn test_fragment_state_get_nr() {
        let qual = vec![37; 100];
        let rec = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let state = FragmentState::from_record(rec, 42);
        assert_eq!(state.get_nr(), 42);
    }

    #[test]
    fn test_fragment_state_equality() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1, state2);
    }

    #[test]
    fn test_fragment_state_inequality() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let rec2 = create_record(b"read2", vec![Cigar::Match(100)], &qual, "100", false);
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_ne!(state1, state2);
    }

    #[test]
    fn test_fragment_state_multiple_records() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let mut state = FragmentState::from_record(rec1, 0);
        state.records.push(rec2);
        assert_eq!(state.records.len(), 2);
        assert_eq!(state.first_qname(), b"read1");
    }

    #[test]
    fn test_fragment_state_order_mates_multiple_records() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", true);
        let mut state = FragmentState::from_record(rec1, 0);
        state.records.push(rec2);
        let order = state.order_mates();
        assert_eq!(order, vec![0, 1]); // Forward read should come before reverse read
    }

    #[test]
    fn test_fragment_state_partial_ord_with_unmapped() {
        let qual = vec![37; 100];
        let mut rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        rec1.set_unmapped();
        let rec2 = create_record(b"read2", vec![Cigar::Match(100)], &qual, "100", false);
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Less));
    }

    #[test]
    fn test_fragment_state_partial_ord_both_unmapped() {
        let qual = vec![37; 100];
        let mut rec1 = create_record(b"read1", vec![], &qual, "100", false);
        rec1.set_unmapped();
        let mut rec2 = create_record(b"read2", vec![], &qual, "100", false);
        rec2.set_unmapped();
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Equal));
    }

    #[test]
    fn test_fragment_state_partial_ord_no_quick_balance() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "90A10", false);
        let rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "80T20", false);
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), None); // No quick balance
    }

    #[test]
    fn test_fragment_state_partial_ord_perfect_vs_imperfect() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "90A10", false);
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Greater)); // Perfect match is better
    }

    #[test]
    fn test_fragment_state_partial_ord_imperfect_vs_perfect() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "90A10", false);
        let rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Less)); // Perfect match is better
    }

    #[test]
    fn test_fragment_state_partial_ord_equal_imperfects() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "90A10", false);
        let rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "90A10", false);
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), None); // Same imperfect matches
    }

    #[test]
    fn test_fragment_state_partial_ord_multiple_records() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "100", false);
        let rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "90A10", false);
        let mut state1 = FragmentState::from_record(rec1, 0);
        state1.records.push(create_record(
            b"read1",
            vec![Cigar::Match(100)],
            &qual,
            "100",
            false,
        ));
        let mut state2 = FragmentState::from_record(rec2, 0);
        state2.records.push(create_record(
            b"read1",
            vec![Cigar::Match(100)],
            &qual,
            "90A10",
            false,
        ));
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Greater)); // Perfect matches are better
    }

    #[test]
    fn test_fragment_state_partial_ord_multiple_records_no_quick_balance() {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "90A10", false);
        let rec2 = create_record(b"read1", vec![Cigar::Match(100)], &qual, "80T20", false);
        let mut state1 = FragmentState::from_record(rec1, 0);
        state1.records.push(create_record(
            b"read1",
            vec![Cigar::Match(100)],
            &qual,
            "85G15",
            false,
        ));
        let mut state2 = FragmentState::from_record(rec2, 0);
        state2.records.push(create_record(
            b"read1",
            vec![Cigar::Match(100)],
            &qual,
            "80T20",
            false,
        ));
        assert_eq!(state1.partial_cmp(&state2), None); // No quick balance
    }
}
