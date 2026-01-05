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
            && (self.records[0].is_mate_unmapped() || !self.records[0].is_paired())
        {
            if other.records[0].is_unmapped()
                && (other.records[0].is_mate_unmapped() || !other.records[0].is_paired())
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
                    #[cfg(test)]
                    let pair_str = format!("{:?}", pair);
                    let (first, second) = pair.are_perfect_match();
                    #[cfg(test)]
                    eprintln!(
                        "Comparing fragment states:first perfect: {first}, second perfect: {second}",
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
    use crate::tests::create_record;
    use anyhow::Result;

    #[test]
    fn test_fragment_state_ordering() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 1);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Equal));
        Ok(())
    }

    #[test]
    fn test_order_mates_function() -> Result<()> {
        let qual = vec![37; 100];
        let mut rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let mut rec2 = create_record(b"read1", "100M", &[], &qual, "100", true)?;
        rec1.set_first_in_template();
        rec2.set_last_in_template();
        let order1 = order_mates(&rec1);
        let order2 = order_mates(&rec2);
        assert!(order1 < order2); // Forward read should come before reverse read
        Ok(())
    }

    #[test]
    fn test_fragment_state_first_qname() -> Result<()> {
        let qual = vec![37; 100];
        let rec = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let state = FragmentState::from_record(rec, 0);
        assert_eq!(state.first_qname(), b"read1");
        Ok(())
    }

    #[test]
    fn test_fragment_state_get_nr() -> Result<()> {
        let qual = vec![37; 100];
        let rec = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let state = FragmentState::from_record(rec, 42);
        assert_eq!(state.get_nr(), 42);
        Ok(())
    }

    #[test]
    fn test_fragment_state_equality() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1, state2);
        Ok(())
    }

    #[test]
    fn test_fragment_state_inequality() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let rec2 = create_record(b"read2", "100M", &[], &qual, "100", false)?;
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_ne!(state1, state2);
        Ok(())
    }

    #[test]
    fn test_fragment_state_multiple_records() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let mut state = FragmentState::from_record(rec1, 0);
        state.records.push(rec2);
        assert_eq!(state.records.len(), 2);
        assert_eq!(state.first_qname(), b"read1");
        Ok(())
    }

    #[test]
    fn test_fragment_state_order_mates_multiple_records() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let rec2 = create_record(b"read1", "100M", &[], &qual, "100", true)?;
        let mut state = FragmentState::from_record(rec1, 0);
        state.records.push(rec2);
        let order = state.order_mates();
        assert_eq!(order, vec![0, 1]); // Forward read should come before reverse read
        Ok(())
    }

    #[test]
    fn test_fragment_state_partial_ord_with_unmapped() -> Result<()> {
        let qual = vec![37; 100];
        let seq = vec![b'A'; 100];
        let mut rec1 = create_record(b"read1", "", &seq, &qual, "", false)?;
        rec1.set_unmapped();
        let rec2 = create_record(b"read2", "100M", &seq, &qual, "", false)?;
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Less));
        Ok(())
    }

    #[test]
    fn test_fragment_state_partial_ord_both_unmapped() -> Result<()> {
        let qual = vec![37; 100];
        let seq = vec![b'A'; 100];
        let mut rec1 = create_record(b"read1", "", &seq, &qual, "", false)?;
        rec1.set_unmapped();
        let mut rec2 = create_record(b"read2", "", &seq, &qual, "", false)?;
        rec2.set_unmapped();
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Equal));
        Ok(())
    }

    #[test]
    fn test_fragment_state_partial_ord_no_quick_balance() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
        let rec2 = create_record(b"read1", "100M", &[], &qual, "80T20", false)?;
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), None); // No quick balance
        Ok(())
    }

    #[test]
    fn test_fragment_state_partial_ord_perfect_vs_imperfect() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let rec2 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Greater)); // Perfect match is better
        Ok(())
    }

    #[test]
    fn test_fragment_state_partial_ord_imperfect_vs_perfect() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
        let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Less)); // Perfect match is better
        Ok(())
    }

    #[test]
    fn test_fragment_state_partial_ord_equal_imperfects() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
        let rec2 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
        let state1 = FragmentState::from_record(rec1, 0);
        let state2 = FragmentState::from_record(rec2, 0);
        assert_eq!(state1.partial_cmp(&state2), None); // Same imperfect matches
        Ok(())
    }

    #[test]
    fn test_fragment_state_partial_ord_multiple_records() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
        let rec2 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
        let mut state1 = FragmentState::from_record(rec1, 0);
        state1
            .records
            .push(create_record(b"read1", "100M", &[], &qual, "100", false)?);
        let mut state2 = FragmentState::from_record(rec2, 0);
        state2
            .records
            .push(create_record(b"read1", "100M", &[], &qual, "90A10", false)?);
        assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Greater)); // Perfect matches are better
        Ok(())
    }

    #[test]
    fn test_fragment_state_partial_ord_multiple_records_no_quick_balance() -> Result<()> {
        let qual = vec![37; 100];
        let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
        let rec2 = create_record(b"read1", "100M", &[], &qual, "80T20", false)?;
        let mut state1 = FragmentState::from_record(rec1, 0);
        state1
            .records
            .push(create_record(b"read1", "100M", &[], &qual, "85G15", false)?);
        let mut state2 = FragmentState::from_record(rec2, 0);
        state2
            .records
            .push(create_record(b"read1", "100M", &[], &qual, "80T20", false)?);
        assert_eq!(state1.partial_cmp(&state2), None); // No quick balance
        Ok(())
    }
}
