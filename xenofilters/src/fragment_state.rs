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
        /*#[cfg(test)]
        eprintln!(
            "Comparing FragmentStates:\nself unmapped: {}, other unmapped: {}, self single/mate_unmapped: {}, other single/mate_unmapped: {}",
            self.records[0].is_unmapped(),
            other.records[0].is_unmapped(),
            self.records[0].is_paired() == false || self.records[0].is_mate_unmapped(),
            other.records[0].is_paired() == false || other.records[0].is_mate_unmapped(),
        );*/
        // If the record is unmapped, it is always the first record.
        if self.records[0].is_unmapped()
            && (!self.records[0].is_paired() || self.records[0].is_mate_unmapped())
        {
            if other.records[0].is_unmapped()
                && (!other.records[0].is_paired() || other.records[0].is_mate_unmapped())
            {
                //#[cfg(test)]
                //eprintln!("both all reads unmapped");
                return Some(Ordering::Equal);
            }
            //#[cfg(test)]
            //eprintln!("self totally unmapped, other not");
            return Some(Ordering::Less);
        }

        if other.records[0].is_unmapped()
            && (!other.records[0].is_paired() || other.records[0].is_mate_unmapped())
        {
            //#[cfg(test)]
            //eprintln!("other totally unmapped, self not");
            return Some(Ordering::Greater);
        }

        if self.records[0].is_unmapped() && other.records[0].is_mate_unmapped() {
            //#[cfg(test)]
            //eprintln!("distinct records are mapped");
            return None;
        }
        if self.records[0].is_mate_unmapped() && other.records[0].is_unmapped() {
            //#[cfg(test)]
            //eprintln!("distinct records are mapped");
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
                    /*#[cfg(test)]
                    eprintln!(
                        "Comparing fragment states:first perfect: {first}, second perfect: {second}",
                    );*/
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
                        // both mapped but not perfect means no quick winner.
                        return None;
                    }
                }
                Err(e) => panic!("Error during prepared alignment pair iteration: {:?}", e),
            }
        }
        balance
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
