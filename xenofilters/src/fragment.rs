use std::cmp::Ordering;
use rust_htslib::bam::record::Record;
use smallvec::{SmallVec, smallvec};
use anyhow::Result;
use crate::MAX_Q;
use crate::alignment::{PreparedAlignment, PreparedAlignmentIter, PrepareError};

type PreparedAlignmentIterBox<'a> = Box<dyn Iterator<Item = Result<PreparedAlignment<'a>, PrepareError>> + 'a>;

pub enum Evaluation<'a> {
    Greater,
    Less,
    Equal,
    NeedsScoring(PreparedAlignmentIterBox<'a>, PreparedAlignmentIterBox<'a>),
}

pub struct FragmentState {
    records: SmallVec<[Record; 2]>,
    log_likelihood_mismatch: [f64; MAX_Q + 2],
}
// TODO2: if there are supplementary alignments, modify sigar to skip the clipped sections.
//        may require a penalty for the placement of a secondary, elsewhere.


impl FragmentState {
    pub fn from_record(r: Record, log_likelihood_mismatch: [f64; MAX_Q + 2]) -> Self {
        FragmentState {
            records: smallvec![r],
            log_likelihood_mismatch,
        }
    }
    pub fn first_qname(&self) -> &[u8] {
        self.records.first().map_or(b"", |r| r.qname())
    }

    pub fn add_record(&mut self, r: Record) {
        if r.is_secondary() || r.is_supplementary() {
            self.records.push(r);
        } else if r.is_first_in_template() {
            self.records.insert(0, r);
        } else if r.is_last_in_template() && !self.records.is_empty() {
            self.records.insert(1, r);
        } else {
            self.records.push(r);
        }
    }
    pub fn drain(&mut self) -> SmallVec<[Record; 2]> {
        self.records.drain(..).collect()
    }

    fn try_needs_score<'a>(&self, iter: PreparedAlignmentIterBox<'a>, other_iter: PreparedAlignmentIterBox<'a>) -> (Option<f64>, Option<f64>) {
        let mut ret1 = Some(0.0);
        for p in iter.filter_map(|x| x.ok()) {
            match p.score(&self.log_likelihood_mismatch) {
                Ok(sc) => ret1 = ret1.map(|n| n + sc),
                Err(_) => {
                    ret1 = None;
                    break;
                }
            }
        }
        let mut ret2 = Some(0.0);
        for o in other_iter.filter_map(|x| x.ok()) {
            match o.score(&self.log_likelihood_mismatch) {
                Ok(sc) => ret2 = ret2.map(|n| n + sc),
                Err(_) => {
                    ret2 = None;
                    break;
                }
            }
        }
        (ret1, ret2)
    }
    fn try_quick_cmp<'a>(&'a self, other: &'a Self) -> Evaluation<'a> {
        if self.records[0].is_unmapped() && self.records[0].is_mate_unmapped() {
            if other.records[0].is_unmapped() && other.records[0].is_mate_unmapped() {
                return Evaluation::Equal;
            } else {
                return Evaluation::Less;
            }
        }
        if other.records[0].is_unmapped() && other.records[0].is_mate_unmapped() {
            return Evaluation::Greater;
        }

        let self_iter = PreparedAlignmentIter::new(&self.records);
        let other_iter = PreparedAlignmentIter::new(&other.records);

        if self.records[0].is_unmapped() && other.records[0].is_mate_unmapped() {
            return Evaluation::NeedsScoring(Box::new(self_iter), Box::new(other_iter));
        }
        if other.records[0].is_unmapped() && self.records[0].is_mate_unmapped() {
            return Evaluation::NeedsScoring(Box::new(self_iter), Box::new(other_iter));
        }

        let mut balance = None;
        for (s, o) in self_iter.zip(other_iter) {
            balance = match (s, o) {
                (Ok(s), Ok(o)) => {
                    if s.is_perfect_match() {
                        if !o.is_perfect_match() {
                            match balance {
                                Some(Evaluation::Less) => break,
                                None => Some(Evaluation::Greater),
                                _ => return Evaluation::Greater,
                            }
                        } else {
                            Some(Evaluation::Equal)
                        }
                    } else if o.is_perfect_match() {
                        match balance {
                            Some(Evaluation::Greater) => break,
                            None => Some(Evaluation::Less),
                            _ => return Evaluation::Less,
                        }
                    } else {
                        // both mapped but not perfect means no quick balance.
                        break;
                    }
                }
                (Ok(_s), Err(_)) => match balance {
                    Some(Evaluation::Less) => break,
                    None => Some(Evaluation::Greater),
                    _ => return Evaluation::Greater,
                }
                (Err(_), Ok(_o)) => match balance {
                    Some(Evaluation::Greater) => break,
                    None => Some(Evaluation::Less),
                    _ => return Evaluation::Less,
                }
                (Err(_), Err(_)) => { // both unmapped
                    if let Some(balance) = balance {
                        return balance;
                    }
                    Some(Evaluation::Equal)
                }
            }
        }
        Evaluation::NeedsScoring(
            Box::new(PreparedAlignmentIter::new(&self.records)),
            Box::new(PreparedAlignmentIter::new(&other.records)),
        )
    }
}

impl Ord for FragmentState {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.try_quick_cmp(other) {
            Evaluation::Greater => Ordering::Greater,
            Evaluation::Less => Ordering::Less,
            Evaluation::Equal => Ordering::Equal,
            Evaluation::NeedsScoring(s, o) => {
                match self.try_needs_score(s, o) {
                    (Some(s), Some(o)) => s.partial_cmp(&o).unwrap_or(Ordering::Equal),
                    (Some(_), None) => Ordering::Greater,
                    (None, Some(_)) => Ordering::Less,
                    (None, None) => Ordering::Equal,
                }
            }
        }
    }
}

impl PartialEq for FragmentState {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl PartialOrd for FragmentState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Eq for FragmentState {}
