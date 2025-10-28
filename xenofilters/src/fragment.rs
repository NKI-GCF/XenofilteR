use anyhow::Result;
use rust_htslib::bam::record::Record;
use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;

use crate::MAX_Q;
use crate::alignment::{PrepareError, PreparedAlignmentPair, PreparedAlignmentPairIter};

type PreparedAlignmentPairIterBox<'a> =
    Box<dyn Iterator<Item = Result<PreparedAlignmentPair<'a>, PrepareError>> + 'a>;

pub enum Evaluation<'a> {
    Greater,
    Less,
    Equal,
    NeedsScoring(PreparedAlignmentPairIterBox<'a>),
}

pub struct FragmentState {
    records: SmallVec<[Record; 2]>,
    log_likelihood_mismatch: [f64; MAX_Q + 2],
}
// TODO: if there are supplementary alignments, modify cigar to skip the clipped sections.
//       may require a penalty for the placement of a secondary, elsewhere.

impl FragmentState {
    #[must_use]
    pub fn from_record(r: Record, log_likelihood_mismatch: [f64; MAX_Q + 2]) -> Self {
        FragmentState {
            records: smallvec![r],
            log_likelihood_mismatch,
        }
    }
    #[must_use]
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
    fn try_needs_score(&self, iter: PreparedAlignmentPairIterBox<'_>) -> Option<f64> {
        let mut total_score_diff = Some(0.0);

        // We filter_map to skip pairs that had a PrepareError,
        // as they won't contribute to the score.
        for pair_result in iter.filter_map(std::result::Result::ok) {
            match pair_result.score(&self.log_likelihood_mismatch) {
                Ok(score_diff) => total_score_diff = total_score_diff.map(|n| n + score_diff),
                Err(_) => return None,
            }
        }
        total_score_diff
    }
    fn try_quick_cmp<'a>(&'a self, other: &'a Self) -> Evaluation<'a> {
        if self.records[0].is_unmapped() && self.records[0].is_mate_unmapped() {
            if other.records[0].is_unmapped() && other.records[0].is_mate_unmapped() {
                return Evaluation::Equal;
            }
            return Evaluation::Less;
        }
        if other.records[0].is_unmapped() && other.records[0].is_mate_unmapped() {
            return Evaluation::Greater;
        }

        let iter = PreparedAlignmentPairIter::new(&self.records, &other.records);

        if self.records[0].is_unmapped() && other.records[0].is_mate_unmapped() {
            return Evaluation::NeedsScoring(Box::new(iter));
        }
        if other.records[0].is_unmapped() && self.records[0].is_mate_unmapped() {
            return Evaluation::NeedsScoring(Box::new(iter));
        }

        let mut balance = None;
        for pair_result in iter {
            balance = match pair_result {
                Ok(pair) => {
                    let s_perfect = pair.is_perfect_match(true);
                    let o_perfect = pair.is_perfect_match(false);

                    if s_perfect {
                        if o_perfect {
                            Some(Evaluation::Equal)
                        } else {
                            match balance {
                                Some(Evaluation::Less) => break,
                                None => Some(Evaluation::Greater),
                                _ => return Evaluation::Greater,
                            }
                        }
                    } else if o_perfect {
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
                Err(_) => {
                    // A PrepareError occurred (e.g., No MD tag).
                    // We can't do a quick check, so break and go to full scoring.
                    break;
                }
            }
        }
        Evaluation::NeedsScoring(Box::new(PreparedAlignmentPairIter::new(
            &self.records,
            &other.records,
        )))
    }
}

impl Ord for FragmentState {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.try_quick_cmp(other) {
            Evaluation::Greater => Ordering::Greater,
            Evaluation::Less => Ordering::Less,
            Evaluation::Equal => Ordering::Equal,
            Evaluation::NeedsScoring(pair_iter) => {
                match self.try_needs_score(pair_iter) {
                    Some(score_diff) => {
                        // Compare the final score difference to 0.0
                        score_diff.partial_cmp(&0.0).unwrap_or(Ordering::Equal)
                    }
                    None => Ordering::Equal, // Error case
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
