use std::cmp::Ordering;
use rust_htslib::bam::record::Record;
use smallvec::{SmallVec, smallvec};
use anyhow::Result;
use crate::MAX_Q;
use crate::alignment::{PreparedAlignment, PreparedAlignmentIter, PrepareError};

pub enum Evaluation<'a> {
    PerfectMatch,
    BothUnmapped,
    NeedsScoring(Box<dyn Iterator<Item = PreparedAlignment<'a>> + 'a>),
}

pub struct FragmentState {
    alns: SmallVec<[Record; 2]>,
    log_likelihood_mismatch: [f64; MAX_Q + 2],
}

impl FragmentState {
    pub fn from_record(r: Record, log_likelihood_mismatch: [f64; MAX_Q + 2]) -> Self {
        FragmentState {
            alns: smallvec![r],
            log_likelihood_mismatch,
        }
    }
    pub fn first_qname(&self) -> &[u8] {
        self.alns.first().map_or(b"", |r| r.qname())
    }

    pub fn add_record(&mut self, r: Record) {
        self.alns.push(r);
    }
    pub fn drain(&mut self) -> SmallVec<[Record; 2]> {
        self.alns.drain(..).collect()
    }
    fn pre_eval_fragment(&self) -> Result<Evaluation<'_>, PrepareError> {
        if self.alns[0].is_unmapped() && self.alns[0].is_mate_unmapped() {
            return Ok(Evaluation::BothUnmapped);
        }

        // First pass: check for perfect matches and collect errors
        for prepared_result in PreparedAlignmentIter::new(&self.alns) {
            if prepared_result?.is_perfect_match() {
                return Ok(Evaluation::PerfectMatch);
            }
        }

        // Second pass: create iterator for scoring
        let iter = PreparedAlignmentIter::new(&self.alns)
            .filter_map(|r| r.ok());

        Ok(Evaluation::NeedsScoring(Box::new(iter)))
    }
    /// Returns the score or None if there was an error (logged to stderr)
    fn try_score(&self) -> Option<f64> {
        match self.pre_eval_fragment() {
            Err(e) => {
                eprintln!("Error preparing fragment: {e}");
                None
            }
            Ok(Evaluation::PerfectMatch) => Some(f64::INFINITY),
            Ok(Evaluation::BothUnmapped) => Some(f64::NEG_INFINITY),
            Ok(Evaluation::NeedsScoring(iter)) => {
                let mut score = 0.0;
                for p in iter {
                    match p.score(&self.log_likelihood_mismatch) {
                        Ok(sc) => score += sc,
                        Err(e) => {
                            eprintln!("Error scoring: {}", e);
                            return None;
                        }
                    }
                }
                Some(score)
            }
        }
    }
}

impl Ord for FragmentState {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self.try_score(), other.try_score()) {
            (Some(s), Some(o)) => s.partial_cmp(&o).unwrap_or(Ordering::Equal),
            (None, None) => Ordering::Equal,
            (None, Some(_)) => Ordering::Less,
            (Some(_), None) => Ordering::Greater,
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
