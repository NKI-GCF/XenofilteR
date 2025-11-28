use anyhow::Result;
use rust_htslib::bam::record::{Aux, Record};

use super::{
    AlignmentCompareIterator, AlignmentError, AlnCmpOp, MdOpIterator,
    PrepareError,
};
use crate::{LOG_LIKELIHOOD_MISMATCH, CONFIG, LOG_LIKELIHOOD_MATCH};
use crate::alignment::{UnifiedOpIterator, UnifiedOp};

#[allow(dead_code)]
pub fn print_req(i: usize, rec: &Record) {
    let qname = String::from_utf8_lossy(rec.qname());
    let cigar = rec.cigar().to_string();
    eprint!("{i}:{qname}\t{cigar}");
    if let Ok(Aux::String(md)) = rec.aux(b"MD") {
        eprint!("\tMD:Z:{md}");
    }
    eprintln!();
}

pub struct PreparedAlignmentPair<'a> {
    iter1: UnifiedOpIterator<'a>, // host
    iter2: UnifiedOpIterator<'a>, // graft
    qual: &'a [u8],
}

impl<'a> PreparedAlignmentPair<'a> {
    #[must_use]
    pub fn is_perfect_match(&mut self, first: bool) -> bool {
        if first {
            match self.iter1.peek() {
                Some(UnifiedOp::Match(len)) => *len == self.qual.len(),
                _ => false,
            }
        } else {
            match self.iter2.peek() {
                Some(UnifiedOp::Match(len)) => *len == self.qual.len(),
                _ => false,
            }
        }
    }
    fn score_op(&self, op: &UnifiedOp, read_i: &mut usize) -> Result<f64, PrepareError> {
        let log_likelihood_mismatch = LOG_LIKELIHOOD_MISMATCH.get().unwrap();
        let config = CONFIG.get().unwrap();
        let mut op_score = 0.0;

        match op {
            UnifiedOp::Match(len) => {
                for _ in 0..*len {
                    let q = *self
                        .qual
                        .get(*read_i)
                        .ok_or(AlignmentError::QualIndexOutOfBounds)?;
                    op_score += LOG_LIKELIHOOD_MATCH[q as usize];
                    *read_i += 1;
                }
            }
            UnifiedOp::Mis(len) => {
                for _ in 0..*len {
                    let q = *self
                        .qual
                        .get(*read_i)
                        .ok_or(AlignmentError::QualIndexOutOfBounds)?;
                    op_score += log_likelihood_mismatch[q as usize];
                    *read_i += 1;
                }
            },
            UnifiedOp::Ins(len) => {
                op_score += config.gap_open + config.gap_extend * (*len as f64);
                *read_i += *len as usize;
            }
            UnifiedOp::Del(len) => op_score += config.gap_open + config.gap_extend * (*len as f64),
            UnifiedOp::Deletion(seq) => op_score += config.gap_open + config.gap_extend * (seq.len() as f64),
            UnifiedOp::RefSkip(_) => {},
            UnifiedOp::Trans(score) | UnifiedOp::Translocate { score, ..} => op_score += score,
            UnifiedOp::Mate | UnifiedOp::MateSwitch {..} => {},
            UnifiedOp::Mismatch(seq) => {
                for _ in 0..seq.len() {
                    let q = *self
                        .qual
                        .get(*read_i)
                        .ok_or(AlignmentError::QualIndexOutOfBounds)?;
                    op_score += log_likelihood_mismatch[q as usize];
                    *read_i += 1;
                }
            }
            UnifiedOp::Insertion(seq) => {
                let len = seq.len();
                op_score += config.gap_open + config.gap_extend * (len as f64);
                *read_i += len;
            }
        }
        Ok(op_score)
    }

    fn score_aln_op_against_mm(
        &self,
        iter: UnifiedOpIterator<'a>,
    ) -> Result<f64, PrepareError> {
        let mut read_i = 0;
        let mut tot_score = 0.0;
        let log_likelihood_mismatch = LOG_LIKELIHOOD_MISMATCH.get().unwrap();
        let config = CONFIG.get().unwrap();

        for op in iter {
            let mut len = read_i;
            tot_score += self.score_op(&op, &mut read_i)?;
            len = read_i - len;
            for i in 0..len {
                let q = *self
                    .qual
                    .get(read_i - len + i)
                    .ok_or(AlignmentError::QualIndexOutOfBounds)?;
                tot_score -= log_likelihood_mismatch[q as usize];
            }

            match op {
                UnifiedOp::Match(len) => {
                    for _ in 0..*len {
                        let q = *self
                            .qual
                            .get(read_i)
                            .ok_or(AlignmentError::QualIndexOutOfBounds)?;
                        tot_score += LOG_LIKELIHOOD_MATCH[q as usize] - log_likelihood_mismatch[q as usize];
                        read_i += 1;
                    }
                }
                UnifiedOp::Mis(len) => read_i += *len as usize,
                UnifiedOp::Ins(len) => {
                    tot_score += config.gap_open + config.gap_extend * (*len as f64);
                    for _ in 0..*len {
                        let q = *self
                            .qual
                            .get(read_i)
                            .ok_or(AlignmentError::QualIndexOutOfBounds)?;
                        read_i += 1;
                        tot_score -= log_likelihood_mismatch[q as usize];
                    }
                },
                UnifiedOp::Del(len) => tot_score += config.gap_open + config.gap_extend * (*len as f64),
                UnifiedOp::Mismatch(seq) => read_i += seq.len(), // mismatch against unmapped: no score
                UnifiedOp::Insertion(sv) => {
                    let len = sv.len();
                    tot_score += config.gap_open + config.gap_extend * (len as f64);
                    for _ in 0..len {
                        let q = *self
                            .qual
                            .get(read_i)
                            .ok_or(AlignmentError::QualIndexOutOfBounds)?;
                        read_i += 1;
                        tot_score -= log_likelihood_mismatch[q as usize];
                    }
                }
                UnifiedOp::Deletion(sv) => tot_score += config.gap_open + config.gap_extend * (sv.len() as f64),
                UnifiedOp::RefSkip(_) => {},
                UnifiedOp::Trans(score) | UnifiedOp::Translocate { score, ..} => tot_score += *score,
                UnifiedOp::Mate | UnifiedOp::MateSwitch {..} => {},
            }
        }
        Ok(tot_score)
    }
    fn score_one( &self, iter: UnifiedOpIterator<'a>) -> Result<f64, PrepareError> {
        let mut read_i = 0;
        let mut score = 0.0;

        for op in iter {
            score += self.score_op(&op, &mut read_i)?;
        }
        Ok(score)
    }

    /// Scores the alignment using the provided log likelihood (mis)match arrays.
    /// If unmapped read, it sums the mismatch scores for each base quality.
    ///
    /// Higher scores are better
    ///
    pub fn score(mut self) -> Result<f64> {
        let mut score = 0.0_f64;

        let aln_iter = AlignmentCompareIterator::new(self.iter1, self.iter2, self.qual);

        for pair in aln_iter {
            match pair? {
                AlnCmpOp::Equal => {},
                AlnCmpOp::UnEqual(op1, op2, q) => match (op1, op2) {
                    (ref x, ref y) if x.same_variant(y) => {}
                    (x, y) => {
                        score += x.score(q, &mut None) - y.score(q, &mut None);
                    }
                },
                AlnCmpOp::Left(op1, q) => {
                    score += op1.score(q, &mut None);
                },
                AlnCmpOp::Right(op2, q) => {
                    score -= op2.score(q, &mut None);
                }
            }
        }

        /*match (self.cigar1.take(), self.cigar2.take()) {
            (Some(cigar1), Some(cigar2)) => {
                let mut indel_gap1 = None;
                let mut indel_gap2 = None;

                let md_iter1 = self.md_iter1.take().ok_or(PrepareError::NoMdTag)?;
                let md_iter2 = self.md_iter2.take().ok_or(PrepareError::NoMdTag)?;
                let aln_iter = AlignmentCompareIterator::new(
                    cigar1.iter(),
                    md_iter1,
                    cigar2.iter(),
                    md_iter2,
                    self.qual,
                );

                for pair in aln_iter {
                    match pair? {
                        AlnCmpOp::Equal => {
                            indel_gap1 = None;
                            indel_gap2 = None;
                        }
                        AlnCmpOp::UnEqual(op1, op2, q) => match (op1, op2) {
                            (ref x, ref y) if x.same_variant(y) => {}
                            (AlignmentOp::Mismatch, AlignmentOp::SoftClip)
                            | (AlignmentOp::SoftClip, AlignmentOp::Mismatch) => {}
                            (x, y) => {
                                score += x.score(q, &mut indel_gap1) - y.score(q, &mut indel_gap2);
                            }
                        },
                        AlnCmpOp::Left(op1, q) => {
                            score += op1.score(q, &mut indel_gap1);
                        }
                        AlnCmpOp::Right(op2, q) => {
                            score -= op2.score(q, &mut indel_gap2);
                        }
                    }
                }
            }
            (Some(cigar1), None) => {
                let md_iter1 = self.md_iter1.take().ok_or(PrepareError::NoMdTag)?;
                score += self.score_aln_op_against_mm(cigar1, md_iter1)?;
            }
            (None, Some(cigar2)) => {
                let md_iter2 = self.md_iter2.take().ok_or(PrepareError::NoMdTag)?;
                score -= self.score_aln_op_against_mm(cigar2, md_iter2)?;
            }
            (None, None) => {}
        } */
        Ok(score)
    }
}

pub struct PreparedAlignmentPairIter<'a> {
    records1: std::slice::Iter<'a, Record>,
    records2: std::slice::Iter<'a, Record>,
}

impl<'a> PreparedAlignmentPairIter<'a> {
    #[must_use]
    pub fn new(alns1: &'a [Record], alns2: &'a [Record]) -> Self {
        Self {
            records1: alns1.iter(),
            records2: alns2.iter(),
        }
    }
}

impl<'a> Iterator for PreparedAlignmentPairIter<'a> {
    type Item = Result<PreparedAlignmentPair<'a>, PrepareError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut read_host: Option<&Record> = None;
        let mut read_graft: Option<&Record> = None;
        loop {
            if read_host.is_none() {
                read_host = match self.records1.next() {
                    Some(r) if r.is_secondary() => continue,
                    Some(r) => Some(r),
                    None => return None,
                };
            }
            if read_graft.is_none() {
                read_graft = match self.records2.next() {
                    Some(r) if r.is_secondary() => continue,
                    Some(r) => Some(r),
                    None => return None,
                };
            }
            let read_host = read_host.take().unwrap();
            let qual = read_host.qual();
            let read_graft = read_graft.take().unwrap();
            assert_eq!(qual, read_graft.qual());

            if read_host.is_unmapped() {
                return Some(Ok(PreparedAlignmentPair {
                    cigar1: None,
                    cigar2: read_graft.is_unmapped().then_some(read_graft.cigar()),
                    md_iter1: None,
                    md_iter2: None,
                    qual,
                }));
            }
            if read_graft.is_unmapped() {
                return Some(Ok(PreparedAlignmentPair {
                    cigar1: Some(read_host.cigar()),
                    cigar2: None,
                    md_iter1: None,
                    md_iter2: None,
                    qual,
                }));
            }

            let md_iter1 = match read_host.aux(b"MD") {
                Ok(Aux::String(md)) => Some(MdOpIterator::new(md)),
                Ok(_) => return Some(Err(PrepareError::NoMdTag)),
                Err(e) => return Some(Err(PrepareError::AuxError(e.to_string()))),
            };
            let md_iter2 = match read_graft.aux(b"MD") {
                Ok(Aux::String(md)) => Some(MdOpIterator::new(md)),
                Ok(_) => return Some(Err(PrepareError::NoMdTag)),
                Err(e) => return Some(Err(PrepareError::AuxError(e.to_string()))),
            };
            return Some(Ok(PreparedAlignmentPair {
                cigar1: Some(read_host.cigar()),
                cigar2: Some(read_graft.cigar()),
                md_iter1,
                md_iter2,
                qual,
            }));
        }
    }
}
