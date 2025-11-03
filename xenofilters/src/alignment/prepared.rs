use anyhow::Result;
use rust_htslib::bam::record::{Aux, Cigar, CigarStringView, Record};

use super::{
    AlignmentCompareIterator, AlignmentError, AlignmentIterator, AlnCmpOp, MdOp, MdOpIterator,
    PrepareError,
};
use crate::{AlignmentOp, MAX_Q};

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
    cigar1: Option<CigarStringView>,
    cigar2: Option<CigarStringView>,
    md_iter1: Option<MdOpIterator<'a>>,
    md_iter2: Option<MdOpIterator<'a>>,
    qual: &'a [u8],
}

impl<'a> PreparedAlignmentPair<'a> {
    #[must_use]
    pub fn is_perfect_match(&self, first: bool) -> bool {
        if first {
            if let Some(cigar) = &self.cigar1
                && cigar
                    .iter()
                    .all(|c| matches!(c, Cigar::Match(_) | Cigar::Equal(_)))
                && let Some(md_iter) = &self.md_iter1
            {
                let mut md_iter = md_iter.clone();
                return matches!(md_iter.next(), Some(MdOp::Match(_))) && md_iter.next().is_none();
            }
        } else if let Some(cigar) = &self.cigar2
            && cigar
                .iter()
                .all(|c| matches!(c, Cigar::Match(_) | Cigar::Equal(_)))
            && let Some(md_iter) = &self.md_iter2
        {
            let mut md_iter = md_iter.clone();
            return matches!(md_iter.next(), Some(MdOp::Match(_))) && md_iter.next().is_none();
        }
        false
    }

    fn score_aln_op_against_mm(
        &self,
        cigar: CigarStringView,
        md_iter: MdOpIterator<'a>,
        log_likelihood_mismatch: &[f64; MAX_Q + 2],
    ) -> Result<f64, PrepareError> {
        let mut read_i = 0;
        let mut indel_gap = None;
        let mut score = 0.0;

        for op in AlignmentIterator::new(cigar.take().iter(), md_iter) {
            let q = *self.qual.get(read_i).ok_or(AlignmentError::QualIndexOutOfBounds)?;
            let op = op?;
            if !matches!(op, AlignmentOp::Deletion | AlignmentOp::RefSkip(_)) {
                read_i += 1;
            }
            score += op.score(q, &mut indel_gap, log_likelihood_mismatch)
                - log_likelihood_mismatch[q as usize];
        }
        Ok(score)
    }

    /// Scores the alignment using the provided log likelihood (mis)match arrays.
    /// If unmapped read, it sums the mismatch scores for each base quality.
    ///
    /// Higher scores are better
    ///
    pub fn score(mut self, log_likelihood_mismatch: &[f64; MAX_Q + 2]) -> Result<f64> {
        let mut score = 0.0_f64;

        match (self.cigar1.take(), self.cigar2.take()) {
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
                                score += x.score(q, &mut indel_gap1, log_likelihood_mismatch)
                                    - y.score(q, &mut indel_gap2, log_likelihood_mismatch);
                            }
                        },
                        AlnCmpOp::Left(op1, q) => {
                            score += op1.score(q, &mut indel_gap1, log_likelihood_mismatch);
                        }
                        AlnCmpOp::Right(op2, q) => {
                            score -= op2.score(q, &mut indel_gap2, log_likelihood_mismatch);
                        }
                    }
                }
            }
            (Some(cigar1), None) => {
                let md_iter1 = self.md_iter1.take().ok_or(PrepareError::NoMdTag)?;
                score += self.score_aln_op_against_mm(cigar1, md_iter1, log_likelihood_mismatch)?;
            }
            (None, Some(cigar2)) => {
                let md_iter2 = self.md_iter2.take().ok_or(PrepareError::NoMdTag)?;
                score -= self.score_aln_op_against_mm(cigar2, md_iter2, log_likelihood_mismatch)?;
            }
            (None, None) => {}
        }
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
        let mut r1: Option<&Record> = None;
        let mut r2: Option<&Record> = None;
        loop {
            if r1.is_none() {
                r1 = match self.records1.next() {
                    Some(r) if r.is_secondary() => continue,
                    Some(r) => Some(r),
                    None => return None,
                };
            }
            if r2.is_none() {
                r2 = match self.records2.next() {
                    Some(r) if r.is_secondary() => continue,
                    Some(r) => Some(r),
                    None => return None,
                };
            }
            let r1 = r1.take().unwrap();
            let qual = r1.qual();
            let r2 = r2.take().unwrap();
            assert_eq!(qual, r2.qual());

            if r1.is_unmapped() {
                return Some(Ok(PreparedAlignmentPair {
                    cigar1: None,
                    cigar2: r2.is_unmapped().then_some(r2.cigar()),
                    md_iter1: None,
                    md_iter2: None,
                    qual,
                }));
            }
            if r2.is_unmapped() {
                return Some(Ok(PreparedAlignmentPair {
                    cigar1: Some(r1.cigar()),
                    cigar2: None,
                    md_iter1: None,
                    md_iter2: None,
                    qual,
                }));
            }

            let md_iter1 = match r1.aux(b"MD") {
                Ok(Aux::String(md)) => Some(MdOpIterator::new(md)),
                Ok(_) => return Some(Err(PrepareError::NoMdTag)),
                Err(e) => return Some(Err(PrepareError::AuxError(e.to_string()))),
            };
            let md_iter2 = match r2.aux(b"MD") {
                Ok(Aux::String(md)) => Some(MdOpIterator::new(md)),
                Ok(_) => return Some(Err(PrepareError::NoMdTag)),
                Err(e) => return Some(Err(PrepareError::AuxError(e.to_string()))),
            };
            return Some(Ok(PreparedAlignmentPair {
                cigar1: Some(r1.cigar()),
                cigar2: Some(r2.cigar()),
                md_iter1,
                md_iter2,
                qual,
            }));
        }
    }
}
