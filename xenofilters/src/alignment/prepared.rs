use rust_htslib::bam::record::{Aux, Cigar, CigarStringView, Record};
use anyhow::Result;
use crate::{MAX_Q, LOG_LIKELIHOOD_MATCH, AlignmentOp};
use super::{MdOp, MdOpIterator, PrepareError, AlignmentIterator};

pub struct PreparedAlignment<'a> {
    cigar: Option<CigarStringView>,
    md_iter: Option<MdOpIterator<'a>>,
    qual: &'a [u8],
}

impl<'a> PreparedAlignment<'a> {
    pub fn is_perfect_match(&mut self) -> bool {
        if let Some(cigar) = &self.cigar
            && cigar.iter().all(|c| matches!(c, Cigar::Match(_) | Cigar::Equal(_)))
                && let Some(md_iter) = &mut self.md_iter {
            return md_iter.all(|m| matches!(m, MdOp::Match(_)));
        }
        false
    }

    /// Scores the alignment using the provided log likelihood mismatch array.
    /// If no CIGAR is present (unmapped read), it sums the mismatch scores for each base quality.
    /// If a CIGAR is present, it uses the CIGAR and MD tag to compute the score.
    /// The log_likelihood_mismatch array must have length MAX_Q + 2, where the last two entries
    /// are used for gap open and gap extension penalties.
    ///
    /// Higher scores are better, with perfect matches scoring f64::INFINITY.
    ///
    pub fn score(mut self, log_likelihood_mismatch: &[f64; MAX_Q + 2]) -> Result<f64> {
        if let Some(cigar) = self.cigar.take() {
            let mut score = 0.0_f64;
            let gap_open = log_likelihood_mismatch[MAX_Q];
            let gap_ext = log_likelihood_mismatch[MAX_Q + 1];

            let mut indel_gap = None; // Some(true)=insertion, Some(false)=deletion

            let md_iter = self.md_iter.take().ok_or(PrepareError::NoMdTag)?;
            let aln_iter = AlignmentIterator::new(cigar, md_iter, self.qual);

            for op in aln_iter {
                match op? {
                    AlignmentOp::Match(q) => {
                        score += LOG_LIKELIHOOD_MATCH[q as usize];
                        indel_gap = None;
                    },
                    AlignmentOp::Mismatch(q) | AlignmentOp::SoftClip(q) => {
                        score += log_likelihood_mismatch[q as usize];
                        indel_gap = None;
                    },
                    AlignmentOp::Insertion(_q) => {
                        if indel_gap != Some(true) {
                            score += gap_open;
                            indel_gap = Some(true)
                        }
                        score += gap_ext;
                    },
                    AlignmentOp::Deletion => {
                        if indel_gap == Some(false) {
                            score += gap_open;
                            indel_gap = Some(false)
                        }
                        score += gap_ext;
                    },
                }
            }
            Ok(score)
        } else {
            Ok(self.qual.iter().map(|&q| log_likelihood_mismatch[q as usize]).sum())
        }
    }
}

pub struct PreparedAlignmentIter<'a> {
    records: std::slice::Iter<'a, Record>,
}

impl<'a> PreparedAlignmentIter<'a> {
    pub fn new(alns: &'a [Record]) -> Self {
        Self { records: alns.iter() }
    }
}

impl<'a> Iterator for PreparedAlignmentIter<'a> {
    type Item = Result<PreparedAlignment<'a>, PrepareError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let r = self.records.next()?;

            if r.is_secondary() {
                continue;
            }

            if r.is_unmapped() {
                return Some(Ok(PreparedAlignment {
                    cigar: None,
                    md_iter: None,
                    qual: r.qual(),
                }));
            }

            let md_iter = match r.aux(b"MD") {
                Ok(Aux::String(md)) => Some(MdOpIterator::new(md)),
                Ok(_) => return Some(Err(PrepareError::NoMdTag)),
                Err(e) => return Some(Err(PrepareError::AuxError(e.to_string()))),
            };

            return Some(Ok(PreparedAlignment {
                cigar: Some(r.cigar()),
                md_iter,
                qual: r.qual(),
            }));
        }
    }
}
