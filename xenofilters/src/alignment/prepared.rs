use rust_htslib::bam::record::{Aux, Cigar, CigarStringView, Record};
use anyhow::Result;
use crate::{MAX_Q, LOG_LIKELIHOOD_MATCH, AlignmentOp};
use super::{MdOp, parse_md, PrepareError, AlignmentIterator};

pub struct PreparedAlignment<'a> {
    cigar: Option<CigarStringView>,
    md_ops: Vec<MdOp>,
    qual: &'a [u8],
}

impl<'a> PreparedAlignment<'a> {
    pub fn is_perfect_match(&self) -> bool {
        if let Some(cigar) = &self.cigar {
            cigar.iter().all(|c| matches!(c, Cigar::Match(_) | Cigar::Equal(_)))
                && self.md_ops.iter().all(|m| matches!(m, MdOp::Match(_)))
        } else {
            false
        }
    }

    pub fn score(mut self, log_likelihood_mismatch: &[f64; MAX_Q + 2]) -> Result<f64> {
        if let Some(cigar) = self.cigar.take() {
            let mut score = 0.0_f64;
            let gap_open = log_likelihood_mismatch[MAX_Q];
            let gap_ext = log_likelihood_mismatch[MAX_Q + 1];

            // --- State for Affine Gap Penalties ---
            let mut indel_gap = None; // None, Some(true)=insertion, Some(false)=deletion
                                  //
            let aln_iter = AlignmentIterator::new(cigar, self.md_ops.clone(), self.qual);

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
                            score += gap_open; // Add open penalty only for the first base in the gap.
                            indel_gap = Some(true)
                        }
                        score += gap_ext; // Add extend penalty for every base.
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
                // An unmapped read is a valid item to yield for potential scoring.
                return Some(Ok(PreparedAlignment {
                    cigar: None,
                    md_ops: vec![],
                    qual: r.qual(),
                }));
            }

            let md_str = match r.aux(b"MD") {
                Ok(Aux::String(md)) => md,
                Ok(_) => return Some(Err(PrepareError::NoMdTag)),
                Err(e) => return Some(Err(PrepareError::AuxError(e.to_string()))),
            };

            let md_ops = parse_md(md_str);
            let cigar = Some(r.cigar());

            return Some(Ok(PreparedAlignment {
                cigar,
                md_ops,
                qual: r.qual(),
            }));
        }
    }
}
