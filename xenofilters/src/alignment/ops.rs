use crate::alignment::PrepareError;
use crate::alignment::{AlignmentError, MdOp, MdOpIterator};
use rust_htslib::bam::record::{Cigar, Record};
use std::cmp::min;

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum UnifiedOp {
    /// A pure reference match
    Match(u32),

    /// A refskip from CIGAR
    RefSkip(u32),

    // simple variants, not variant-aware
    Mis(u32),
    Ins(u32),
    Del(u32),
    Relocate {
        pos: i64,
        penalty_score: f64,
    },
}
impl TryFrom<Cigar> for UnifiedOp {
    type Error = AlignmentError;

    fn try_from(cigar: Cigar) -> Result<Self, Self::Error> {
        match cigar {
            Cigar::Ins(len) => Ok(UnifiedOp::Ins(len)),
            Cigar::RefSkip(len) => Ok(UnifiedOp::RefSkip(len)),
            _ => Err(AlignmentError::UnImplemented),
        }
    }
}

#[cfg_attr(test, derive(Debug))]
pub(crate) struct UnifiedOpIterator<'a> {
    cigar_iter: std::vec::IntoIter<Cigar>,
    md_iter: MdOpIterator<'a>,
    next_md_op: Option<MdOp>,
    next_cig: Option<Cigar>,
    next_op: Option<UnifiedOp>,
}

impl<'a> UnifiedOpIterator<'a> {
    pub(crate) fn new(rec: &'a Record) -> Result<Self, PrepareError> {
        Ok(Self {
            cigar_iter: rec.cigar().to_vec().into_iter(),
            md_iter: MdOpIterator::new(rec)?,
            next_op: None,
            next_md_op: None,
            next_cig: None,
        })
    }

    pub(crate) fn empty() -> Self {
        Self {
            cigar_iter: Vec::new().into_iter(),
            md_iter: MdOpIterator::empty(),
            next_op: None,
            next_md_op: None,
            next_cig: None,
        }
    }

    fn match_md_op(&mut self, md_op: MdOp, cig_len: u32) -> Result<UnifiedOp, AlignmentError> {
        match md_op {
            MdOp::Match(md_len) => {
                let op_len = min(cig_len, md_len);

                if md_len > op_len {
                    self.next_md_op = Some(MdOp::Match(md_len - op_len));
                } else if cig_len > op_len {
                    self.next_cig = Some(Cigar::Match(cig_len - op_len));
                }
                Ok(UnifiedOp::Match(op_len))
            }
            MdOp::Mismatch(md_len) => {
                if md_len as u32 != 1 || cig_len != 1 {
                    self.next_cig = Some(Cigar::Match(cig_len - 1));
                }
                Ok(UnifiedOp::Mis(1))
            }
            MdOp::Deletion(_) => Err(AlignmentError::MdCigarMismatch),
        }
    }

    pub(crate) fn is_single_match(&self) -> bool {
        self.cigar_iter.len() == 1 && self.md_iter.is_single_operation()
    }
}

impl<'a> Iterator for UnifiedOpIterator<'a> {
    type Item = Result<UnifiedOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(op) = self.next_op.take() {
                return Some(Ok(op));
            }

            let next_cig = self.next_cig.take().or_else(|| self.cigar_iter.next());

            match next_cig {
                Some(Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len)) => {
                    let next_md = self.next_md_op.take().map(Ok).or_else(|| self.md_iter.next());
                    if let Some(Ok(md_op)) = next_md {
                        return Some(self.match_md_op(md_op, len));
                    }
                    return Some(Err(AlignmentError::MdCigarMismatch));
                }
                Some(Cigar::SoftClip(len)) => return Some(Ok(UnifiedOp::Mis(len))),
                Some(Cigar::Del(len)) => {
                    let next_md = self.next_md_op.take().map(Ok).or_else(|| self.md_iter.next());
                    match next_md {
                        Some(Ok(MdOp::Deletion(d))) if d.len() as u32 == len => {
                            return Some(Ok(UnifiedOp::Del(len)));
                        }
                        _ => return Some(Err(AlignmentError::MdCigarMismatch)),
                    }
                }
                Some(Cigar::HardClip(_) | Cigar::Pad(_)) => {
                    continue; // Replaces recursive call to avoid stack depth and overhead
                }
                Some(x) => return Some(UnifiedOp::try_from(x)),
                None => {
                    let next_md = self.next_md_op.take().map(Ok).or_else(|| self.md_iter.next());
                    if let Some(Ok(md_op)) = next_md {
                        match md_op {
                            MdOp::Match(_) | MdOp::Mismatch(_) => {
                                return Some(Err(AlignmentError::MdCigarMismatch))
                            }
                            MdOp::Deletion(_) => return Some(Err(AlignmentError::MissingMdDeletion)),
                        }
                    }
                    return None;
                }
            }
        }
    }
}

#[cfg(test)]
pub mod tests;
