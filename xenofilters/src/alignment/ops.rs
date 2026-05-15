use crate::alignment::PrepareError;
use crate::alignment::{AlignmentError, MdOp, MdOpIterator};
use noodles::bam::record::Record;
use std::cmp::min;
use noodles::sam::alignment::record::cigar::{op::Kind, Op};
use std::io;


#[derive(Debug, Clone, PartialEq)]
pub(crate) enum UnifiedOp {
    /// A pure reference match
    Match(usize),

    /// A refskip from CIGAR
    RefSkip(usize),

    // simple variants, not variant-aware
    Mis(usize),
    Ins(usize),
    Del(usize),
    Relocate {
        pos: usize,
        penalty_score: f64,
    },
}
impl TryFrom<(Kind, usize)> for UnifiedOp {
    type Error = AlignmentError;

    fn try_from(t: (Kind, usize)) -> Result<Self, Self::Error> {
        match t.0 {
            Kind::Insertion => Ok(UnifiedOp::Ins(t.1)),
            Kind::Skip => Ok(UnifiedOp::RefSkip(t.1)),
            _ => Err(AlignmentError::UnImplemented),
        }
    }
}

#[cfg_attr(test, derive(Debug))]
pub(crate) struct UnifiedOpIterator<'a> {
    cigar: Option<Box<dyn Iterator<Item = io::Result<Op>> + 'a>>,
    cigar_len: Option<usize>,
    md_iter: MdOpIterator<'a>,
    next_md_op: Option<MdOp>,
    next_cig: Option<Op>,
    next_op: Option<UnifiedOp>,
}

impl<'a> UnifiedOpIterator<'a> {
    pub(crate) fn new(rec: &'a Record) -> Result<Self, PrepareError> {
        let cig = rec.cigar();
        let cigar_len = cig.len();
        let cigar_vec: Vec<_> = cig.iter().collect();
        let cigar: Box<dyn Iterator<Item = io::Result<Op>> + 'a> = Box::new(cigar_vec.into_iter());
        Ok(Self {
            cigar: Some(cigar),
            cigar_len: Some(cigar_len),
            md_iter: MdOpIterator::new(rec)?,
            next_op: None,
            next_md_op: None,
            next_cig: None,
        })
    }

    pub(crate) fn empty() -> Self {
        Self {
            cigar: None,
            cigar_len: None,
            md_iter: MdOpIterator::empty(),
            next_op: None,
            next_md_op: None,
            next_cig: None,
        }
    }

    fn match_md_op(&mut self, md_op: MdOp, cig_len: usize) -> Result<UnifiedOp, AlignmentError> {
        match md_op {
            MdOp::Match(md_len) => {
                let op_len = min(cig_len, md_len);

                if md_len > op_len {
                    self.next_md_op = Some(MdOp::Match(md_len - op_len));
                } else if cig_len > op_len {
                    self.next_cig = Some(Op::new(Kind::Match, cig_len - op_len));
                }
                Ok(UnifiedOp::Match(op_len))
            }
            MdOp::Mismatch(md_len) => {
                if md_len != 1 || cig_len != 1 {
                    self.next_cig = Some(Op::new(Kind::Match, cig_len - 1));
                }
                Ok(UnifiedOp::Mis(1))
            }
            MdOp::Deletion(_) => Err(AlignmentError::MdCigarMismatch),
        }
    }

    pub(crate) fn is_single_match(&self) -> bool {
        match self.cigar_len {
            Some(1) => !self.md_iter.is_single_operation(),
            _ => false,
        }
    }
}

impl<'a> Iterator for UnifiedOpIterator<'a> {
    type Item = Result<UnifiedOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(op) = self.next_op.take() {
                return Some(Ok(op));
            }

            let c = match self.next_cig.take() {
                Some(c) => c,
                None if self.cigar.is_none() => return None,
                None => {
                    match self.cigar.as_mut().and_then(|c| c.next()) {
                        Some(Ok(op)) => op,
                        _ => {
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
            };

            match c.kind() {
                Kind::Match| Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let next_md = self.next_md_op.take().map(Ok).or_else(|| self.md_iter.next());
                    if let Some(Ok(md_op)) = next_md {
                        return Some(self.match_md_op(md_op, c.len()));
                    }
                    return Some(Err(AlignmentError::MdCigarMismatch));
                }
                Kind::SoftClip => return Some(Ok(UnifiedOp::Mis(c.len()))),
                Kind::Deletion => {
                    let next_md = self.next_md_op.take().map(Ok).or_else(|| self.md_iter.next());
                    match next_md {
                        Some(Ok(MdOp::Deletion(d))) if d.len() == c.len() => {
                            return Some(Ok(UnifiedOp::Del(c.len())));
                        }
                        _ => return Some(Err(AlignmentError::MdCigarMismatch)),
                    }
                }
                Kind::HardClip | Kind::Pad => continue,
                x => return Some(UnifiedOp::try_from((x, c.len()))),
            }
        }
    }
}

#[cfg(test)]
pub(crate) mod tests;
