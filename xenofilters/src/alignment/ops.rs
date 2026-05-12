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
        let md_iter = MdOpIterator::new(rec)?;
        let cigar_iter = rec.cigar().to_vec().into_iter();
        Ok(UnifiedOpIterator {
            cigar_iter,
            md_iter,
            next_op: None,
            next_md_op: None,
            next_cig: None,
        })
    }
    pub(crate) fn empty() -> Self {
        UnifiedOpIterator {
            cigar_iter: vec![].into_iter(),
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
        let cig_ct: usize = self.cigar_iter.as_slice().len();
        cig_ct == 1 && self.md_iter.is_single_operation()
    }
}

impl<'a> Iterator for UnifiedOpIterator<'a> {
    type Item = Result<UnifiedOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_op.is_some() {
            return Some(Ok(self.next_op.take().unwrap()));
        }
        let next_cig = if self.next_cig.is_some() {
            self.next_cig.take()
        } else {
            self.cigar_iter.next()
        };
        //#[cfg(test)]
        //eprintln!("Next CIGAR: {:?}", next_cig);
        match next_cig {
            Some(Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len)) => {
                let next_md_op = if self.next_md_op.is_some() {
                    self.next_md_op.take().map(Ok)
                } else {
                    self.md_iter.next()
                };
                //#[cfg(test)]
                //eprintln!("Next MD op: {:?}", next_md_op);
                if let Some(Ok(md_op)) = next_md_op {
                    //#[cfg(test)]
                    //eprintln!("Matching CIGAR {:?} with MD {:?}", next_cig, &md_op);
                    return Some(self.match_md_op(md_op, len));
                };
                #[cfg(test)]
                eprintln!("No MD op to match CIGAR {:?}", next_cig);
                Some(Err(AlignmentError::MdCigarMismatch))
            }
            Some(Cigar::SoftClip(len)) => {
                //#[cfg(test)]
                //eprintln!("Handling CIGAR SoftClip {:?}", len);
                // a sofclip has no matching MD operation
                Some(Ok(UnifiedOp::Mis(len)))
            }
            Some(Cigar::Del(len)) => {
                let next_md_op = if self.next_md_op.is_some() {
                    self.next_md_op.take().map(Ok)
                } else {
                    self.md_iter.next()
                };
                //#[cfg(test)]
                //eprintln!("Next MD op: {:?}", next_md_op);
                match next_md_op {
                    Some(Ok(MdOp::Deletion(d))) if d.len() as u32 == len => {
                        //#[cfg(test)]
                        //eprintln!("Matching CIGAR Deletion {:?} with MD Deletion {:?}", len, d);
                        Some(Ok(UnifiedOp::Del(len)))
                    }
                    _ => {
                        #[cfg(test)]
                        eprintln!("CIGAR Deletion {:?} has no matching MD Deletion", len);
                        Some(Err(AlignmentError::MdCigarMismatch))
                    }
                }
            }
            Some(Cigar::HardClip(_) | Cigar::Pad(_)) => {
                //#[cfg(test)]
                //eprintln!("Skipping CIGAR HardClip/Pad {:?}", next_cig);
                self.next()
            }
            Some(x) => {
                //#[cfg(test)]
                //eprintln!("Handling Refskip/Ins CIGAR {:?}", x);
                Some(UnifiedOp::try_from(x))
            } // RefSkip, Ins
            None => {
                //#[cfg(test)]
                //eprintln!("CIGAR operations exhausted, checking MD iterator");
                let next_md_op = if self.next_md_op.is_some() {
                    self.next_md_op.take().map(Ok)
                } else {
                    self.md_iter.next()
                };
                //#[cfg(test)]
                //eprintln!("Next MD op: {:?}", next_md_op);
                if let Some(Ok(md_op)) = next_md_op {
                    #[cfg(test)]
                    eprintln!("Excess MD operation after CIGAR end: {:?}", md_op);
                    match md_op {
                        MdOp::Match(_) | MdOp::Mismatch(_) => {
                            Some(Err(AlignmentError::MdCigarMismatch))
                        }
                        MdOp::Deletion(_) => Some(Err(AlignmentError::MissingMdDeletion)),
                    }
                } else {
                    None
                }
            }
        }
    }
}

#[cfg(test)]
pub mod tests;
