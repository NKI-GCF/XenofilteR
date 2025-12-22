use crate::alignment::PrepareError;
use crate::alignment::{AlignmentError, MdOp, MdOpIterator};
use rust_htslib::bam::record::{Cigar, Record};
use std::cmp::min;

// TODO: for paired-end introduce another operation for read switching

#[derive(Debug, Clone)]
pub enum UnifiedOp {
    /// A pure reference match
    Match(u32),

    /// A refskip from CIGAR
    RefSkip(u32),

    // simple variants, not variant-aware
    Mis(u32),
    Ins(u32),
    Del(u32),
    Trans(f64),
    /* TODO: variant-aware operations
    /// A mismatch or a softclip, with sequence
    Mismatch(SmallVec<[u8; 1]>), // The reference base

    /// An insertion with bases
    Insertion(SmallVec<[u8; 1]>),

    /// A deletion with ref bases
    Deletion(SmallVec<[u8; 1]>),

    /// A translocation (stitched CIGAR)
    Translocate {
        new_tid: i32,
        new_pos: i64,
        new_strand: bool,
        score: f64,
    },*/
}

impl UnifiedOp {
    #[must_use]
    pub fn len(&self) -> u32 {
        match self {
            UnifiedOp::Match(len)
            | UnifiedOp::Mis(len)
            | UnifiedOp::Ins(len)
            | UnifiedOp::Del(len)
            | UnifiedOp::RefSkip(len) => *len,
            //UnifiedOp::Mismatch(seq) | UnifiedOp::Insertion(seq) | UnifiedOp::Deletion(seq) => seq.len(),
            _ => 0,
        }
    }
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

pub struct UnifiedOpIterator<'a> {
    cigar_iter: std::vec::IntoIter<Cigar>,
    md_iter: MdOpIterator<'a>,
    next_op: Option<UnifiedOp>,
    is_rev: bool,
}

impl<'a> UnifiedOpIterator<'a> {
    pub fn new(rec: &'a Record) -> Result<Self, PrepareError> {
        let md_iter = MdOpIterator::new(rec)?;
        let cigar_iter = rec.cigar().to_vec().into_iter();
        Ok(UnifiedOpIterator {
            cigar_iter,
            md_iter,
            next_op: None,
            is_rev: rec.is_reverse(),
        })
    }
    pub fn empty(is_rev: bool) -> Self {
        UnifiedOpIterator {
            cigar_iter: vec![].into_iter(),
            md_iter: MdOpIterator::empty(),
            next_op: None,
            is_rev,
        }
    }

    pub fn peek(&'a mut self) -> Option<&'a UnifiedOp> {
        if self.next_op.is_none() {
            match self.next() {
                Some(Ok(op)) => self.next_op = Some(op),
                Some(Err(_)) => return None,
                None => return None,
            }
        }
        self.next_op.as_ref()
    }

    fn process_match(&mut self, op: UnifiedOp, len: u32) -> UnifiedOp {
        let md_len = op.len();
        let op_len = min(len, md_len);

        if md_len > op_len {
            self.next_op = Some(UnifiedOp::Match(md_len - op_len));
        }
        UnifiedOp::Match(op_len)
    }
    fn match_md_op(&mut self, md_op: MdOp, len: u32) -> Result<UnifiedOp, AlignmentError> {
        match md_op {
            MdOp::Match(md_len) => {
                let op_len = std::cmp::min(len, md_len);

                if md_len > op_len {
                    self.next_op = Some(UnifiedOp::Match(md_len - op_len));
                }
                Ok(UnifiedOp::Match(op_len))
            }
            MdOp::Mismatch(_) => Ok(UnifiedOp::Mis(1)),
            MdOp::Deletion(_) => Err(AlignmentError::MdCigarMismatch),
        }
    }
}

impl<'a> Iterator for UnifiedOpIterator<'a> {
    type Item = Result<UnifiedOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(op) = self.next_op.take() {
            return Some(Ok(op));
        }
        let next_cig = if self.is_rev {
            self.cigar_iter.next_back()
        } else {
            self.cigar_iter.next()
        };
        match next_cig? {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                if let Some(op) = self.next_op.take() {
                    return Some(Ok(self.process_match(op, len)));
                }
                let next_md_op = if self.is_rev {
                    self.md_iter.next_back()
                } else {
                    self.md_iter.next()
                };
                if let Some(Ok(md_op)) = next_md_op {
                    return Some(self.match_md_op(md_op, len));
                };
                Some(Err(AlignmentError::MdCigarMismatch))
            }
            Cigar::SoftClip(len) => Some(Ok(UnifiedOp::Mis(len))),
            Cigar::Del(len) => match self.md_iter.next() {
                Some(Ok(MdOp::Deletion(d))) if d.len() as u32 == len => {
                    Some(Ok(UnifiedOp::Del(len)))
                }
                _ => Some(Err(AlignmentError::MdCigarMismatch)),
            },
            Cigar::HardClip(_) | Cigar::Pad(_) => self.next(),
            x => Some(UnifiedOp::try_from(x)), // RefSkip, Ins
        }
    }
}
