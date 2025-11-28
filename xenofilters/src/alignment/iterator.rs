use super::{AlignmentError, AlnCmpOp};
use anyhow::Result;

use crate::alignment::{UnifiedOpIterator, UnifiedOp};

pub struct AlignmentCompareIterator<'a> {
    iter1: UnifiedOpIterator<'a>,
    iter2: UnifiedOpIterator<'a>,
    qual: &'a [u8],
    read_i: usize,
}

impl<'a> AlignmentCompareIterator<'a> {
    #[must_use]
    pub fn new(
        iter1: UnifiedOpIterator<'a>,
        iter2: UnifiedOpIterator<'a>,
        qual: &'a [u8],
    ) -> Self {
        AlignmentCompareIterator {
            iter1,
            iter2,
            qual,
            read_i: 0,
        }
    }
    fn get_qual(&mut self) -> Result<u8, AlignmentError> {
        if let Some(&q) = self.qual.get(self.read_i) {
            self.read_i += 1;
            Ok(q)
        } else {
            Err(AlignmentError::QualIndexOutOfBounds)
        }
    }
    /// Handles identical operations that consume read (Ins, SoftClip).
    fn handle_equal_ops_read(&'a mut self, len1: u32, len2: u32, op_fn: fn(u32) -> UnifiedOp) -> AlnCmpOp {
        let len = len1.min(len2);
        self.read_i += len as usize;

        self.set_cigar_ops(op_fn(len1 - len), op_fn(len2 - len));
        AlnCmpOp::Equal
    }

    fn set_cigar_ops(&'a mut self, op1: UnifiedOp, op2: UnifiedOp) {
        self.iter1.set_op(op1);
        self.iter2.set_op(op2);
    }

    /*fn next_md_ops(&mut self) -> (Option<MdOp>, Option<MdOp>) {
        (self.iter1.next_md_op(), self.iter2.next_md_op())
    }
    fn next_cigar_ops(&mut self) -> (Option<UnifiedOp>, Option<UnifiedOp>) {
        (self.iter1.next_cigar_op(), self.iter2.next_cigar_op())
    }

    /// Handles identical (Del, Del) operations.
    fn handle_equal_ops_del(&mut self, len1: u32, len2: u32) -> AlnCmpOp {
        let len = len1.min(len2);
        if let Some(p) = self.iter1.peek_mut() {
            *p = UnifiedOp::Deletion(len1 - len);
        }
        if let Some(p) = self.iter2.peek_mut() {
            *p = UnifiedOp::Deletion(len2 - len);
        }
        AlnCmpOp::Equal
    }

    /// Handles the "slow path" where only iter1 has an operation.
    fn handle_differing_left(&mut self, cigar1: UnifiedOp) -> Option<Result<AlnCmpOp, AlignmentError>> {
        self.iter1
            .get_next_op(cigar1)
            .map(|res| res.and_then(|op| self.get_qual().map(|q| AlnCmpOp::Left(op, q))))
            .or(self
                .iter1
                .next_cigar_op()
                .and_then(|cigar| self.handle_differing_left(cigar)))
    }

    /// Handles the "slow path" where only iter2 has an operation.
    fn handle_differing_right(
        &mut self,
        cigar2: UnifiedOp,
    ) -> Option<Result<AlnCmpOp, AlignmentError>> {
        self.iter2
            .get_next_op(cigar2)
            .map(|res| res.and_then(|op| self.get_qual().map(|q| AlnCmpOp::Right(op, q))))
            .or(self
                .iter2
                .next_cigar_op()
                .and_then(|cigar| self.handle_differing_right(cigar)))
    }

    /// Handles the "slow path" where both iters have differing operations.
    fn handle_differing_both(
        &mut self,
        cigar1: UnifiedOp,
        cigar2: UnifiedOp,
    ) -> Option<Result<AlnCmpOp, AlignmentError>> {
        self.iter1
            .get_next_op(cigar1)
            .and_then(|res1| {
                res1.ok().and_then(|op1| {
                    self.iter2
                        .get_next_op(cigar2)
                        .map(|res2| {
                            res2.and_then(|op2| {
                                self.get_qual().map(|q| AlnCmpOp::UnEqual(op1, op2, q))
                            })
                        })
                        .or(self
                            .iter2
                            .next_cigar_op()
                            .map_or(self.handle_differing_left(cigar1), |cigar| {
                                self.handle_differing_both(cigar1, cigar)
                            }))
                })
            })
            .or(self
                .iter1
                .next_cigar_op()
                .map_or(self.handle_differing_right(cigar2), |cigar| {
                    self.handle_differing_both(cigar, cigar2)
                }))
    }*/
}

impl<'a> Iterator for AlignmentCompareIterator<'a> {
    type Item = Result<AlnCmpOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        match (self.iter1.peek(), self.iter2.peek()) {
            (Some(UnifiedOp::Match(l1)), Some(UnifiedOp::Match(l2))) => {
                let skip = *l1.min(l2);
                self.read_i += skip as usize;
                self.iter1.set_op(UnifiedOp::Match(l1 - skip));
                self.iter2.set_op(UnifiedOp::Match(l2 - skip));
                Some(Ok(AlnCmpOp::Equal))
            }
            (Some(UnifiedOp::Del(l1)), Some(UnifiedOp::Del(l2))) => {
                let skip = *l1.min(l2);
                self.iter1.set_op(UnifiedOp::Del(l1 - skip));
                self.iter2.set_op(UnifiedOp::Del(l2 - skip));
                Some(Ok(AlnCmpOp::Equal))
            }
            (Some(UnifiedOp::Ins(l1)), Some(UnifiedOp::Ins(l2))) => {
                let skip = *l1.min(l2);
                self.read_i += skip as usize;
                self.iter1.set_op(UnifiedOp::Ins(l1 - skip));
                self.iter2.set_op(UnifiedOp::Ins(l2 - skip));
                Some(Ok(AlnCmpOp::Equal))
            }
            (Some(UnifiedOp::Mis(l1)), Some(UnifiedOp::Mis(l2))) => {
                let skip = *l1.min(l2);
                self.read_i += skip as usize;
                self.iter1.set_op(UnifiedOp::Mis(l1 - skip));
                self.iter2.set_op(UnifiedOp::Mis(l2 - skip));
                Some(Ok(AlnCmpOp::Equal))
            }
            (Some(UnifiedOp::Insertion(sv1)), Some(UnifiedOp::Insertion(sv2))) => {
                //TODO: lookup variants to apply adapted score. Only if neither match, dothis:

                let skip = sv1.len().min(sv2.len());
                self.read_i += skip as usize;
                self.iter1.set_op(UnifiedOp::Ins(l1 - skip));
                self.iter2.set_op(UnifiedOp::Ins(l2 - skip));
                Some(Ok(AlnCmpOp::Equal))
            }
            (Some(UnifiedOp::Deletion(l1)), Some(UnifiedOp::Deletion(l2))) => {

                Some(Ok(self.handle_equal_ops_del(l1, l2)))
            }

            (Some(UnifiedOp::HardClip(_) | UnifiedOp::Pad(_) | UnifiedOp::RefSkip(_)), Some(cigar2)) => {
                self.iter2.set_cigar_op(cigar2); // Put c2 back, consume c1
                self.next()
            }
            (Some(cigar1), Some(UnifiedOp::HardClip(_) | UnifiedOp::Pad(_) | UnifiedOp::RefSkip(_))) => {
                self.iter1.set_cigar_op(cigar1); // Put c1 back, consume c2
                self.next()
            }

            (Some(cigar1), None) => self.handle_differing_left(cigar1),
            (None, Some(cigar2)) => self.handle_differing_right(cigar2),
            (Some(cigar1), Some(cigar2)) => self.handle_differing_both(cigar1, cigar2),

            (None, None) => None,
        }
    }
}




/*pub struct AlignmentIterator<'a> {
    cigar_iter: Iter<'a, UnifiedOp>,
    md_iter: MdOpIterator<'a>,
    current_cigar_op: Option<UnifiedOp>,
    current_md_op: Option<MdOp>,
}

impl<'a> AlignmentIterator<'a> {
    #[must_use]
    pub fn new(cigar_iter: Iter<'a, UnifiedOp>, md_iter: MdOpIterator<'a>) -> Self {
        AlignmentIterator {
            cigar_iter,
            md_iter,
            current_cigar_op: None,
            current_md_op: None,
        }
    }
    pub fn next_md_op(&mut self) -> Option<MdOp> {
        self.current_md_op
            .take()
            .filter(|o| !o.is_empty())
            .or(self.md_iter.next())
    }
    pub fn next_cigar_op(&mut self) -> Option<UnifiedOp> {
        self.current_cigar_op
            .take()
            .filter(|o| !o.is_empty())
            .or(self.cigar_iter.next().copied())
    }
    pub fn set_cigar_op(&mut self, cigar: UnifiedOp) {
        self.current_cigar_op = Some(cigar);
    }
    pub fn set_md_op(&mut self, md_op: MdOp) {
        self.current_md_op = Some(md_op);
    }
    pub fn consume_md_op(&mut self) {
        if let Some(MdOp::Match(len)) = self.next_md_op() {
            self.set_md_op(MdOp::Match(len - 1));
        }
    }

    pub fn consume_md_ops(&mut self, mut remaining: u32) {
        while remaining > 0 {
            match self.next_md_op() {
                Some(MdOp::Match(len)) if len > remaining => {
                    self.set_md_op(MdOp::Match(len - remaining));
                    break;
                }
                Some(MdOp::Match(len)) => remaining -= len,
                Some(_) => remaining -= 1, // Mismatch or Deletion consumes 1
                None => break,
            }
        }
    }
    pub fn get_next_op(&mut self, cigar: UnifiedOp) -> Option<Result<UnifiedOp, AlignmentError>> {
        match cigar {
            UnifiedOp::Match(len) | UnifiedOp::Equal(len) | UnifiedOp::Diff(len) => {
                // Match, Equal and Diff are handled similarly.
                self.set_cigar_op(UnifiedOp::Match(len - 1));
                self.consume_md_op();

                match self.current_md_op {
                    Some(MdOp::Match(_)) => Some(Ok(UnifiedOp::Match)),
                    _ => Some(Ok(UnifiedOp::Mismatch)),
                }
            }
            // Insertions consume read bases but not reference bases.
            UnifiedOp::Insertion(len) => {
                self.set_cigar_op(UnifiedOp::Insertion(len - 1));
                Some(Ok(UnifiedOp::Insertion))
            }
            // Deletions consume reference bases but not read bases.
            UnifiedOp::Deletion(len) => {
                self.set_cigar_op(UnifiedOp::Deletion(len - 1));

                // Assume the deletion is valid even if MD tag is inconsistent.
                self.consume_md_op();
                Some(Ok(UnifiedOp::Deletion))
            }
            // Soft clips consume read bases but not reference bases.
            UnifiedOp::SoftClip(len) => {
                self.set_cigar_op(UnifiedOp::SoftClip(len - 1));
                Some(Ok(UnifiedOp::SoftClip))
            }
            UnifiedOp::RefSkip(len) => Some(Ok(UnifiedOp::RefSkip(len))),

            // These operations consume neither read nor reference bases. skip.
            UnifiedOp::HardClip(_) | UnifiedOp::Pad(_) => None,
        }
    }
}

impl Iterator for AlignmentIterator<'_> {
    type Item = Result<UnifiedOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(cigar) = self.next_cigar_op() {
            if let Some(res) = self.get_next_op(cigar) {
                return Some(res);
            }
        }
        None
    }
}*/


