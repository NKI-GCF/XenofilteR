use super::{AlignmentError, AlignmentOp, AlnCmpOp, MdOp, MdOpIterator};
use anyhow::Result;
use rust_htslib::bam::record::Cigar;
use std::slice::Iter;

pub struct AlignmentIterator<'a> {
    cigar_iter: Iter<'a, Cigar>,
    md_iter: MdOpIterator<'a>,
    current_cigar_op: Option<Cigar>,
    current_md_op: Option<MdOp>,
}

impl<'a> AlignmentIterator<'a> {
    #[must_use]
    pub fn new(cigar_iter: Iter<'a, Cigar>, md_iter: MdOpIterator<'a>) -> Self {
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
    pub fn next_cigar_op(&mut self) -> Option<Cigar> {
        self.current_cigar_op
            .take()
            .filter(|o| !o.is_empty())
            .or(self.cigar_iter.next().copied())
    }
    pub fn set_cigar_op(&mut self, cigar: Cigar) {
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
}

impl Iterator for AlignmentIterator<'_> {
    type Item = Result<AlignmentOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        // This loop continues until an AlignmentOp is produced or the iterators are exhausted.
        // It allows us to consume CIGAR ops like HardClip or RefSkip that don't yield an op.
        while let Some(cigar) = self.next_cigar_op() {
            match cigar {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    // Match, Equal and Diff are handled similarly.
                    self.set_cigar_op(Cigar::Match(len - 1));
                    self.consume_md_op();
                    match self.current_md_op {
                        Some(MdOp::Match(_)) => return Some(Ok(AlignmentOp::Match)),
                        _ => return Some(Ok(AlignmentOp::Mismatch)),
                    }
                }
                // Insertions consume read bases but not reference bases.
                Cigar::Ins(len) => {
                    self.set_cigar_op(Cigar::Ins(len - 1));
                    return Some(Ok(AlignmentOp::Insertion));
                }
                // Deletions consume reference bases but not read bases.
                Cigar::Del(len) => {
                    self.set_cigar_op(Cigar::Del(len - 1));

                    // Assume the deletion is valid even if MD tag is inconsistent.
                    self.consume_md_op();
                    return Some(Ok(AlignmentOp::Deletion));
                }
                // Soft clips consume read bases but not reference bases.
                Cigar::SoftClip(len) => {
                    self.set_cigar_op(Cigar::SoftClip(len - 1));
                    return Some(Ok(AlignmentOp::SoftClip));
                }

                // These operations consume neither read nor reference bases. skip.
                Cigar::HardClip(_) | Cigar::Pad(_) | Cigar::RefSkip(_) => {}
            }
        }
        None
    }
}

pub struct AlignmentCompareIterator<'a> {
    iter1: AlignmentIterator<'a>,
    iter2: AlignmentIterator<'a>,
    qual: &'a [u8],
    read_i: usize,
}

impl<'a> AlignmentCompareIterator<'a> {
    #[must_use]
    pub fn new(
        cigar_iter1: Iter<'a, Cigar>,
        md_iter1: MdOpIterator<'a>,
        cigar_iter2: Iter<'a, Cigar>,
        md_iter2: MdOpIterator<'a>,
        qual: &'a [u8],
    ) -> Self {
        AlignmentCompareIterator {
            iter1: AlignmentIterator::new(cigar_iter1, md_iter1),
            iter2: AlignmentIterator::new(cigar_iter2, md_iter2),
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
    fn next_md_ops(&mut self) -> (Option<MdOp>, Option<MdOp>) {
        (self.iter1.next_md_op(), self.iter2.next_md_op())
    }
    fn next_cigar_ops(&mut self) -> (Option<Cigar>, Option<Cigar>) {
        (self.iter1.next_cigar_op(), self.iter2.next_cigar_op())
    }

    /// Handles the (Match/Eq/Diff, Match/Eq/Diff) fast-path case.
    fn handle_match_match(&mut self, len1: u32, len2: u32) -> Result<AlnCmpOp, AlignmentError> {
        let len = len1.min(len2);
        let mut skip = 1;

        let res = match self.next_md_ops() {
            // Both are matches, fast-skip
            (Some(MdOp::Match(n1)), Some(MdOp::Match(n2))) => {
                let n = n1.min(n2);
                skip = len.min(n);
                self.read_i += skip as usize;
                self.iter1.set_md_op(MdOp::Match(n1 - skip));
                self.iter2.set_md_op(MdOp::Match(n2 - skip));
                Ok(AlnCmpOp::Equal)
            }
            // 1 is match, 2 is mismatch
            (Some(MdOp::Match(n1)), _) => {
                self.iter1.set_md_op(MdOp::Match(n1 - 1));
                self.get_qual()
                    .map(|q| AlnCmpOp::UnEqual(AlignmentOp::Match, AlignmentOp::Mismatch, q))
            }
            // 1 is mismatch, 2 is match
            (_, Some(MdOp::Match(n2))) => {
                self.iter2.set_md_op(MdOp::Match(n2 - 1));
                self.get_qual()
                    .map(|q| AlnCmpOp::UnEqual(AlignmentOp::Mismatch, AlignmentOp::Match, q))
            }
            // Both are mismatches
            (_, _) => {
                self.read_i += 1;
                Ok(AlnCmpOp::Equal)
            }
        };

        self.iter1.set_cigar_op(Cigar::Match(len1 - skip));
        self.iter2.set_cigar_op(Cigar::Match(len2 - skip));
        res
    }

    /// Handles identical operations that consume read (Ins, SoftClip).
    fn handle_equal_ops_read(&mut self, len1: u32, len2: u32, op_fn: fn(u32) -> Cigar) -> AlnCmpOp {
        let len = len1.min(len2);
        self.read_i += len as usize;
        self.iter1.set_cigar_op(op_fn(len1 - len));
        self.iter2.set_cigar_op(op_fn(len2 - len));
        AlnCmpOp::Equal
    }

    /// Handles identical (Del, Del) operations.
    fn handle_equal_ops_del(&mut self, len1: u32, len2: u32) -> AlnCmpOp {
        let len = len1.min(len2);
        self.iter1.set_cigar_op(Cigar::Del(len1 - len));
        self.iter2.set_cigar_op(Cigar::Del(len2 - len));
        self.iter1.consume_md_ops(len);
        self.iter2.consume_md_ops(len);
        AlnCmpOp::Equal
    }

    /// Handles the "slow path" where only iter1 has an operation.
    fn handle_differing_left(&mut self, cigar1: Cigar) -> Option<Result<AlnCmpOp, AlignmentError>> {
        self.iter1.set_cigar_op(cigar1); // undo to redo:
        match self.iter1.next() {
            Some(Ok(op)) => Some(self.get_qual().map(|q| AlnCmpOp::Left(op, q))),
            Some(Err(e)) => Some(Err(e)),
            None => None, // Signals the loop to `break`
        }
    }

    /// Handles the "slow path" where only iter2 has an operation.
    fn handle_differing_right(
        &mut self,
        cigar2: Cigar,
    ) -> Option<Result<AlnCmpOp, AlignmentError>> {
        self.iter2.set_cigar_op(cigar2); // undo to redo:
        match self.iter2.next() {
            Some(Ok(op)) => Some(self.get_qual().map(|q| AlnCmpOp::Right(op, q))),
            Some(Err(e)) => Some(Err(e)),
            None => None, // Signals the loop to `break`
        }
    }

    /// Handles the "slow path" where both iters have differing operations.
    fn handle_differing_both(
        &mut self,
        cigar1: Cigar,
        cigar2: Cigar,
    ) -> Option<Result<AlnCmpOp, AlignmentError>> {
        self.iter1.set_cigar_op(cigar1);
        self.iter2.set_cigar_op(cigar2); // undo both to redo:
        match self.iter1.next() {
            Some(Ok(op1)) => match self.iter2.next() {
                Some(Ok(op2)) => Some(self.get_qual().map(|q| AlnCmpOp::UnEqual(op1, op2, q))),
                Some(Err(e)) => Some(Err(e)),
                None => None, // iter2 ran out, break
            },
            Some(Err(e)) => Some(Err(e)),
            None => None, // iter1 ran out, break
        }
    }
}

impl<'a> Iterator for AlignmentCompareIterator<'a> {
    type Item = Result<AlnCmpOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.next_cigar_ops() {
                (
                    Some(Cigar::Match(l1) | Cigar::Equal(l1) | Cigar::Diff(l1)),
                    Some(Cigar::Match(l2) | Cigar::Equal(l2) | Cigar::Diff(l2)),
                ) => {
                    return Some(self.handle_match_match(l1, l2));
                }
                (Some(Cigar::Ins(l1)), Some(Cigar::Ins(l2))) => {
                    return Some(Ok(self.handle_equal_ops_read(l1, l2, Cigar::Ins)));
                }
                (Some(Cigar::Del(l1)), Some(Cigar::Del(l2))) => {
                    return Some(Ok(self.handle_equal_ops_del(l1, l2)));
                }
                (Some(Cigar::SoftClip(l1)), Some(Cigar::SoftClip(l2))) => {
                    return Some(Ok(self.handle_equal_ops_read(l1, l2, Cigar::SoftClip)));
                }

                (Some(Cigar::HardClip(_) | Cigar::Pad(_) | Cigar::RefSkip(_)), Some(cigar2)) => {
                    self.iter2.set_cigar_op(cigar2); // Put c2 back, consume c1
                }
                (Some(cigar1), Some(Cigar::HardClip(_) | Cigar::Pad(_) | Cigar::RefSkip(_))) => {
                    self.iter1.set_cigar_op(cigar1); // Put c1 back, consume c2
                }

                (Some(cigar1), None) => {
                    if let Some(item) = self.handle_differing_left(cigar1) {
                        return Some(item);
                    }
                    break; // iter1 ran out
                }
                (None, Some(cigar2)) => {
                    if let Some(item) = self.handle_differing_right(cigar2) {
                        return Some(item);
                    }
                    break; // iter2 ran out
                }
                (Some(cigar1), Some(cigar2)) => {
                    if let Some(item) = self.handle_differing_both(cigar1, cigar2) {
                        return Some(item);
                    }
                    break; // One of the iters ran out
                }

                (None, None) => break,
            }
        }
        None
    }
}
