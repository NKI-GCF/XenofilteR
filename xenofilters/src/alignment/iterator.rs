use rust_htslib::bam::record::Cigar;
use anyhow::Result;
use super::{AlignmentOp, MdOp, MdOpIterator, AlignmentError};
use std::slice::Iter;

pub struct AlignmentIterator<'a> {
    cigar_iter: Iter<'a, Cigar>,
    md_iter: MdOpIterator<'a>,
    qual: &'a [u8],
    read_i: usize,
    current_cigar_op: Option<Cigar>,
    current_md_op: Option<MdOp>,
}

impl<'a> AlignmentIterator<'a> {
    pub fn new(cigar_iter: Iter<'a, Cigar>, md_iter: MdOpIterator<'a>, qual: &'a [u8]) -> Self {
        AlignmentIterator {
            cigar_iter,
            md_iter,
            qual,
            read_i: 0,
            current_cigar_op: None,
            current_md_op: None,
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
    fn next_md_op(&mut self) -> Option<MdOp> {
        if self.current_md_op.is_none() {
            self.current_md_op = self.md_iter.next();
        }
        self.current_md_op.take()
    }
}

impl<'a> Iterator for AlignmentIterator<'a> {
    type Item = Result<AlignmentOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        // This loop continues until an AlignmentOp is produced or the iterators are exhausted.
        // It allows us to consume CIGAR ops like HardClip or RefSkip that don't yield an op.
        loop {
            if self.current_cigar_op.is_none() {
                self.current_cigar_op = self.cigar_iter.next().map(|&cigar| cigar);
            }

            if let Some(cigar) = self.current_cigar_op.take() {
                match cigar {
                    Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                        // MAtch Equal and Diff are handled similarly.
                        self.current_cigar_op = (len > 1).then_some(Cigar::Match(len - 1));
                        match self.next_md_op().take() {
                            Some(MdOp::Match(n)) => {
                                if n > 1 {
                                    self.current_md_op = Some(MdOp::Match(n - 1));
                                }
                                return Some(self.get_qual().map(AlignmentOp::Match));
                            }
                            _ => return Some(self.get_qual().map(AlignmentOp::Mismatch)),
                        }
                    },
                    // Insertions consume read bases but not reference bases.
                    Cigar::Ins(len) => {
                        self.current_cigar_op = (len > 1).then_some(Cigar::Ins(len - 1));
                        return Some(self.get_qual().map(AlignmentOp::Insertion));
                    }
                    // Deletions consume reference bases but not read bases.
                    Cigar::Del(len) => {
                        self.current_cigar_op = (len > 1).then_some(Cigar::Del(len - 1));
                        if let Some(MdOp::Deletion) = self.next_md_op() {
                            // It matches, so consume it.
                            self.current_md_op = None;
                        }
                        // Assume the deletion is valid even if MD tag is inconsistent.
                        return Some(Ok(AlignmentOp::Deletion));
                    }
                    // Soft clips consume read bases but not reference bases.
                    Cigar::SoftClip(len) => {
                        self.current_cigar_op = (len > 1).then_some(Cigar::SoftClip(len - 1));
                        return Some(self.get_qual().map(AlignmentOp::SoftClip));
                    }

                    // These operations consume neither read nor reference bases. skip.
                    Cigar::HardClip(_) | Cigar::Pad(_) | Cigar::RefSkip(_) => {},
                }
            } else {
                return None;
            }
        }
    }
}
