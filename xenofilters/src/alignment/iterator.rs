use rust_htslib::bam::record::Cigar;
use anyhow::Result;
use super::{AlignmentOp, MdOp, MdOpIterator, AlignmentError};
use std::slice::Iter;

pub struct AlignmentIterator<'a> {
    cigar_iter: Iter<'a, Cigar>,
    md_iter: MdOpIterator<'a>,
    qual: &'a [u8],
    read_i: usize,
    current_cigar_op: Option<(Cig, u32)>,
    current_md_op: Option<MdOp>,
}

enum Cig {
    Match,
    Ins,
    Del,
    SoftClip,
    Equal,
    Diff,
    HardClip,
    Pad,
    RefSkip,
}

impl From<Cigar> for Cig {
    fn from(cigar: Cigar) -> Self {
        match cigar {
            Cigar::Match(_) => Cig::Match,
            Cigar::Ins(_) => Cig::Ins,
            Cigar::Del(_) => Cig::Del,
            Cigar::SoftClip(_) => Cig::SoftClip,
            Cigar::Equal(_) => Cig::Equal,
            Cigar::Diff(_) => Cig::Diff,
            Cigar::HardClip(_) => Cig::HardClip,
            Cigar::Pad(_) => Cig::Pad,
            Cigar::RefSkip(_) => Cig::RefSkip,
        }
    }
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
}

impl<'a> Iterator for AlignmentIterator<'a> {
    type Item = Result<AlignmentOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        // This loop continues until an AlignmentOp is produced or the iterators are exhausted.
        // It allows us to consume CIGAR ops like HardClip or RefSkip that don't yield an op.
        loop {
            if self.current_cigar_op.is_none() {
                match self.cigar_iter.next() {
                    Some(&cigar_view) => {
                        let cigar_len = cigar_view.len();
                        self.current_cigar_op = Some((Cig::from(cigar_view), cigar_len));
                    }
                    None => return None,
                }
            }

            if let Some((ref cigar, ref mut len)) = self.current_cigar_op {
                let op_result = match cigar {
                    Cig::Match | Cig::Equal | Cig::Diff => {
                        if self.current_md_op.is_none() {
                            self.current_md_op = self.md_iter.next();
                        }

                        let q = *self.qual.get(self.read_i).unwrap();
                        self.read_i += 1;

                        match self.current_md_op.take() {
                            Some(MdOp::Match(n)) => {
                                if n > 1 {
                                    self.current_md_op = Some(MdOp::Match(n - 1));
                                }
                                Ok(AlignmentOp::Match(q))
                            }
                            _ => Ok(AlignmentOp::Mismatch(q)),
                        }
                    }
                    // Insertions consume read bases but not reference bases.
                    Cig::Ins => {
                        let q = *self.qual.get(self.read_i).unwrap();
                        self.read_i += 1;
                        Ok(AlignmentOp::Insertion(q))
                    }
                    // Deletions consume reference bases but not read bases.
                    Cig::Del => {
                        if self.current_md_op.is_none() {
                            self.current_md_op = self.md_iter.next();
                        }
                        if let Some(MdOp::Deletion) = self.current_md_op {
                            // It matches, so consume it.
                            self.current_md_op = None;
                        }
                        // Assume the deletion is valid even if MD tag is inconsistent.
                        Ok(AlignmentOp::Deletion)
                    }
                    // Soft clips consume read bases but not reference bases.
                    Cig::SoftClip => {
                        let q = *self.qual.get(self.read_i).unwrap();
                        self.read_i += 1;
                        Ok(AlignmentOp::SoftClip(q))
                    }

                    // These operations consume neither read nor reference bases. skip.
                    Cig::HardClip | Cig::Pad | Cig::RefSkip => {
                        *len = 0;
                        self.current_cigar_op = None;
                        continue;
                    },
                };

                *len -= 1;
                if *len == 0 {
                    self.current_cigar_op = None;
                }

                // Wrap the successful or failed result in Some and return.
                return Some(op_result);
            }
        }
    }
}
