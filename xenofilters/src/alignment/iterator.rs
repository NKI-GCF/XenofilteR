use rust_htslib::bam::record::{Cigar, CigarStringView};
use anyhow::Result;
use super::{AlignmentOp, MdOp, MdOpIterator, AlignmentError};

pub struct AlignmentIterator<'a> {
    cigar: CigarStringView,
    md_iter: MdOpIterator<'a>,
    qual: &'a [u8],
    read_i: usize,
    current_cigar_op: Option<(Cigar, u32)>,
    current_md_op: Option<MdOp>,
}

impl<'a> AlignmentIterator<'a> {
    pub fn new(cigar: CigarStringView, md_iter: MdOpIterator<'a>, qual: &'a [u8]) -> Self {
        AlignmentIterator {
            cigar,
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
                match self.cigar.iter().next() {
                    Some(&cigar_view) => {
                        let cigar_len = cigar_view.len();
                        self.current_cigar_op = Some((cigar_view, cigar_len));
                    }
                    None => return None,
                }
            }

            if let Some((ref cigar, ref mut len)) = self.current_cigar_op {
                let op_result = match cigar {
                    Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                        if self.current_md_op.is_none() {
                            self.current_md_op = self.md_iter.next();
                        }

                        let q = self.qual[self.read_i];
                        self.read_i += 1;

                        match self.current_md_op.take() {
                            Some(MdOp::Match(n)) => {
                                if n > 1 {
                                    self.current_md_op = Some(MdOp::Match(n - 1));
                                }
                                Ok(AlignmentOp::Match(q))
                            }
                            Some(MdOp::Mismatch) => Ok(AlignmentOp::Mismatch(q)),
                            // An MD Deletion cannot occur during a CIGAR Match.
                            Some(MdOp::Deletion) => Err(AlignmentError::MismatchedDeletion),
                            None => Err(AlignmentError::UnexpectedMdEnd),
                        }
                    }

                    Cigar::Ins(_) => {
                        let q = self.qual[self.read_i];
                        self.read_i += 1;
                        Ok(AlignmentOp::Insertion(q))
                    }

                    Cigar::Del(_) => {
                        if self.current_md_op.is_none() {
                            self.current_md_op = self.md_iter.next();
                        }
                        match self.current_md_op.take() {
                            Some(MdOp::Deletion) => Ok(AlignmentOp::Deletion),
                            _ => Err(AlignmentError::MissingMdDeletion),
                        }
                    }

                    Cigar::SoftClip(_) => {
                        let q = self.qual[self.read_i];
                        self.read_i += 1;
                        Ok(AlignmentOp::SoftClip(q))
                    }

                    // These operations consume neither read nor reference bases. skip.
                    Cigar::HardClip(_) | Cigar::Pad(_) | Cigar::RefSkip(_) => {
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
