use crate::alignment::PrepareError;
use crate::alignment::{AlignmentError, MdOp, MdOpIterator};
use rust_htslib::bam::record::{Cigar, Record};
use std::cmp::min;

// TODO: for paired-end introduce another operation for read switching

#[derive(Debug, Clone, PartialEq)]
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
    next_md_op: Option<MdOp>,
    next_cig: Option<Cigar>,
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
            next_md_op: None,
            next_cig: None,
            is_rev: rec.is_reverse(),
        })
    }
    pub fn empty(is_rev: bool) -> Self {
        UnifiedOpIterator {
            cigar_iter: vec![].into_iter(),
            md_iter: MdOpIterator::empty(),
            next_op: None,
            next_md_op: None,
            next_cig: None,
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
}

impl<'a> Iterator for UnifiedOpIterator<'a> {
    type Item = Result<UnifiedOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_op.is_some() {
            return Some(Ok(self.next_op.take().unwrap()));
        }
        let next_cig = if self.next_cig.is_some() {
            self.next_cig.take()
        } else if self.is_rev {
            self.cigar_iter.next_back()
        } else {
            self.cigar_iter.next()
        };
        match next_cig? {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let next_md_op = if self.next_md_op.is_some() {
                    self.next_md_op.take().map(Ok)
                } else if self.is_rev {
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

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Aux, CigarString};

    pub fn create_record(
        qname: &[u8],
        vec_cig: Vec<Cigar>,
        qual: &[u8],
        md: &str,
        is_rev: bool,
    ) -> rust_htslib::bam::Record {
        let mut record = Record::new();
        let cigar_string = CigarString(vec_cig);

        // Set basic fields
        record.set(qname, Some(&cigar_string), &vec![b'A'; qual.len()], qual);
        if is_rev {
            record.set_reverse();
        }

        // HTSlib requires the MD tag to be set manually for tests
        record.push_aux(b"MD", Aux::String(md)).unwrap();
        record
    }

    #[test]
    fn simple_match_cig_10m_md_10() {
        let rec = create_record(b"read1", vec![Cigar::Match(10)], &vec![30; 10], "10", false);
        let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
        let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
        assert_eq!(ops, vec![UnifiedOp::Match(10)]);

        let rec = create_record(b"read2", vec![Cigar::Match(10)], &vec![30; 10], "10", true);
        let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
        let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
        assert_eq!(ops, vec![UnifiedOp::Match(10)]);
    }

    #[test]
    fn mismatch_cig_10m_md_5a4() {
        let rec = create_record(
            b"read1",
            vec![Cigar::Match(10)],
            &vec![30; 10],
            "5A4",
            false,
        );
        let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
        let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
        assert_eq!(
            ops,
            vec![UnifiedOp::Match(5), UnifiedOp::Mis(1), UnifiedOp::Match(4)]
        );

        let rec = create_record(b"read2", vec![Cigar::Match(10)], &vec![30; 10], "5A4", true);
        let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
        let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
        assert_eq!(
            ops,
            vec![UnifiedOp::Match(4), UnifiedOp::Mis(1), UnifiedOp::Match(5)]
        );
    }

    #[test]
    fn indels_cig_5m2i3m_md_10() -> Result<(), PrepareError> {
        let rec = create_record(
            b"read3",
            vec![Cigar::Match(5), Cigar::Ins(2), Cigar::Match(3)],
            &vec![30; 10],
            "10",
            false,
        );
        let uop_iter = UnifiedOpIterator::new(&rec)?;
        let ops: Vec<UnifiedOp> = uop_iter
            .map(|r| r)
            .collect::<Result<Vec<UnifiedOp>, AlignmentError>>()?;
        assert_eq!(
            ops,
            vec![UnifiedOp::Match(5), UnifiedOp::Ins(2), UnifiedOp::Match(3)]
        );
        Ok(())
    }

    #[test]
    fn soft_clip_cig_5s5m_md_5() {
        let rec = create_record(
            b"read4",
            vec![Cigar::SoftClip(5), Cigar::Match(5)],
            &vec![30; 10],
            "5",
            false,
        );
        let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
        let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
        assert_eq!(ops, vec![UnifiedOp::Mis(5), UnifiedOp::Match(5)]);
    }

    #[test]
    fn deletion_cig_5m3d5m_md_5daaa5() -> Result<(), PrepareError> {
        let rec = create_record(
            b"read5",
            vec![Cigar::Match(5), Cigar::Del(3), Cigar::Match(5)],
            &vec![30; 10],
            "5^AAA5",
            false,
        );
        let uop_iter = UnifiedOpIterator::new(&rec)?;
        let ops: Vec<UnifiedOp> = uop_iter
            .map(|r| r)
            .collect::<Result<Vec<UnifiedOp>, AlignmentError>>()?;
        assert_eq!(
            ops,
            vec![UnifiedOp::Match(5), UnifiedOp::Del(3), UnifiedOp::Match(5)]
        );
        Ok(())
    }

    #[test]
    fn test_length_invariants() {
        // "50M10I40M"
        let rec = create_record(
            b"read1",
            vec![Cigar::Match(50), Cigar::Ins(10), Cigar::Match(40)],
            &[30; 100],
            "90",
            false,
        );
        let iter = UnifiedOpIterator::new(&rec).unwrap();

        let total_len: u32 = iter.map(|r| r.unwrap().len()).sum();
        // Sum of Match/Mis/Ins/Del should match the logical alignment footprint
        assert_eq!(total_len, 100);
    }
}
