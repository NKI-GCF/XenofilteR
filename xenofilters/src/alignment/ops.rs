use crate::alignment::AlignmentError;
use noodles::sam::alignment::record::cigar::{op::Kind, Op};
use smallvec::SmallVec;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum MdOp {
    Match(usize),
    Mismatch(u8),
    Deletion(SmallVec<[u8; 4]>), // slightly larger stack buffer
    Error(u8),
}

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

pub(crate) struct UnifiedOpIterator<'a> {
    cigar: &'a [Op],
    cig_at: usize,
    md: &'a [u8],
    md_at: usize,
    match_remain: isize,
}

impl<'a> UnifiedOpIterator<'a> {
    pub(crate) fn new(cigar: &'a [Op], md: &'a [u8]) -> Self {
        Self {
            cigar,
            cig_at: 0,
            md,
            md_at: 0,
            match_remain: 0,
        }
    }

    fn cig_match(&mut self, cig_len: usize) -> Result<UnifiedOp, AlignmentError> {
        match self.md_next() {
            Some(MdOp::Match(md_len)) => {
                self.match_remain = (md_len - cig_len) as isize;
                Ok(UnifiedOp::Match(if self.match_remain > 0 { cig_len } else { md_len }))
            }
            Some(MdOp::Mismatch(md_len)) => {
                if md_len != 1 || cig_len != 1 {
                    self.match_remain = cig_len as isize - 1;
                }
                Ok(UnifiedOp::Mis(1))
            }
            md_op => Err(AlignmentError::MdCigMis(Some(Op::new(Kind::Match, cig_len)), md_op)),
        }
    }

    fn md_next(&mut self) -> Option<MdOp> {
        let b = self.md.get(self.md_at)?;
        self.md_at += 1;

        match *b {
            b'A' | b'C' | b'G' | b'T' | b'N' => Some(MdOp::Mismatch(*b)),

            b'^' => {
                let mut deleted = SmallVec::new();
                while let Some(&b) = self.md.get(self.md_at) {
                    if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'N') {
                        break;
                    }
                    deleted.push(b);
                    self.md_at += 1;
                }
                Some(MdOp::Deletion(deleted))
            }

            n if n.is_ascii_digit() => {
                let mut num = (n - b'0') as usize;
                while let Some(&b) = self.md.get(self.md_at) && b.is_ascii_digit() {
                     num = num * 10 + (b - b'0') as usize;
                     self.md_at += 1;
                }
                Some(MdOp::Match(num))
            }
            x => Some(MdOp::Error(x)),
        }
    }
}

impl<'a> Iterator for UnifiedOpIterator<'a> {
    type Item = Result<UnifiedOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.match_remain < 0 {
            return Some(self.cig_match(-self.match_remain as usize));
        }
        while let Some(c) = self.cigar.get(self.cig_at).map(|c| (c.kind(), c)) {
            self.cig_at += 1;
            let res = match c {
                (Kind::HardClip | Kind::Pad, _) => continue,
                (Kind::Match| Kind::SequenceMatch | Kind::SequenceMismatch, c) => {
                    let cig_len = c.len();
                    if self.match_remain <= 0 {
                        self.cig_match(cig_len)
                    } else {
                        let md_len = self.match_remain as usize;
                        self.match_remain -= cig_len as isize;
                        Ok(UnifiedOp::Match(if self.match_remain > 0 { cig_len } else { md_len }))
                    }
                },
                (Kind::SoftClip, c) => Ok(UnifiedOp::Mis(c.len())),
                (Kind::Deletion, c) if self.match_remain <= 0 => {
                    match self.md_next() {
                        Some(MdOp::Deletion(del_seq)) if del_seq.len() == c.len() => {
                            Ok(UnifiedOp::Del(c.len()))
                        }
                        other => Err(AlignmentError::MdCigMis(Some(*c), other)),
                    }
                },
                (Kind::Deletion, c) => {
                    Err(AlignmentError::MdCigMis(Some(*c), Some(MdOp::Match(self.match_remain as usize))))
                },
                (kind, c) => UnifiedOp::try_from((kind, c.len())),
            };
            return Some(res);
        }
        if self.match_remain > 0 {
            Some(Err(AlignmentError::MdCigMis(None, Some(MdOp::Match(self.match_remain as usize)))))
        } else {
            let md_op = self.md_next()?;
            Some(Err(AlignmentError::MdCigMis(None, Some(md_op))))
        }
    }
}

#[cfg(test)]
impl<'a> std::fmt::Debug for UnifiedOpIterator<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("UnifiedOpIterator")
            .field("cig_match_remain", &self.cig_match_remain)
            .field("cigar", &self.cigar)
            .finish()
    }
}

#[cfg(test)]
pub(crate) mod tests;
