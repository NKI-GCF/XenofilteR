use crate::alignment::AlignmentError;
use noodles::sam::alignment::record::cigar::op::{Op, Kind};
use crate::alignment::MdCigFlags;

pub(crate) enum BaseOp {
    Match,
    Mis,
    Del(usize),   // still worth grouping — no per-base work
    Ins(usize),
    RefSkip(usize),
    Relocate {
        pos: usize,
        penalty_score: f64,
    },
}

pub(crate) struct ScoreOpIter<'a> {
    md: &'a [u8],
    cigar: Box<dyn Iterator<Item = Result<Op, std::io::Error>> + 'a>,
    md_at: usize,
    md_match_remain: usize,
    /// Remaining length in the current CIGAR M-op we haven't emitted yet.
    cig_m_remain: usize,
}

impl<'a> ScoreOpIter<'a> {
    pub(crate) fn new(md_cig: &'a MdCigFlags) -> Self {
        let md = md_cig.get_md();
        let cigar = md_cig.get_cigar().iter();
        Self { md, cigar, md_at: 0, md_match_remain: 0, cig_m_remain: 0 }
    }
}

impl<'a> Iterator for ScoreOpIter<'a> {
    type Item = Result<BaseOp, AlignmentError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Drain remaining bases of the current M/X/= CIGAR op.
        if self.cig_m_remain == 0 {
            let op = self.cigar.next().and_then(|c| c.ok())?;
            match op.kind() {
                Kind::HardClip | Kind::Pad => self.next(),
                Kind::SoftClip => Some(Ok(BaseOp::Mis)), // 1 soft-clipped base
                Kind::Insertion => Some(Ok(BaseOp::Ins(op.len()))),
                Kind::Skip => Some(Ok(BaseOp::RefSkip(op.len()))),
                Kind::Deletion => Some(self.skip_md_deletion(op.len()).map(|()| BaseOp::Del(op.len()))),
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    self.cig_m_remain = op.len() - 1;
                    Some(self.next_md_base())
                }
            }
        } else {
            self.cig_m_remain -= 1;
            Some(self.next_md_base())
        }
    }
}

impl<'a> ScoreOpIter<'a> {
    fn next_md_base(&mut self) -> Result<BaseOp, AlignmentError> {
        if self.md_match_remain != 0 {
            let md = self.md.get(self.md_at);
            self.md_at += 1;
            match md {
                Some(b'A' | b'C' | b'G' | b'T' | b'N') => Ok(BaseOp::Mis),
                Some(n) if n.is_ascii_digit() => {
                    let mut num = (n - b'0') as usize;
                    while let Some(&d) = self.md.get(self.md_at) {
                        if !d.is_ascii_digit() { break; }
                        num = num * 10 + (d - b'0') as usize;
                        self.md_at += 1;
                    }
                    // num includes this base; consume 1 now, buffer the rest.
                    self.md_match_remain = num.saturating_sub(1);
                    Ok(BaseOp::Match)
                }
                other => Err(AlignmentError::MdCigMis(None, other.copied())),
            }
        } else {
            self.md_match_remain -= 1;
            Ok(BaseOp::Match)
        }
    }

    fn skip_md_deletion(&mut self, mut cig_remain: usize) -> Result<(), AlignmentError> {
        match self.md.get(self.md_at) {
            Some(b'^') => {
                self.md_at += 1;
                cig_remain += self.md_at;
                while let Some(b) = self.md.get(self.md_at) {
                    if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'N') { break; }
                    self.md_at += 1;
                }
                if cig_remain < self.md_at {
                    Ok(())
                } else {
                    Err(AlignmentError::MdCigMis(None, None))
                }
            }
            other => Err(AlignmentError::MdCigMis(None, other.copied())),
        }
    }
}

#[cfg(test)]
impl<'a> std::fmt::Debug for ScoreOpIter<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ScoreOpIter")
            .field("cig_m_remain", &self.cig_m_remain)
            .field("cigar", &self.cigar)
            .finish()
    }
}

#[cfg(test)]
pub(crate) mod tests;
