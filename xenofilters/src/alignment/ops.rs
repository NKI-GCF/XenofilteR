use crate::alignment::MdCigFlags;
use crate::Error;
use noodles::sam::alignment::record::cigar::op::{Kind, Op};

pub(crate) enum BaseOp {
    Match,
    Mis,
    Del(usize), // still worth grouping — no per-base work
    Ins(usize),
    Clip(usize),
    RefSkip(usize),
    /// C→T on the forward strand or G→A on the reverse strand.
    /// Caused by bisulfite conversion of unmethylated cytosines; not a
    /// true substitution error.  Scored with zero penalty (Penalty::log_likelihood_bisulfite).
    BisulfiteConversion,
}

pub(crate) struct ScoreOpIter<'a> {
    md: &'a [u8],
    cigar: Box<dyn Iterator<Item = Result<Op, std::io::Error>> + 'a>,
    md_at: usize,
    md_match_remain: usize,
    /// Remaining length in the current CIGAR M-op we haven't emitted yet.
    cig_m_remain: usize,
    // Bisulfite support
    seq: Option<&'a [u8]>,  // read sequence; None → bisulfite detection disabled
    is_reverse: bool,              // true = reverse-complemented read
    read_pos: usize,               // current position in read sequence
}

impl<'a> ScoreOpIter<'a> {
    pub(crate) fn new(md_cig: &'a MdCigFlags) -> Self {
        let md = md_cig.get_md();
        let cigar = md_cig.get_cigar().iter();
        Self {
            md,
            cigar,
            md_at: 0,
            md_match_remain: 0,
            cig_m_remain: 0,
            seq: md_cig.seq(),
            is_reverse: md_cig.is_reverse_complemented(),
            read_pos: 0,
        }
    }

    /// Classify a mismatch: bisulfite conversion or genuine substitution.
    /// `ref_base`: the reference base from the MD tag (ASCII uppercase).
    fn classify_mismatch(&mut self, ref_base: u8) -> BaseOp {
        let Some(seq) = self.seq else { return BaseOp::Mis; };
        let read_base = seq.get(self.read_pos).copied().map(|b| b.to_ascii_uppercase());
        let Some(rb) = read_base else { return BaseOp::Mis; };

        let is_conversion = if !self.is_reverse {
            // Forward strand: bisulfite converts C → T.
            ref_base == b'C' && rb == b'T'
        } else {
            // Reverse strand: bisulfite converts C → T on the original strand,
            // which appears as G → A in the read-as-sequenced orientation.
            ref_base == b'G' && rb == b'A'
        };

        if is_conversion {
            BaseOp::BisulfiteConversion
        } else {
            BaseOp::Mis
        }
    }
}

impl<'a> Iterator for ScoreOpIter<'a> {
    type Item = Result<BaseOp, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        // Drain remaining bases of the current M/X/= CIGAR op.
        if self.cig_m_remain == 0 {
            let op = self.cigar.next().and_then(|c| c.ok())?;
            match op.kind() {
                Kind::HardClip | Kind::Pad => self.next(),
                Kind::SoftClip => {
                    self.read_pos += op.len();
                    Some(Ok(BaseOp::Clip(op.len())))
                }
                Kind::Insertion => {
                    self.read_pos += op.len();
                    Some(Ok(BaseOp::Ins(op.len())))
                }
                Kind::Skip => Some(Ok(BaseOp::RefSkip(op.len()))),
                Kind::Deletion => Some(
                    self.skip_md_deletion(op.len())
                        .map(|()| BaseOp::Del(op.len())),
                ),
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
    fn next_md_base(&mut self) -> Result<BaseOp, Error> {
        if self.md_match_remain == 0 {
            let md = self.md.get(self.md_at);
            self.md_at += 1;
            match md {
                Some(b'A' | b'C' | b'G' | b'T' | b'N') => {
                    let ref_base = md.copied().unwrap();
                    self.read_pos += 1;
                    Ok(self.classify_mismatch(ref_base))
                }
                Some(n) if n.is_ascii_digit() => {
                    let mut num = (n - b'0') as usize;
                    while let Some(&d) = self.md.get(self.md_at) {
                        if !d.is_ascii_digit() {
                            break;
                        }
                        num = num * 10 + (d - b'0') as usize;
                        self.md_at += 1;
                    }
                    // num includes this base; consume 1 now, buffer the rest.
                    self.md_match_remain = num.saturating_sub(1);
                    self.read_pos += 1;
                    Ok(BaseOp::Match)
                }
                other => Err(Error::MdCigMis(None, other.copied())),
            }
        } else {
            self.md_match_remain -= 1;
            self.read_pos += 1;
            Ok(BaseOp::Match)
        }
    }

    fn skip_md_deletion(&mut self, cig_remain: usize) -> Result<(), Error> {
        match self.md.get(self.md_at) {
            Some(b'^') => {
                self.md_at += 1;
                let md_start = self.md_at;
                while let Some(b) = self.md.get(self.md_at) {
                    if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'N') {
                        break;
                    }
                    self.md_at += 1;
                }
                if self.md_at - md_start == cig_remain {
                    Ok(())
                } else {
                    Err(Error::MdCigMis(None, None))
                }
            }
            other => Err(Error::MdCigMis(None, other.copied())),
        }
    }
}

#[cfg(test)]
pub(crate) mod tests;
