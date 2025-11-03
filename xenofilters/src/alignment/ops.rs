use crate::{LOG_LIKELIHOOD_MATCH, MAX_Q};

#[derive(Debug, Clone)]
pub enum AlignmentOp {
    Match,
    Mismatch,
    Insertion,
    Deletion,
    SoftClip,
    RefSkip(u32),
}

impl AlignmentOp {
    #[must_use]
    pub fn same_variant(&self, other: &AlignmentOp) -> bool {
        matches!(
            (self, other),
            (AlignmentOp::Match, AlignmentOp::Match)
                | (AlignmentOp::Mismatch, AlignmentOp::Mismatch)
                | (AlignmentOp::Insertion, AlignmentOp::Insertion)
                | (AlignmentOp::Deletion, AlignmentOp::Deletion)
                | (AlignmentOp::SoftClip, AlignmentOp::SoftClip)
        )
    }
    pub fn score(
        self,
        q: u8,
        indel_gap: &mut Option<bool>,
        log_likelihood_mismatch: &[f64; MAX_Q + 2],
    ) -> f64 {
        let gap_open = log_likelihood_mismatch[MAX_Q];
        let gap_ext = log_likelihood_mismatch[MAX_Q + 1];
        match self {
            AlignmentOp::Match => {
                *indel_gap = None;
                LOG_LIKELIHOOD_MATCH[q as usize]
            }
            AlignmentOp::Mismatch | AlignmentOp::SoftClip => {
                *indel_gap = None;
                log_likelihood_mismatch[q as usize]
            }
            AlignmentOp::Insertion => {
                if *indel_gap == Some(true) {
                    gap_ext
                } else {
                    *indel_gap = Some(true);
                    gap_open + gap_ext
                }
            }
            AlignmentOp::Deletion => {
                if *indel_gap == Some(false) {
                    gap_ext
                } else {
                    *indel_gap = Some(false);
                    gap_open + gap_ext
                }
            }
            AlignmentOp::RefSkip(_) => 0.0,
        }
    }
    pub fn ref_consumed(&self) -> usize {
        match self {
            AlignmentOp::Match | AlignmentOp::Mismatch | AlignmentOp::Deletion => 1,
            AlignmentOp::RefSkip(n) => *n as usize,
            AlignmentOp::Insertion | AlignmentOp::SoftClip => 0,
        }
    }
    pub fn consumes_read(&self) -> bool {
        matches!(
            self,
            AlignmentOp::Match
                | AlignmentOp::Mismatch
                | AlignmentOp::Insertion
                | AlignmentOp::SoftClip
        )
    }
}

pub enum AlnCmpOp {
    Equal,
    UnEqual(AlignmentOp, AlignmentOp, u8),
    Left(AlignmentOp, u8),
    Right(AlignmentOp, u8),
}

#[derive(Debug, Clone)]
pub enum MdOp {
    Match(u32),
    Mismatch,
    Deletion,
}

impl MdOp {
    #[must_use]
    pub fn is_empty(&self) -> bool {
        match self {
            MdOp::Match(n) => *n == 0,
            _ => false,
        }
    }
}

#[derive(Clone)]
pub struct MdOpIterator<'a> {
    chars: std::str::Chars<'a>,
    // State is simpler: just the number being built and a flag.
    num: u32,
    in_deletion: bool,
}

impl<'a> MdOpIterator<'a> {
    #[must_use]
    pub fn new(md: &'a str) -> Self {
        MdOpIterator {
            chars: md.chars(),
            num: 0,
            in_deletion: false,
        }
    }
}

impl<'a> Iterator for MdOpIterator<'a> {
    type Item = MdOp;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let Some(c) = self.chars.next() else {
                if self.num > 0 {
                    let match_op = MdOp::Match(self.num);
                    self.num = 0;
                    return Some(match_op);
                }
                return None;
            };

            if c.is_ascii_digit() {
                self.num = (self.num * 10) + (c as u32 - '0' as u32);
                self.in_deletion = false;
                continue;
            }

            if self.num > 0 {
                let match_op = MdOp::Match(self.num);
                self.num = 0;
                if c == '^' {
                    self.in_deletion = true;
                }
                return Some(match_op);
            }

            if c == '^' {
                self.in_deletion = true;
                continue;
            }

            return self
                .in_deletion
                .then_some(MdOp::Deletion)
                .or(Some(MdOp::Mismatch));
        }
    }
}
