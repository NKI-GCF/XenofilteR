use crate::alignment::AlignmentError;
use crate::{CONFIG, LOG_LIKELIHOOD_MATCH, LOG_LIKELIHOOD_MISMATCH};
use rust_htslib::bam::record::{Cigar, CigarStringView};
use smallvec::SmallVec;

use std::iter::{Peekable, once};
use std::str::Chars;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RichMdOp {
    Match(u32),
    Mismatch(u8),      // A single base mismatch
    Deletion(SmallVec<[u8; 1]>), // One or more deleted bases
}

#[derive(Clone)]
pub struct RichMdOpIterator<'a> {
    chars: Peekable<Chars<'a>>,
    match_remain: u32,
}

impl<'a> RichMdOpIterator<'a> {
    pub fn new(md: &'a str) -> Self {
        Self {
            chars: md.chars().peekable(),
            match_remain: 0,
        }
    }
}

impl<'a> Iterator for RichMdOpIterator<'a> {
    type Item = RichMdOp;

    fn next(&mut self) -> Option<Self::Item> {
        if self.match_remain > 0 {
            return Some(RichMdOp::Match({
                let len = self.match_remain;
                self.match_remain = 0;
                len
            }));
        }
        match self.chars.peek()? {
            n if n.is_ascii_digit() => {
                // Handle Match
                let mut num = *n as u32 - '0' as u32;
                self.chars.next(); // Consume first digit
                while let Some(&digit) = self.chars.peek().filter(|c| c.is_ascii_digit()) {
                    self.chars.next(); // Consume digit
                    num = (num * 10) + (digit as u32 - '0' as u32);
                }
                Some(RichMdOp::Match(num))
            }
            '^' => {
                // Handle Deletion
                self.chars.next(); // Consume '^'

                let mut deleted_bases: SmallVec<[u8; 1]> = SmallVec::new();
                while let Some(&base) = self.chars.peek().filter(|c| c.is_ascii_alphabetic()) {
                    self.chars.next(); // Consume base
                    deleted_bases.push(base as u8);
                }
                Some(RichMdOp::Deletion(deleted_bases))
            }
            _ => {
                // Handle Mismatch
                let c = self.chars.next().unwrap(); // Consume base
                Some(RichMdOp::Mismatch(c as u8))
            }
        }
    }
}

// TODO: for paired-end introduce another operation for read switching

#[derive(Debug, Clone)]
pub enum UnifiedOp {
    /// A pure reference match
    Match(usize),

    /// A refskip from CIGAR
    RefSkip(usize),

    // simple variants, not variant-aware
    Mis(usize),
    Ins(usize),
    Del(usize),
    Trans(f64),
    Mate,

    // variant-aware operations
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
    },
    MateSwitch {
        new_tid: i32,
        new_pos: i64,
        new_strand: bool,
    },
}

impl UnifiedOp {
    #[must_use]
    pub fn ref_consumed(&self) -> usize {
        match self {
            UnifiedOp::Match(len) | UnifiedOp::Mis(len) | UnifiedOp::Del(len) | UnifiedOp::RefSkip(len)=> *len,
            UnifiedOp::Mismatch(seq) | UnifiedOp::Deletion(seq) => seq.len(),
            _ => 0,
        }
    }
    #[must_use]
    pub fn read_consumed(&self) -> usize {
        match self {
            UnifiedOp::Match(len) | UnifiedOp::Mis(len) | UnifiedOp::Ins(len) => *len,
            UnifiedOp::Mismatch(seq) | UnifiedOp::Insertion(seq) => seq.len(),
            _ => 0,
        }
    }
    pub fn same_variant(&self, other: &UnifiedOp) -> bool {
        match (self, other) {
            (UnifiedOp::Match(l1), UnifiedOp::Match(l2)) | (UnifiedOp::Mis(l1), UnifiedOp::Mis(l2)) | (UnifiedOp::Ins(l1), UnifiedOp::Ins(l2)) | (UnifiedOp::Del(l1), UnifiedOp::Del(l2)) => l1 == l2,
            (UnifiedOp::Mate, UnifiedOp::Mate) | (UnifiedOp::MateSwitch{..}, UnifiedOp::MateSwitch{..}) | (UnifiedOp::RefSkip(_), UnifiedOp::RefSkip(_)) => true,
            _ => false,
        }
    }
}

pub struct UnifiedOpIterator<'a> {
    ops: Peekable<std::slice::Iter<'a, UnifiedOp>>,
    next_op: Option<UnifiedOp>,
}

impl<'a> UnifiedOpIterator<'a> {
    #[must_use]
    pub fn new(ops: &'a [UnifiedOp]) -> Self {
        Self {
            ops: ops.iter().peekable(),
            next_op: None,
        }
    }

    #[must_use]
    pub fn peek(&mut self) -> Option<&&'a UnifiedOp> {
        self.ops.peek()
    }
    pub fn set_op(&'a mut self, new_op: UnifiedOp) {
        if new_op.ref_consumed() == 0 && new_op.read_consumed() == 0 {
            return;
        }
        self.next_op = Some(new_op);
        if let Some(op) = self.ops.peek_mut() {
            *op = &self.next_op.as_ref().unwrap();
        }
    }
}

impl<'a> Iterator for UnifiedOpIterator<'a> {
    type Item = &'a UnifiedOp;

    fn next(&mut self) -> Option<Self::Item> {
        self.ops.next()
    }
}

pub enum AlnCmpOp {
    Equal,
    UnEqual(UnifiedOp, UnifiedOp, u8),
    Left(UnifiedOp, u8),
    Right(UnifiedOp, u8),
}

pub fn lift_alignment_ops(
    cigar: CigarStringView,
    md: &str,
) -> Result<Vec<UnifiedOp>, AlignmentError> {
    let mut unified_ops = Vec::new();
    let cigar_iter = cigar.iter().peekable();
    let mut md_iter: Box<dyn Iterator<Item = RichMdOp>> =
        Box::new(RichMdOpIterator::new(md).peekable());

    for cigar_op in cigar_iter {
        match *cigar_op {
            // --- CIGAR Ops that consume CIGAR and MD ---
            Cigar::Match(mut len) | Cigar::Equal(mut len) | Cigar::Diff(mut len) => {
                while len > 0 {
                    let Some(md_op) = (*md_iter).next() else {
                        return Err(AlignmentError::MdCigarMismatch);
                    };

                    match md_op {
                        RichMdOp::Match(md_len) => {
                            let op_len = usize::min(len, md_len);
                            unified_ops.push(UnifiedOp::Match(op_len));

                            // Put back the remainder if MD was longer
                            if md_len > op_len {
                                md_iter = Box::new(
                                    once(RichMdOp::Match(md_len - op_len))
                                        .chain(md_iter)
                                        .peekable(),
                                );
                            }
                            len -= op_len;
                        }
                        RichMdOp::Mismatch(base) => {
                            unified_ops.push(UnifiedOp::Mismatch(base));
                            len -= 1;
                        }
                        RichMdOp::Deletion(bases) => {
                            // This is complex. A CIGAR 'M' can contain an MD 'D'.
                            // This implies the CIGAR was wrong, but we follow MD.
                            unified_ops.push(UnifiedOp::Deletion(bases));
                            len -= 1; // Or bases.len()? This implies a spec violation.
                            // Let's assume MD matches CIGAR.
                        }
                    }
                }
            }

            // --- CIGAR Ops that consume MD only ---
            Cigar::Del(len) => {
                let Some(RichMdOp::Deletion(bases)) = (*md_iter).next() else {
                    return Err(AlignmentError::MdCigarMismatch);
                };
                if bases.len() != len as usize {
                    return Err(AlignmentError::MdCigarMismatch);
                }
                unified_ops.push(UnifiedOp::Deletion(bases));
            }

            // --- CIGAR Ops that consume CIGAR only ---
            Cigar::Ins(len) => unified_ops.push(UnifiedOp::Insertion(len)),
            Cigar::SoftClip(len) => unified_ops.push(UnifiedOp::SoftClip(len)),
            Cigar::RefSkip(len) => unified_ops.push(UnifiedOp::RefSkip(len)),

            // --- Ignored Ops ---
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    // (You will need a `put_back` method on your `RichMdOpIterator`)
    // (This code is a sketch; the logic is complex but this is the pattern)
    Ok(unified_ops)
}

#[derive(Debug, Clone)]
pub enum AlignmentOp {
    Match,
    Mismatch,
    Insertion,
    Deletion,
    SoftClip,
    #[allow(dead_code)]
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
    pub fn score(self, q: u8, indel_gap: &mut Option<bool>) -> f64 {
        let log_likelihood_mismatch = LOG_LIKELIHOOD_MISMATCH.get().unwrap();
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
                let config = CONFIG.get().unwrap();
                if *indel_gap == Some(true) {
                    config.gap_extend
                } else {
                    *indel_gap = Some(true);
                    config.gap_open + config.gap_extend
                }
            }
            AlignmentOp::Deletion => {
                let config = CONFIG.get().unwrap();
                if *indel_gap == Some(false) {
                    config.gap_extend
                } else {
                    *indel_gap = Some(false);
                    config.gap_open + config.gap_extend
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
