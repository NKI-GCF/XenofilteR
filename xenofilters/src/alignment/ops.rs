use crate::alignment::AlignmentError;
use crate::{CONFIG, LOG_LIKELIHOOD_MATCH, LOG_LIKELIHOOD_MISMATCH};
use rust_htslib::bam::record::{Cigar, CigarStringView};

use std::iter::{Peekable, once};
use std::str::Chars;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RichMdOp {
    Match(u32),
    Mismatch(u8),      // A single base mismatch
    Deletion(Vec<u8>), // One or more deleted bases
    Translocate {
        // Your new op, parsed from a custom MD tag
        new_tid: u32,
        new_pos: i64,
        new_strand: bool,
    },
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
        let &c = self.chars.peek()?;

        if c.is_ascii_digit() {
            // 1. Handle Match
            let mut num = 0;
            while let Some(&digit) = self.chars.peek().filter(|c| c.is_ascii_digit()) {
                self.chars.next(); // Consume digit
                num = (num * 10) + (digit as u32 - '0' as u32);
            }
            return Some(RichMdOp::Match(num));
        } else if c == '^' {
            // 2. Handle Deletion or Translocate
            self.chars.next(); // Consume '^'

            // Check for *your* custom translocation tag
            // e.g., ^TID,POS,STRAND
            if self.chars.peek() == Some(&'T') {
                // This is a custom parser for your op.
                // ... logic to parse TID,POS,STRAND ...
                // return Some(RichMdOp::Translocate { ... });
                unimplemented!("Translocate parsing not yet implemented");
            }

            // Standard Deletion
            let mut deleted_bases = Vec::new();
            while let Some(&base) = self.chars.peek().filter(|c| c.is_ascii_alphabetic()) {
                self.chars.next(); // Consume base
                deleted_bases.push(base as u8);
            }
            return Some(RichMdOp::Deletion(deleted_bases));
        } else if c.is_ascii_alphabetic() {
            // 3. Handle Mismatch
            self.chars.next(); // Consume base
            return Some(RichMdOp::Mismatch(c as u8));
        }

        // Should be unreachable
        self.chars.next(); // Consume unknown char
        None
    }
}

#[derive(Debug, Clone)]
pub enum UnifiedOp {
    /// A pure reference match
    Match(u32),

    /// A mismatch found in the MD tag
    Mismatch(u8), // The reference base

    /// An insertion from CIGAR (MD-agnostic)
    Insertion(u32),

    /// A deletion from CIGAR, with ref bases from MD
    Deletion(Vec<u8>),

    /// A softclip from CIGAR
    SoftClip(u32),

    /// A refskip from CIGAR
    RefSkip(u32),

    /// A translocation (from your stitched CIGAR idea)
    Translocate {
        new_tid: i32,
        new_pos: i64,
        new_strand: bool,
        gap_score: f64,
    },
}

pub fn lift_alignment_ops(
    cigar: CigarStringView,
    md: &str,
) -> Result<Vec<UnifiedOp>, AlignmentError> {
    let mut unified_ops = Vec::new();
    let mut cigar_iter = cigar.iter().peekable();
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
                            let op_len = u32::min(len, md_len);
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
                        RichMdOp::Translocate { .. } => {
                            return Err(AlignmentError::UnexpectedTranslocate);
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
