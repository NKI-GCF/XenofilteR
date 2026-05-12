use crate::alignment::MdOpIteratorError;
use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::errors::Error as HtslibError;
use smallvec::SmallVec;
use std::str::Chars;


#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MdOp {
    Match(u32),
    Mismatch(u8),
    Deletion(SmallVec<[u8; 4]>), // slightly larger stack buffer
}

#[cfg_attr(test, derive(Debug))]
pub struct MdOpIterator<'a> {
    chars: Chars<'a>,
    peeked: Option<char>,
}

impl<'a> MdOpIterator<'a> {
    pub fn new(rec: &'a Record) -> Result<Self, MdOpIteratorError> {
        let chars = match rec.aux(b"MD") {
            Ok(Aux::String(md)) => md.chars(),
            Ok(_) => return Err(MdOpIteratorError::BadMdTag),
            Err(HtslibError::BamAuxTagNotFound) => "".chars(),
            Err(e) => return Err(MdOpIteratorError::Aux(e)),
        };
        Ok(Self { chars, peeked: None })
    }

    pub fn empty() -> Self {
        Self { chars: "".chars(), peeked: None }
    }

    /// Very fast path used in hot `partial_cmp`
    pub fn is_single_operation(&self) -> bool {
        for c in self.chars.clone() {
            if !c.is_ascii_digit() {
                return false;
            }
        }
        true
    }
}

impl<'a> Iterator for MdOpIterator<'a> {
    type Item = Result<MdOp, MdOpIteratorError>;

    fn next(&mut self) -> Option<Self::Item> {
        let c = self.peeked.take().or_else(|| self.chars.next())?;

        match c {
            'A' | 'C' | 'G' | 'T' | 'N' => Some(Ok(MdOp::Mismatch(c as u8))),

            '^' => {
                let mut deleted = SmallVec::new();
                for ch in self.chars.by_ref() {
                    if matches!(ch, 'A' | 'C' | 'G' | 'T' | 'N') {
                        deleted.push(ch as u8);
                    } else {
                        self.peeked = Some(ch);
                        break;
                    }
                }
                Some(Ok(MdOp::Deletion(deleted)))
            }

            n if n.is_ascii_digit() => {
                let mut num = (n as u32) - b'0' as u32;
                for ch in self.chars.by_ref() {
                    if ch.is_ascii_digit() {
                        num = num * 10 + (ch as u32 - b'0' as u32);
                    } else {
                        self.peeked = Some(ch);
                        break;
                    }
                }
                Some(Ok(MdOp::Match(num)))
            }
            x => Some(Err(MdOpIteratorError::MdParse(x))),
        }
    }
}

#[cfg(test)]
mod tests;
