use crate::alignment::MdOpIteratorError;
use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::errors::Error as HtslibError;
use smallvec::SmallVec;
use std::str::Chars;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MdOp {
    Match(u32),
    Mismatch(u8),
    Deletion(SmallVec<[u8; 1]>),
}

#[cfg_attr(test, derive(Debug))]
pub struct MdOpIterator<'a> {
    chars: Chars<'a>,
    peeked: String,
}

impl<'a> MdOpIterator<'a> {
    pub fn new(rec: &'a Record) -> Result<Self, MdOpIteratorError> {
        let chars = match rec.aux(b"MD") {
            Ok(Aux::String(md)) => md.chars(),
            Ok(_) => return Err(MdOpIteratorError::BadMdTag),
            Err(HtslibError::BamAuxTagNotFound) => "".chars(),
            Err(e) => return Err(MdOpIteratorError::Aux(e)),
        };
        Ok(MdOpIterator {
            chars,
            peeked: String::new(),
        })
    }
    pub fn empty() -> Self {
        MdOpIterator {
            chars: "".chars(),
            peeked: String::new(),
        }
    }
    pub fn is_single_operation(&self) -> bool {
        for c in self.chars.clone() {
            //#[cfg(test)]
            //eprintln!("Checking char: {}", c);
            match c {
                n if n.is_ascii_digit() => continue,
                _ => return false,
            }
        }
        true
    }
}

impl<'a> Iterator for MdOpIterator<'a> {
    type Item = Result<MdOp, MdOpIteratorError>;

    fn next(&mut self) -> Option<Self::Item> {
        let current_char = if let Some(c) = self.peeked.pop() {
            c
        } else {
            self.chars.next()?
        };
        match current_char {
            c @ ('A' | 'C' | 'G' | 'T' | 'N') => Some(Ok(MdOp::Mismatch(c as u8))),
            '^' => {
                let mut deleted_bases: SmallVec<[u8; 1]> = SmallVec::new();
                for c in self.chars.by_ref() {
                    if matches!(c, 'A' | 'C' | 'G' | 'T' | 'N') {
                        deleted_bases.push(c as u8);
                    } else {
                        self.peeked = String::from(c);
                        break;
                    }
                }
                Some(Ok(MdOp::Deletion(deleted_bases)))
            }
            n if n.is_ascii_digit() => {
                let mut num = n as u32 - '0' as u32;
                for c in self.chars.by_ref() {
                    if c.is_ascii_digit() {
                        num = (c as u32 - '0' as u32) + (num * 10);
                    } else {
                        self.peeked = String::from(c);
                        break;
                    }
                }
                Some(Ok(MdOp::Match(num)))
            }
            x => Some(Err(MdOpIteratorError::MdParse(x))),
        }
    }
}

impl DoubleEndedIterator for MdOpIterator<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let current_char = if let Some(c) = self.peeked.pop() {
            c
        } else {
            self.chars.next_back()?
        };

        match current_char {
            n if n.is_ascii_digit() => {
                let mut num = n as u32 - '0' as u32;
                let mut power_of_10 = 10;
                for c in self.chars.by_ref().rev() {
                    if c.is_ascii_digit() {
                        num += (c as u32 - '0' as u32) * power_of_10;
                        power_of_10 *= 10;
                    } else {
                        self.peeked = String::from(c);
                        break;
                    }
                }
                Some(Ok(MdOp::Match(num)))
            }

            m @ ('A' | 'C' | 'G' | 'T' | 'N') => {
                for c in self.chars.by_ref().rev() {
                    match c {
                        '^' => {
                            self.peeked.push(m);
                            return Some(Ok(MdOp::Deletion(
                                self.peeked.drain(..).map(|b| b as u8).collect(),
                            )));
                        }
                        'A' | 'C' | 'G' | 'T' | 'N' => self.peeked.insert(0, c),
                        _ => {
                            self.peeked.insert(0, c);
                            break;
                        }
                    }
                }
                Some(Ok(MdOp::Mismatch(m as u8)))
            }
            x => Some(Err(MdOpIteratorError::MdParse(x))),
        }
    }
}

#[cfg(test)]
mod tests;
