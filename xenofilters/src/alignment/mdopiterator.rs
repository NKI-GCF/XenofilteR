use rust_htslib::bam::record::{Aux, Record};
use smallvec::SmallVec;
use std::str::Chars;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum MdOpIteratorError {
    #[error("No MD tag found")]
    NoMdTag,

    #[error("Aux error {0}")]
    AuxError(String),

    #[error("MD parsing error: invalid character '{0}'")]
    MdParseError(char),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MdOp {
    Match(u32),
    Mismatch(u8),
    Deletion(SmallVec<[u8; 1]>),
}

pub struct MdOpIterator<'a> {
    chars: Chars<'a>,
    peeked: String,
}

impl<'a> MdOpIterator<'a> {
    pub fn new(rec: &'a Record) -> Result<Self, MdOpIteratorError> {
        let chars = match rec.aux(b"MD") {
            Ok(Aux::String(md)) => md.chars(),
            Ok(_) => return Err(MdOpIteratorError::NoMdTag),
            Err(e) => return Err(MdOpIteratorError::AuxError(e.to_string())),
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
            x => Some(Err(MdOpIteratorError::MdParseError(x))),
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
            x => Some(Err(MdOpIteratorError::MdParseError(x))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;

    fn process_md_forward(md_string: &str) -> Vec<MdOp> {
        let md_iter = MdOpIterator {
            chars: md_string.chars(),
            peeked: String::new(),
        };
        md_iter.map(|r| r.unwrap()).collect()
    }
    fn process_md_reverse(md_string: &str) -> Vec<MdOp> {
        let md_iter = MdOpIterator {
            chars: md_string.chars(),
            peeked: String::new(),
        };
        md_iter.rev().map(|r| r.unwrap()).collect()
    }

    #[test]
    fn test_forward_mdop_10a5() {
        let md_ops = process_md_forward("10A5");
        assert_eq!(
            md_ops,
            vec![MdOp::Match(10), MdOp::Mismatch(b'A'), MdOp::Match(5),]
        );
    }

    #[test]
    fn test_forward_mdop_tga() {
        let md_ops = process_md_forward("TGA");
        assert_eq!(
            md_ops,
            vec![
                MdOp::Mismatch(b'T'),
                MdOp::Mismatch(b'G'),
                MdOp::Mismatch(b'A'),
            ]
        );
    }

    #[test]
    fn test_forward_mdop_adt20g() {
        let md_ops = process_md_forward("A^T20G");
        assert_eq!(
            md_ops,
            vec![
                MdOp::Mismatch(b'A'),
                MdOp::Deletion(smallvec![b'T']),
                MdOp::Match(20),
                MdOp::Mismatch(b'G'),
            ]
        );
    }

    #[test]
    fn test_forward_mdop_5datc3() {
        let md_ops = process_md_forward("5^ATC3");
        assert_eq!(
            md_ops,
            vec![
                MdOp::Match(5),
                MdOp::Deletion(smallvec![b'A', b'T', b'C']),
                MdOp::Match(3),
            ]
        );
    }

    #[test]
    fn test_forward_mdop_1dn5() {
        let md_ops = process_md_forward("1^N5");
        assert_eq!(
            md_ops,
            vec![
                MdOp::Match(1),
                MdOp::Deletion(smallvec![b'N']),
                MdOp::Match(5),
            ]
        );
    }

    #[test]
    fn test_forward_mdop_999g1() {
        let md_ops = process_md_forward("999G1");
        assert_eq!(
            md_ops,
            vec![MdOp::Match(999), MdOp::Mismatch(b'G'), MdOp::Match(1),]
        );
    }

    #[test]
    #[should_panic]
    fn test_forward_mdop_invalid_char() {
        let _md_ops = process_md_forward("10A5X");
    }

    #[test]
    #[should_panic]
    fn test_reverse_mdop_invalid_char() {
        let _md_ops = process_md_reverse("10A5X");
    }

    // --- Reverse Tests ---
    #[test]
    fn test_reverse_mdop_10a5() {
        let md_ops = process_md_reverse("10A5");
        assert_eq!(
            md_ops,
            vec![MdOp::Match(5), MdOp::Mismatch(b'A'), MdOp::Match(10),]
        );
    }

    #[test]
    fn test_reverse_mdop_tga() {
        let md_ops = process_md_reverse("TGA");
        assert_eq!(
            md_ops,
            vec![
                MdOp::Mismatch(b'A'),
                MdOp::Mismatch(b'G'),
                MdOp::Mismatch(b'T'),
            ]
        );
    }

    #[test]
    fn test_reverse_mdop_5datc3() {
        let md_ops = process_md_reverse("5^ATC3");
        assert_eq!(
            md_ops,
            vec![
                MdOp::Match(3),
                MdOp::Deletion(smallvec![b'A', b'T', b'C']),
                MdOp::Match(5),
            ]
        );
    }

    #[test]
    fn test_reverse_mdop_adt20g() {
        let md_ops = process_md_reverse("A^T20G");
        assert_eq!(
            md_ops,
            vec![
                MdOp::Mismatch(b'G'),
                MdOp::Match(20),
                MdOp::Deletion(smallvec![b'T']),
                MdOp::Mismatch(b'A'),
            ]
        );
    }

    #[test]
    fn test_reverse_mdop_100() {
        let md_ops = process_md_reverse("100");
        assert_eq!(md_ops, vec![MdOp::Match(100),]);
    }
}
