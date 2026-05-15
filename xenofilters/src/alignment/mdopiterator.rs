use crate::alignment::MdOpIteratorError;
use noodles::bam::record::Record;
use noodles::sam::alignment::record::data::field::Value;
use smallvec::SmallVec;
use std::slice::Iter;


#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum MdOp {
    Match(usize),
    Mismatch(u8),
    Deletion(SmallVec<[u8; 4]>), // slightly larger stack buffer
}

#[cfg_attr(test, derive(Debug))]
pub(crate) struct MdOpIterator<'a> {
    bytes: Option<Iter<'a, u8>>,
    peeked: Option<u8>,
}

impl<'a> MdOpIterator<'a> {
    pub(crate) fn new(rec: &'a Record) -> Result<Self, MdOpIteratorError> {
        let bytes = match rec.data().get(b"MD") {
            Some(Ok(Value::String(md))) => {
                let r: &[u8] = md.as_ref();
                Some(r.iter())
            },
            Some(Ok(_)) => return Err(MdOpIteratorError::BadMdTag),
            Some(Err(e)) => return Err(MdOpIteratorError::MdError(e)),
            None => None,
        };
        Ok(Self { bytes, peeked: None })
    }

    pub(crate) fn empty() -> Self {
        Self { bytes: None, peeked: None }
    }

    pub(crate) fn is_single_operation(&self) -> bool {
        if let Some (bytes) = &self.bytes {
            for b in bytes.clone() {
                if !b.is_ascii_digit() {
                    return false;
                }
            }
        }
        true
    }
}

impl<'a> Iterator for MdOpIterator<'a> {
    type Item = Result<MdOp, MdOpIteratorError>;

    fn next(&mut self) -> Option<Self::Item> {
        let b = self.peeked.take().or_else(|| self.bytes.as_mut()?.next().copied())?;
            //or_else(|| self.bytes.next())?;

        match b {
            b'A' | b'C' | b'G' | b'T' | b'N' => Some(Ok(MdOp::Mismatch(b))),

            b'^' => {
                let mut deleted = SmallVec::new();
                if let Some(bytes) = self.bytes.as_mut() {
                    for byte in bytes.by_ref() {
                        if matches!(byte, b'A' | b'C' | b'G' | b'T' | b'N') {
                            deleted.push(*byte);
                        } else {
                            self.peeked = Some(*byte);
                            break;
                        }
                    }
                }
                Some(Ok(MdOp::Deletion(deleted)))
            }

            n if n.is_ascii_digit() => {
                let mut num = (n - b'0') as usize;
                if let Some(bytes) = self.bytes.as_mut() {
                    for byte in bytes.by_ref() {
                        if byte.is_ascii_digit() {
                            num = num * 10 + (byte - b'0') as usize;
                        } else {
                            self.peeked = Some(*byte);
                            break;
                        }
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
