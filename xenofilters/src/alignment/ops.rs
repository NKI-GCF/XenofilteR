#[derive(Debug, Clone)]
pub enum AlignmentOp {
    Match(u8),
    Mismatch(u8),
    Insertion(u8),
    Deletion,
    SoftClip(u8),
}

#[derive(Debug, Clone)]
pub enum MdOp {
    Match(u64),
    Mismatch,
    Deletion,
}

pub struct MdOpIterator<'a> {
    chars: std::str::Chars<'a>,
    // State is simpler: just the number being built and a flag.
    num: u64,
    in_deletion: bool,
}

impl<'a> MdOpIterator<'a> {
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
            let c = match self.chars.next() {
                Some(c) => c,
                None => {
                    if self.num > 0 {
                        let match_op = MdOp::Match(self.num);
                        self.num = 0;
                        return Some(match_op);
                    }
                    return None;
                }
            };

            if c.is_ascii_digit() {
                self.num = (self.num * 10) + (c as u64 - '0' as u64);
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

            if self.in_deletion {
                self.in_deletion = false;
                return Some(MdOp::Deletion);
            } else {
                return Some(MdOp::Mismatch);
            }
        }
    }
}
