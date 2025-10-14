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
    Deletion(Vec<u8>),
}

pub fn parse_md(md: &str) -> Vec<MdOp> {
    let mut ops = Vec::new();
    let mut num = String::new();
    let mut chars = md.chars().peekable();

    while let Some(c) = chars.next() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            if !num.is_empty() {
                ops.push(MdOp::Match(num.parse().unwrap()));
                num.clear();
            }
            if c == '^' {
                let mut del = Vec::new();
                while let Some(&d) = chars.peek() {
                    if d.is_ascii_alphabetic() {
                        del.push(d as u8);
                        chars.next();
                    } else {
                        break;
                    }
                }
                ops.push(MdOp::Deletion(del));
            } else {
                ops.push(MdOp::Mismatch);
            }
        }
    }
    if !num.is_empty() {
        ops.push(MdOp::Match(num.parse().unwrap()));
    }
    ops
}
