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

pub fn parse_md(md: &str) -> Vec<MdOp> {
    let mut ops = Vec::new();
    let mut num = 0;
    let mut chars = md.chars();
    let mut in_deletion = false;

    while let Some(c) = chars.next() {
        if c.is_ascii_digit() {
            num = (num * 10) + (c as u64 - '0' as u64);
            in_deletion = false;
        } else if c == '^' {
            in_deletion = true;
        } else {
            if num > 0 {
                ops.push(MdOp::Match(num));
                num = 0;
            }
            if in_deletion {
                ops.push(MdOp::Deletion);
            } else {
                ops.push(MdOp::Mismatch);
            }
        }
    }
    if num > 0 {
        ops.push(MdOp::Match(num));
    }
    ops
}
