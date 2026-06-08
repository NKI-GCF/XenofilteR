mod errors;
mod ops;
mod fragment;
mod fragment_state;

pub(crate) use errors::AlignmentError;
pub(crate) use ops::{BaseOp, ScoreOpIter};
pub(crate) use fragment::Fragment;
pub(crate) use fragment_state::{FragmentState, MdCigFlags};

use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::{Tag, Value};

pub(crate) fn stringify_record<R: Record>(rec: &R) -> String {
    let qname = rec.name();

    let mut cigar = String::new();
    for c in rec.cigar().as_ref().iter() {
        let c = c.unwrap();
        match c.kind() {
            Kind::Match => cigar.push_str(&format!("{}M", c.len())),
            Kind::Insertion => cigar.push_str(&format!("{}I", c.len())),
            Kind::Deletion => cigar.push_str(&format!("{}D", c.len())),
            Kind::SoftClip => cigar.push_str(&format!("{}S", c.len())),
            Kind::HardClip => cigar.push_str(&format!("{}H", c.len())),
            Kind::Skip => cigar.push_str(&format!("{}N", c.len())),
            Kind::Pad => cigar.push_str(&format!("{}P", c.len())),
            Kind::SequenceMatch => cigar.push_str(&format!("{}=", c.len())),
            Kind::SequenceMismatch => cigar.push_str(&format!("{}X", c.len())),
        }
    }
    let mut s = format!("{qname:?}\t{cigar:?}");
    if let Some(Ok(Value::String(md))) = rec.data().as_ref().get(&Tag::MISMATCHED_POSITIONS) {
        s.push_str(&format!("\tMD:Z:{md}"));
    }
    s.push_str(&format!("\treverse:{}", rec.flags().unwrap().is_reverse_complemented()));
    //s.push_str(&format!("\tseq:{}", rec.seq().as_bytes().iter().map(|&b| b as char).collect::<String>()));
    //s.push_str(&format!("\tqual:{}", rec.qual().iter().map(|&q| (q + 33) as char).collect::<String>()));

    s
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    pub(crate) use ops::tests::*;
    pub(crate) use prepared::tests::*;
}
