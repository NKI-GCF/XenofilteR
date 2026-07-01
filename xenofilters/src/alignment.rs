mod fragment;
mod fragment_state;
pub(crate) mod mate_kind;
mod md_cig_flags;
pub(crate) mod ops;
pub(crate) mod pre_assess;
mod read_profile; // private — used only by pre_assess
mod variant_window;

pub(crate) use fragment::{Fragment, SimpleRec};
pub(crate) use fragment_state::FragmentState;
pub(crate) use mate_kind::{mate_slot, segment_id, MateClassifiable, MateKind};
pub(crate) use md_cig_flags::MdCigFlags;
pub(crate) use ops::{BaseOp, ScoreOpIter};
pub(crate) use pre_assess::{pre_assess_alignments, pre_assess_scoring_records, PreAssessResult};
pub(crate) use variant_window::{align_alt_to_read, weighted_ref_score, VariantWindow};

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::Record;

pub(crate) fn stringify_record<R: Record + PartialEq>(rec: &R) -> String {
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
    s.push_str(&format!(
        "\treverse:{}",
        rec.flags().unwrap().is_reverse_complemented()
    ));
    //s.push_str(&format!("\tseq:{}", rec.seq().as_bytes().iter().map(|&b| b as char).collect::<String>()));
    //s.push_str(&format!("\tqual:{}", rec.qual().iter().map(|&q| (q + 33) as char).collect::<String>()));

    s
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::Error;
    pub(crate) use ops::tests::*;

    #[test]
    fn test_stringify_record_includes_qname_cigar_md_and_orientation() -> Result<(), Error> {
        let rec = create_record(b"read1", "5M", &[], &[30; 5], "5", false)?;
        let s = stringify_record(&rec);
        assert!(s.contains("read1"));
        assert!(s.contains("5M"));
        assert!(s.contains("MD:Z:5"));
        assert!(s.contains("reverse:false"));
        Ok(())
    }
}
