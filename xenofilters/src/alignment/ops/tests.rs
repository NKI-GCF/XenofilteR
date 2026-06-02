use crate::alignment::AlignmentError;
use crate::alignment::{UnifiedOp, UnifiedOpIterator};
use crate::tests::{read_len_from_cigar, create_cigar};
use anyhow::Result;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Value;
use noodles::sam::Header;
use noodles::sam::alignment::record_buf::{Sequence, QualityScores};

impl UnifiedOp {
    #[must_use]
    pub(crate) fn len(&self) -> usize {
        match self {
            UnifiedOp::Match(len)
            | UnifiedOp::Mis(len)
            | UnifiedOp::Ins(len)
            | UnifiedOp::Del(len)
            | UnifiedOp::RefSkip(len) => *len,
            //UnifiedOp::Mismatch(seq) | UnifiedOp::Insertion(seq) | UnifiedOp::Deletion(seq) => seq.len(),
            _ => 0,
        }
    }
}

impl UnifiedOpIterator<'_> {
    pub(crate) fn seq_len(&self) -> usize {
        self.cigar_len.unwrap()
    }
}

pub(crate) fn create_record(
    qname: &[u8],
    cig_str: &str,
    seq: &[u8],
    qual: &[u8],
    md: &str,
    is_rev: bool,
) -> Result<RecordBuf> {
    let mut record = RecordBuf::default();
    let is_unmapped = cig_str.is_empty();

    *record.name_mut() = Some(qname.into());
    let read_len = if is_unmapped {
        record.set_flags(4);
        0
    } else {
        if is_rev {
            record.set_flags(16);
        } else {
            record.set_flags(0);
        }
        record.push_aux(b"MD", Aux::String(md))?;
        read_len_from_cigar(cig_str)
    };
    *record.cigar_mut() = create_cigar(cig_str)?;
    *record.sequence_mut() = if seq.is_empty() {
        Sequence::from(&[b'A'; 100][..read_len])
    } else {
        Sequence::from(seq)
    };

    *record.quality_scores_mut() = if qual.is_empty() {
        QualityScores::from(&[30u8; read_len][..read_len])
    } else {
        QualityScores::from(qual)
    };

    Ok(record)
}

#[test]
fn simple_match_cig_10m_md_10() -> Result<()> {
    let rec = create_record(b"read1", "10M", &[], &[30; 10], "10", false)?;
    let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
    let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
    assert_eq!(ops, vec![UnifiedOp::Match(10)]);

    let rec = create_record(b"read2", "10M", &[], &[30; 10], "10", true)?;
    let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
    let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
    assert_eq!(ops, vec![UnifiedOp::Match(10)]);
    Ok(())
}

#[test]
fn mismatch_cig_10m_md_5a4() -> Result<()> {
    let rec = create_record(b"read1", "10M", &[], &[30; 10], "5A4", false)?;
    let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
    let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
    assert_eq!(
        ops,
        vec![UnifiedOp::Match(5), UnifiedOp::Mis(1), UnifiedOp::Match(4)]
    );

    let rec = create_record(b"read2", "10M", &[], &[30; 10], "5A4", true)?;
    let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
    let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
    assert_eq!(
        ops,
        vec![UnifiedOp::Match(5), UnifiedOp::Mis(1), UnifiedOp::Match(4)]
    );
    Ok(())
}

#[test]
fn indels_cig_5m2i3m_md_8() -> Result<()> {
    let rec = create_record(b"read3", "5M2I3M", &[], &[30; 10], "8", false)?;
    let uop_iter = UnifiedOpIterator::new(&rec)?;
    let ops: Vec<UnifiedOp> = uop_iter.collect::<Result<Vec<UnifiedOp>, AlignmentError>>()?;
    assert_eq!(
        ops,
        vec![UnifiedOp::Match(5), UnifiedOp::Ins(2), UnifiedOp::Match(3)]
    );
    Ok(())
}

#[test]
fn soft_clip_cig_5s5m_md_5() -> Result<()> {
    let rec = create_record(b"read4", "5S5M", &[], &[30; 10], "5", false)?;
    let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
    let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
    assert_eq!(ops, vec![UnifiedOp::Mis(5), UnifiedOp::Match(5)]);
    Ok(())
}

#[test]
fn deletion_cig_5m3d5m_md_5daaa5() -> Result<()> {
    let rec = create_record(b"read5", "5M3D5M", &[], &[30; 10], "5^AAA5", false)?;
    let uop_iter = UnifiedOpIterator::new(&rec)?;
    let ops: Vec<UnifiedOp> = uop_iter.collect::<Result<Vec<UnifiedOp>, AlignmentError>>()?;
    assert_eq!(
        ops,
        vec![UnifiedOp::Match(5), UnifiedOp::Del(3), UnifiedOp::Match(5)]
    );
    Ok(())
}

#[test]
fn test_length_invariants() -> Result<()> {
    let rec = create_record(b"read1", "50M10I40M", &[], &[30; 100], "90", false)?;
    let iter = UnifiedOpIterator::new(&rec).unwrap();

    let total_len: u32 = iter.map(|r| r.unwrap().len()).sum();
    assert_eq!(total_len, 100);
    Ok(())
}

#[test]
fn test_create_unmapped_record() -> Result<()> {
    let rec = create_record(b"read_unmapped", "*", &[], &[], "100", false)?;
    assert!(rec.is_unmapped());
    let uop_iter = UnifiedOpIterator::new(&rec).unwrap();
    let ops: Vec<UnifiedOp> = uop_iter.map(|r| r.unwrap()).collect();
    assert_eq!(ops, vec![]);
    Ok(())
}
#[test]
fn test_trailing_ins_mismatch() {
    let rec = create_record(b"INS_TEST", "90M10I", &[], &[], "90", false).unwrap();
    let iter = UnifiedOpIterator::new(&rec).unwrap();
    for res in iter {
        assert!(res.is_ok());
    }
}
