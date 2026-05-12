use crate::alignment::AlignmentError;
use crate::alignment::{UnifiedOp, UnifiedOpIterator};
use crate::tests::read_len_from_cigar;
use anyhow::Result;
use rust_htslib::bam::record::{Aux, Cigar, CigarString, Record};

impl UnifiedOp {
    #[must_use]
    pub(crate) fn len(&self) -> u32 {
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
    pub(crate) fn seq_len(&self) -> u32 {
        self.cigar_iter
            .as_slice()
            .iter()
            .map(|c| match c {
                Cigar::Match(len)
                | Cigar::Equal(len)
                | Cigar::Diff(len)
                | Cigar::Ins(len)
                | Cigar::SoftClip(len) => *len,
                Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => 0,
            })
            .sum()
    }
}

pub(crate) fn create_record(
    qname: &[u8],
    cig_str: &str,
    seq: &[u8],
    qual: &[u8],
    md: &str,
    is_rev: bool,
) -> Result<rust_htslib::bam::Record> {
    let mut record = Record::new();
    let is_unmapped = cig_str == "*";

    let cigar = if is_unmapped {
        CigarString::try_from("")?
    } else {
        CigarString::try_from(cig_str)?
    };
    let read_len = if is_unmapped {
        0
    } else {
        read_len_from_cigar(cig_str)
    };
    let final_seq = if seq.is_empty() {
        vec![b'A'; read_len]
    } else {
        seq.to_vec()
    };

    let final_qual = if qual.is_empty() {
        vec![30u8; read_len]
    } else {
        qual.to_vec()
    };

    record.set(qname, Some(&cigar), &final_seq, &final_qual);

    if is_unmapped {
        record.set_unmapped();
    } else {
        // HTSlib requires the MD tag to be set manually for tests
        record.push_aux(b"MD", Aux::String(md))?;
    }

    if is_rev {
        record.set_reverse();
    }

    /*#[cfg(test)]
    eprintln!(
        "Created record: qname={:?}, cigar={:?}, is_rev={}, md={}",
        std::str::from_utf8(record.qname())?,
        record.cigar(),
        is_rev,
        md
    );*/
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
    // "50M10I40M"
    let rec = create_record(b"read1", "50M10I40M", &[], &[30; 100], "90", false)?;
    let iter = UnifiedOpIterator::new(&rec).unwrap();

    let total_len: u32 = iter.map(|r| r.unwrap().len()).sum();
    // Sum of Match/Mis/Ins/Del should match the logical alignment footprint
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
