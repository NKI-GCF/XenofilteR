use crate::alignment::{BaseOp, MdCigFlags, ScoreOpIter};
use noodles::core::Position;
use noodles::sam::alignment::{
    record::data::field::Tag,
    record_buf::{data::field::Value, Data},
};
use noodles::sam::alignment::{
    record::{
        cigar::{op::Kind, Op},
        Flags,
    },
    record_buf::{Cigar, QualityScores, RecordBuf, Sequence},
};
use std::iter::repeat;
use crate::Error;

pub(crate) fn create_cigar(cigar: &str) -> Result<Cigar, Error> {
    let mut ops = Vec::new();
    let mut num = 0;
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as u8 - b'0') as usize;
        } else {
            let kind = match c {
                'M' => Kind::Match,
                'I' => Kind::Insertion,
                'D' => Kind::Deletion,
                'N' => Kind::Skip,
                'S' => Kind::SoftClip,
                'H' => Kind::HardClip,
                'P' => Kind::Pad,
                '=' => Kind::SequenceMatch,
                'X' => Kind::SequenceMismatch,
                _ => return Err(Error::UnknownCigarOp(c as u32)),
            };
            ops.push(Op::new(kind, num));
            num = 0;
        }
    }
    Ok(Cigar::from(ops))
}

pub(crate) fn read_len_from_cigar(cigar: &str) -> usize {
    cigar
        .chars()
        .fold((0, 0), |(len, acc), c| {
            if c.is_ascii_digit() {
                (len, acc * 10 + (c as u8 - b'0') as usize)
            } else {
                let consumes_read = matches!(c, 'M' | 'I' | 'S' | '=' | 'X');
                (if consumes_read { len + acc } else { len }, 0)
            }
        })
        .0
}

pub(crate) fn create_record(
    qname: &[u8],
    cig_str: &str,
    seq: &[u8],
    qual: &[u8],
    md: &str,
    is_rev: bool,
) -> Result<RecordBuf, Error> {
    let mut record = RecordBuf::default();
    let is_unmapped = cig_str.is_empty() || cig_str == "*";

    *record.name_mut() = Some(qname.into());
    let read_len = if is_unmapped {
        *record.flags_mut() = Flags::from_bits(4).unwrap(); // unmapped
        0
    } else {
        if is_rev {
            *record.flags_mut() = Flags::from_bits(16).unwrap();
        } else {
            *record.flags_mut() = Flags::empty();
        }
        let data: Data = [(Tag::MISMATCHED_POSITIONS, Value::from(md))]
            .into_iter()
            .collect();
        *record.data_mut() = data;
        *record.cigar_mut() = create_cigar(cig_str)?;
        *record.alignment_start_mut() = Some(Position::MIN);
        *record.reference_sequence_id_mut() = Some(0);
        read_len_from_cigar(cig_str)
    };
    *record.sequence_mut() = if seq.is_empty() {
        Sequence::from(repeat(b'A').take(read_len).collect::<Vec<u8>>())
    } else {
        Sequence::from(seq)
    };

    *record.quality_scores_mut() = if qual.is_empty() {
        QualityScores::from_iter(repeat(30u8).take(read_len).collect::<Vec<u8>>())
    } else {
        QualityScores::from_iter(qual.iter().cloned())
    };

    Ok(record)
}

fn op_repr(op: &BaseOp) -> String {
    match op {
        BaseOp::Match => "M".into(),
        BaseOp::Mis => "X".into(),
        BaseOp::Del(n) => format!("D{n}"),
        BaseOp::Ins(n) => format!("I{n}"),
        BaseOp::Clip(n) => format!("C{n}"),
        BaseOp::RefSkip(n) => format!("N{n}"),
    }
}

fn ops_for(cigar: &str, md: &str) -> Result<Vec<String>, Error> {
    let rec = create_record(b"r", cigar, &[], &[], md, false)?;
    let flags = rec.flags();
    let mcf = MdCigFlags::try_from_record(&rec, &flags)?;
    let ops: Result<Vec<BaseOp>, _> = ScoreOpIter::new(&mcf).collect();
    Ok(ops
        .map_err(|e| anyhow::anyhow!("{e}"))?
        .iter()
        .map(op_repr)
        .collect())
}

#[test]
fn test_pure_match() -> Result<(), Error> {
    assert_eq!(ops_for("5M", "5")?, vec!["M", "M", "M", "M", "M"]);
    Ok(())
}

#[test]
fn test_single_mismatch() -> Result<(), Error> {
    assert_eq!(ops_for("5M", "2A2")?, vec!["M", "M", "X", "M", "M"]);
    Ok(())
}

#[test]
fn test_consecutive_mismatches() -> Result<(), Error> {
    // MD "1AC2": 1 match, 2 consecutive mismatches, 2 matches = 5 bases
    assert_eq!(ops_for("5M", "1AC2")?, vec!["M", "X", "X", "M", "M"]);
    Ok(())
}

#[test]
fn test_insertion_not_consumed_by_md_and_grouped() -> Result<(), Error> {
    assert_eq!(ops_for("2M2I2M", "4")?, vec!["M", "M", "I2", "M", "M"]);
    Ok(())
}

#[test]
fn test_deletion_consumes_caret_block() -> Result<(), Error> {
    assert_eq!(ops_for("2M3D2M", "2^AAA2")?, vec!["M", "M", "D3", "M", "M"]);
    Ok(())
}

#[test]
fn test_soft_clip_is_single_grouped_op() -> Result<(), Error> {
    assert_eq!(ops_for("3S5M", "5")?, vec!["C3", "M", "M", "M", "M", "M"]);
    Ok(())
}

#[test]
fn test_ref_skip_not_consumed_by_md() -> Result<(), Error> {
    assert_eq!(ops_for("2M3N2M", "4")?, vec!["M", "M", "N3", "M", "M"]);
    Ok(())
}

#[test]
fn test_hard_clip_and_pad_are_invisible_to_iterator() -> Result<(), Error> {
    assert_eq!(ops_for("2H5M2H", "5")?, vec!["M", "M", "M", "M", "M"]);
    assert_eq!(ops_for("2M2P2M", "4")?, vec!["M", "M", "M", "M"]);
    Ok(())
}

#[test]
fn test_sequence_match_mismatch_kinds_route_through_md() -> Result<(), Error> {
    assert_eq!(ops_for("5=", "5")?, vec!["M", "M", "M", "M", "M"]);
    Ok(())
}

#[test]
fn test_multidigit_md_run_length() -> Result<(), Error> {
    assert_eq!(ops_for("12M", "12")?.len(), 12);
    Ok(())
}

#[test]
fn test_deletion_length_mismatch_is_error() {
    // CIGAR deletes 3 bases, MD caret block only has 2 — must error, not silently desync.
    assert!(ops_for("2M3D2M", "2^AA2").is_err());
}

#[test]
fn test_invalid_md_character_is_error() {
    assert!(ops_for("4M", "2x2").is_err());
}
