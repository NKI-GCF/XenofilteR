use anyhow::Result;
use noodles::sam::alignment::{
    record::{Flags, cigar::{op::Kind, Op}}, record_buf::{Cigar, RecordBuf, Sequence, QualityScores},
};
use std::iter::repeat;
use noodles::sam::{
    alignment::{
        record::data::field::Tag,
        record_buf::{data::field::Value, Data},
    },
};
use noodles::core::Position;

pub(crate) fn create_cigar(cigar: &str) -> Result<Cigar> {
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
                _ => return Err(anyhow::anyhow!("Invalid CIGAR character: {c}")),
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
) -> Result<RecordBuf> {
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

