pub(crate) mod common;
mod exploratory;

use super::*;
#[cfg(test)]
pub(crate) mod property;
#[cfg(test)]
pub(crate) use aln_stream::tests::*;

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
use noodles::sam::Header;
use std::iter::repeat;
use crate::region::load::{load_ambiguous_regions_memory, load_distinct_variants_memory};
use crate::variant::name_to_id::header_name_to_id;

fn create_cigar(cigar: &str) -> Result<Cigar, Error> {
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

fn read_len_from_cigar(cigar: &str) -> usize {
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

pub fn create_record(
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

// Kills mutations in `header_name_to_id` (HashMap::new(), HashMap::from_iter, etc.)
#[test]
fn test_header_name_to_id() {
    // Construct a realistic SAM header to ensure iteration and indexing are correct
    let header_str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200";
    let header: Header = header_str.parse().expect("Failed to parse SAM header");

    let map = header_name_to_id(&header);

    assert_eq!(map.len(), 2);
    assert_eq!(map.get("chr1"), Some(&0));
    assert_eq!(map.get("chr2"), Some(&1));
}

#[test]
fn test_load_ambiguous_regions_ignores_empty_strings() {
    let name_to_id = HashMap::new();
    let specs = vec!["".to_string(), "".to_string()];
    let result = load_ambiguous_regions_memory(&specs, &name_to_id).unwrap();
    assert!(result[0].is_none());
    assert!(result[1].is_none());

    let specs_single = vec!["".to_string()];
    let result_single = load_ambiguous_regions_memory(&specs_single, &name_to_id).unwrap();
    assert!(result_single[0].is_none());
    assert!(result_single[1].is_none());
}

#[test]
fn test_load_distinct_variants_ignores_empty_strings() {
    let name_to_id = HashMap::new();
    let specs = vec!["".to_string(), "".to_string()];
    let result = load_distinct_variants_memory(&specs, &name_to_id).unwrap();
    assert!(result[0].is_none());
    assert!(result[1].is_none());
}

#[test]
fn test_namesorted_sequential_single_alignment() {
    use crate::config::{CommonArgs, NamesortedArgs};
    use crate::config::args::{ChimericArgs, IoArgs, ParallelArgs};

    let args = NamesortedArgs {
        common: CommonArgs {
            io: IoArgs { alignment: vec!["tests/fixtures/dummy1.bam".into()], ..Default::default() },
            ..Default::default()
        },
        parallel: ParallelArgs { threads: 1, score_threads: 1 },
        chimeric: ChimericArgs::default(),
    };
    let _ = run_namesorted(args);
}

#[test]
fn test_namesorted_parallel_dual_alignment() {
    use crate::config::{CommonArgs, NamesortedArgs};
    use crate::config::args::{ChimericArgs, IoArgs, ParallelArgs};

    let args = NamesortedArgs {
        common: CommonArgs {
            io: IoArgs {
                alignment: vec![
                    "tests/fixtures/dummy1.bam".into(),
                    "tests/fixtures/dummy2.bam".into(),
                ],
                ..Default::default()
            },
            ..Default::default()
        },
        parallel: ParallelArgs { threads: 1, score_threads: 2 },
        chimeric: ChimericArgs::default(),
    };
    let _ = run_namesorted(args);
}

#[test]
fn test_get_log_level() {
    assert_eq!(get_log_level(0), "warn");
    assert_eq!(get_log_level(1), "info");
    assert_eq!(get_log_level(2), "debug");
    assert_eq!(get_log_level(5), "debug");
}
