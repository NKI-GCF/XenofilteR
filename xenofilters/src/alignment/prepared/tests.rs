use crate::alignment::{
    PreparedAlignmentPair, PreparedAlignmentPairIter, UnifiedOpIterator, stringify_record,
};
use crate::tests::create_record;
use anyhow::Result;
use rust_htslib::bam::record::Record;

fn print_req(i: usize, rec: &Record) {
    eprintln!("{i}:{}", stringify_record(rec));
}

pub fn read_len_from_cigar(cigar: &str) -> usize {
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

pub fn make_prepared_pair(
    cigar1: &str,
    md1: &str,
    cigar2: &str,
    md2: &str,
) -> Result<PreparedAlignmentPair<'static>> {
    let len1 = read_len_from_cigar(cigar1);
    let rec1 = Box::leak(Box::new(create_record(
        b"read1",
        cigar1,
        &[],
        std::iter::repeat_n(70_u8, len1)
            .collect::<Vec<u8>>()
            .as_slice(),
        md1,
        false,
    )?));

    let len2 = read_len_from_cigar(cigar2);
    let rec2 = Box::leak(Box::new(create_record(
        b"read1",
        cigar2,
        &[],
        std::iter::repeat_n(70_u8, len2)
            .collect::<Vec<u8>>()
            .as_slice(),
        md2,
        false,
    )?));

    let iter1 = UnifiedOpIterator::new(rec1)?;
    let iter2 = UnifiedOpIterator::new(rec2)?;

    Ok(PreparedAlignmentPair { iter1, iter2 })
}
#[test]
fn test_perfect_match() -> Result<()> {
    let mut pair = make_prepared_pair("10M", "10", "10M", "10").unwrap();
    let (is_host_perfect, is_graft_perfect) = pair.host_graft_are_perfect_match();
    assert!(is_host_perfect);
    assert!(is_graft_perfect);

    let mut pair = make_prepared_pair("10M", "10", "9M1I", "9").unwrap();
    let (is_host_perfect, is_graft_perfect) = pair.host_graft_are_perfect_match();
    assert!(is_host_perfect);
    assert!(!is_graft_perfect);

    let mut pair = make_prepared_pair("8M2D", "8^AC", "10M", "10").unwrap();
    let (is_host_perfect, is_graft_perfect) = pair.host_graft_are_perfect_match();
    assert!(!is_host_perfect);
    assert!(is_graft_perfect);

    let mut pair = make_prepared_pair("10M", "10", "10M", "9A1").unwrap();
    let (is_host_perfect, is_graft_perfect) = pair.host_graft_are_perfect_match();
    assert!(is_host_perfect);
    assert!(!is_graft_perfect);
    Ok(())
}
#[test]
fn test_prepared_alignment_pair_iter() -> Result<()> {
    let recs1 = vec![
        create_record(b"read1", "10M", &[], &[30; 10], "10", false)?,
        create_record(b"read2", "8M2I", &[], &[30; 10], "8", false)?,
        create_record(b"read3", "5M5D", &[], &[30; 5], "5", false)?,
    ];
    let recs2 = vec![
        create_record(b"read1", "10M", &[], &[30; 10], "10", false)?,
        create_record(b"read2", "10M", &[], &[30; 10], "10", false)?,
        create_record(b"read3", "5M5D", &[], &[30; 5], "5", false)?,
    ];

    let mut pair_iter = PreparedAlignmentPairIter::new(&recs1, &recs2);
    for _ in 0..3 {
        let pair = pair_iter.next().unwrap().unwrap();
        assert_eq!(pair.iter1.seq_len(), pair.iter2.seq_len());
    }
    assert!(pair_iter.next().is_none());
    Ok(())
}
#[test]
fn test_prepared_alignment_pair_iter_skips_secondary() -> Result<()> {
    let recs1 = vec![
        create_record(b"read1", "10M", &[], &[30; 10], "10", false)?,
        {
            let mut r = create_record(b"read1", "10M", &[], &[30; 10], "10", false)?;
            r.set_flags(r.flags() | 0x100); // secondary
            r
        },
        create_record(b"read2", "8M2I", &[], &[30; 10], "8", false)?,
    ];
    let recs2 = vec![
        create_record(b"read1", "10M", &[], &[30; 10], "10", false)?,
        create_record(b"read2", "10M", &[], &[30; 10], "10", false)?,
    ];

    let mut pair_iter = PreparedAlignmentPairIter::new(&recs1, &recs2);
    for _ in 0..2 {
        let pair = pair_iter.next().unwrap().unwrap();
        assert_eq!(pair.iter1.seq_len(), pair.iter2.seq_len());
    }
    assert!(pair_iter.next().is_none());
    Ok(())
}
#[test]
fn test_prepared_alignment_pair_iter_unequal_lengths() -> Result<()> {
    let recs1 = vec![
        create_record(b"read1", "10M", &[], &[30; 10], "10", false)?,
        create_record(b"read2", "8M2I", &[], &[30; 10], "8", false)?,
    ];
    let recs2 = vec![create_record(b"read1", "10M", &[], &[30; 10], "10", false)?];

    let mut pair_iter = PreparedAlignmentPairIter::new(&recs1, &recs2);
    let pair = pair_iter.next().unwrap().unwrap();
    assert_eq!(pair.iter1.seq_len(), pair.iter2.seq_len());
    assert!(pair_iter.next().is_none());
    Ok(())
}
#[test]
fn test_print_req() -> Result<()> {
    let rec = create_record(b"read1", "10M1I5M", &[], &[30; 16], "15", false)?;
    print_req(1, &rec);
    Ok(())
}
#[test]
fn test_read_len_from_cigar() -> Result<()> {
    assert_eq!(read_len_from_cigar("10M"), 10);
    assert_eq!(read_len_from_cigar("5M2I3M"), 10);
    assert_eq!(read_len_from_cigar("8M2D"), 8);
    assert_eq!(read_len_from_cigar("4S6M"), 10);
    assert_eq!(read_len_from_cigar("3H7M2S"), 9);
    Ok(())
}
#[test]
fn test_make_prepared_pair() -> Result<()> {
    let pair = make_prepared_pair("10M", "10", "10M", "10").unwrap();
    assert_eq!(pair.iter1.seq_len(), 10);
    assert_eq!(pair.iter2.seq_len(), 10);
    Ok(())
}
#[test]
fn test_make_prepared_pair_invalid() -> Result<()> {
    let result = make_prepared_pair("10M", "10", "9M1I", "9");
    assert!(result.is_ok());
    Ok(())
}
#[test]
fn test_make_prepared_pair_invalid_md() -> Result<()> {
    // htslib does not recognize invalid MD tags, so this passes
    let mut prepared = make_prepared_pair("10M", "5X4", "10M", "10")?;
    for op in prepared.iter2.by_ref() {
        assert!(op.is_ok());
    }
    for op in prepared.iter1.by_ref() {
        if op.is_err() {
            return Ok(());
        }
    }
    panic!("Expected error due to invalid MD tag");
}
#[test]
fn test_make_prepared_pair_invalid_cigar() -> Result<()> {
    let result = make_prepared_pair("10Z", "10", "10M", "10");
    assert!(result.is_err());
    Ok(())
}
#[test]
fn test_make_prepared_pair_unmapped() -> Result<()> {
    let rec1 = create_record(b"read1", "10M", &[], &[30; 10], "10", false)?;

    let iter1 = UnifiedOpIterator::new(&rec1).unwrap();
    let iter2 = UnifiedOpIterator::empty(false);

    let mut pair = PreparedAlignmentPair { iter1, iter2 };
    assert_eq!(pair.iter1.seq_len(), 10);
    match pair.iter2.next() {
        None => (),
        Some(_) => panic!("Expected no operations for unmapped alignment"),
    }
    Ok(())
}
#[test]
fn test_make_prepared_pair_unmapped_both() -> Result<()> {
    let iter1 = UnifiedOpIterator::empty(false);
    let iter2 = UnifiedOpIterator::empty(false);

    let mut pair = PreparedAlignmentPair { iter1, iter2 };
    match pair.iter1.next() {
        None => (),
        Some(_) => panic!("Expected no operations for unmapped alignment"),
    }
    match pair.iter2.next() {
        None => (),
        Some(_) => panic!("Expected no operations for unmapped alignment"),
    }
    Ok(())
}
#[test]
fn test_make_prepared_pair_different_lengths() -> Result<()> {
    let result = make_prepared_pair("10M", "10", "9M", "9");
    assert!(result.is_ok());
    Ok(())
}
#[test]
fn test_make_prepared_pair_different_lengths_invalid() -> Result<()> {
    let result = make_prepared_pair("10M", "10", "8M", "8");
    assert!(result.is_ok());
    Ok(())
}
#[test]
fn test_make_prepared_pair_different_lengths_invalid_md() -> Result<()> {
    // htslib does not recognize invalid MD tags, so this passes
    let mut prepared = make_prepared_pair("10M", "10", "9M", "8A1")?;
    for op in prepared.iter1.by_ref() {
        assert!(op.is_ok());
    }
    for op in prepared.iter2.by_ref() {
        if op.is_err() {
            return Ok(());
        }
    }
    panic!("Expected error due to invalid MD tag");
}
#[test]
fn test_make_prepared_pair_different_lengths_invalid_cigar() -> Result<()> {
    let result = make_prepared_pair("10M", "10", "9Z", "9");
    assert!(result.is_err());
    Ok(())
}
