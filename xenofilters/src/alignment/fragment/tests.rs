use crate::Config;
use crate::Penalty;
use crate::alignment::UnifiedOp;
use crate::alignment::fragment::Fragment;
use crate::alignment::fragment::calculate_translocation_penalty;
use crate::alignment::fragment::revcmp_encoded;
use crate::alignment::build_fragment;
use crate::tests::create_record;
use anyhow::Result;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::record::{Cigar, Record};
use smallvec::{SmallVec, smallvec};
use crate::variant::Variant;
use crate::alignment::stringify_record;

pub(crate) fn setup_penalties() -> Penalty {
    let c = Config::default();
    let mut p = c.to_penalties();
    p.log_likelihood_match = [0.0; 93]; // 0 log-lik for match
    p.log_likelihood_mismatch = [-1.0; 93];
    p.gap_open = -2.0;
    p.gap_extend = -0.5;
    p
}

#[test]
fn test_stitched_fragment_scoring() {
    let penalties = setup_penalties();
    let variants = SmallVec::<[SmallVec::<[&dyn Variant; 8]>; 2]>::new();
    let mut qual = SmallVec::<[&[u8]; 2]>::new();
    let mut seq = SmallVec::<[&[u8]; 2]>::new();
    let segment_revcmp = smallvec![false, false];
    for _ in 0..5 {
        qual.push(&[30; 5]);
        seq.push(&[b'A'; 5]);
    }
    let current_segment_idx = 0;
    let local_offset = 0;
    let mut fragment = Fragment {
        seq,
        qual,
        tid: 0,
        pos: 100,
        ops: Box::new(
            vec![
                Ok(UnifiedOp::Match(2)),
                Ok(UnifiedOp::Mis(1)),
                Ok(UnifiedOp::Ins(1)), // Gap open
                Ok(UnifiedOp::Del(2)), // Gap open + 1 extend
            ]
            .into_iter(),
        ),
        current_segment_idx,
        local_offset,
        segment_revcmp,
        penalties: &penalties,
        variants,
    };

    let score = fragment.score().unwrap();
    // Match(2)*0 + Mis(1)*-1 + Ins(1)*-2 + Del(2)*(-2 + -0.5) = -5.5
    assert_eq!(score, -5.5);
}

#[test]
fn test_translocation_penalty_logic() {
    let penalties = setup_penalties();
    // Create a record with 10M (perfect) and 5S (soft clip)
    // High quality on soft clip = high penalty
    // High quality on matches = offsets penalty (lowers it)
    let mut record = Record::new();
    let vec_cig = vec![Cigar::Match(10), Cigar::SoftClip(5)];
    let cigar = CigarString(vec_cig);
    record.set(b"read1", Some(&cigar), &[b'A'; 15], &[30; 15]);

    let penalty = calculate_translocation_penalty(&penalties, &record).unwrap();

    // 5 * |-1.0| (mismatch for softclips) - 10 * 0.0 (match for matches) = 5.0
    assert!(penalty > 0.0);
    assert_eq!(penalty, 5.0);
}

#[test]
fn test_stitched_fragment_creation() -> Result<()> {
    let penalties = setup_penalties();

    let record1 = create_record(b"read1", "5M3S", &[b'A'; 8], &[30; 8], "5", false)?;
    let record2 = create_record(b"read1", "4M4S", &[b'A'; 8], &[30; 8], "4", false)?;

    let records = vec![record1, record2];
    let order = smallvec![0, 1];
    let variants = SmallVec::<[SmallVec::<[&dyn Variant; 8]>; 2]>::new();

    let stitched = build_fragment(&penalties, &records, order, variants).unwrap();

    // Check that the stitched fragment has the correct TID and POS
    assert_eq!(stitched.tid, records[0].tid());
    assert_eq!(stitched.pos, records[0].pos());
    Ok(())
}

#[test]
fn test_revcmp_encoded() {
    // A <-> T
    assert_eq!(revcmp_encoded(1), 8);
    assert_eq!(revcmp_encoded(8), 1);

    // C <-> G
    assert_eq!(revcmp_encoded(2), 4);
    assert_eq!(revcmp_encoded(4), 2);

    // N stays N
    assert_eq!(revcmp_encoded(15), 15);

    // '=' or garbage stays unchanged
    assert_eq!(revcmp_encoded(0), 0);
    assert_eq!(revcmp_encoded(42), 42);
}

#[test]
fn test_calculate_translocation_penalty() -> Result<()> {
    let penalties = setup_penalties();
    let record = create_record(b"read1", "6M4S", &[b'A'; 10], &[30; 10], "6", false)?;

    let penalty = calculate_translocation_penalty(&penalties, &record).unwrap();
    // 4 * |-1.0| - 6 * 0.0 = 4.0
    assert_eq!(penalty, 4.0);
    Ok(())
}

#[test]
fn test_calculate_translocation_penalty_no_softclip() -> Result<()> {
    let penalties = setup_penalties();
    let record = create_record(b"read1", "10M", &[b'A'; 10], &[30; 10], "10", false)?;

    let penalty = calculate_translocation_penalty(&penalties, &record).unwrap();
    // No soft clips, so penalty should be 0.0
    assert_eq!(penalty, 0.0);
    Ok(())
}

#[test]
fn test_calculate_translocation_penalty_high_quality_softclip() -> Result<()> {
    let penalties = setup_penalties();
    let record = create_record(
        b"read1",
        "5M5S",
        &[b'A'; 10],
        &[40; 10], // High quality scores
        "5",
        false,
    )?;

    let penalty = calculate_translocation_penalty(&penalties, &record).unwrap();
    // 5 * |-1.0| - 5 * 0.0 = 5.0
    assert_eq!(penalty, 5.0);
    Ok(())
}
