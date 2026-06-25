use super::*;
use crate::alignment::fragment::SimpleRec;
use crate::alignment::MdCigFlags;
use crate::config::Config;
use crate::penalty::Penalty;
use crate::tests::create_record;
use crate::variant::{Eval, Variant};
use anyhow::Result;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{smallvec, SmallVec};

pub(crate) fn setup_penalties() -> Penalty {
    let c = Config::default();
    let mut p = c.to_penalties();
    p.log_likelihood_match = [0.0; 93]; // 0 log-like for match
    p.log_likelihood_mismatch = [-1.0; 93];
    p.gap_open = -2.0;
    p.gap_extend = -0.5;
    p
}

struct TestVariant {
    pos: usize,
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
    p_variant: f64,
}
impl TestVariant {
    fn new(pos: usize, len: usize, p_variant: f64) -> Self {
        TestVariant {
            pos,
            ref_a: vec![b'A'; len],
            alt_a: vec![b'G'; len],
            p_variant,
        }
    }
}
impl Variant for TestVariant {
    fn pos(&self) -> usize {
        self.pos
    }
    fn ref_allele(&self) -> &[u8] {
        &self.ref_a
    }
    fn alt_allele(&self) -> &[u8] {
        &self.alt_a
    }
    fn p_variant(&self) -> f64 {
        self.p_variant
    }
}

fn make_eval<'a>(v: &'a dyn Variant) -> Eval<'a> {
    let mut e = Eval::new();
    e.set_variant(v);
    e
}
/* conversion to Record not possible
#[test]
fn test_simple_rec_for_bam_record() -> Result<()> {
    use noodles::bam;
    use noodles::sam::Header;

    let header = Header::default();
    let mut buf = RecordBuf::default();
    *buf.reference_sequence_id_mut() = Some(42);
    *buf.quality_scores_mut() =
        noodles::sam::alignment::record_buf::QualityScores::from(vec![30, 31, 32]);

    // Convert to a BAM record
    let bam_record = bam::Record::try_from_alignment_record(&header, &buf)?;

    assert_eq!(bam_record.ref_seq_id().unwrap()?, 42);

    assert_eq!(bam_record.quality_at(0), Some(30));
    assert_eq!(bam_record.quality_at(2), Some(32));
    assert_eq!(bam_record.quality_at(3), None);

    let converted = bam_record.as_record_buf(&header)?;
    assert_eq!(converted.reference_sequence_id(), Some(42));

    Ok(())
}*/

#[test]
fn test_fragment_requires_revcmp() -> Result<()> {
    let rec_fwd = create_record(b"read1", "5M", b"AAAAA", &[30; 5], "5", false)?;
    let mut rec_rev = create_record(b"read1", "5M", b"AAAAA", &[30; 5], "10", false)?;
    *rec_rev.flags_mut() = Flags::from_bits(0x10).unwrap(); // Reverse complemented

    let flags_fwd = rec_fwd.flags();
    let flags_rev = rec_rev.flags();

    let p = setup_penalties();
    let mut fragment = Fragment::new(
        &p,
        smallvec![&rec_fwd, &rec_rev],
        smallvec![
            MdCigFlags::try_from_record(&rec_fwd, &flags_fwd)?,
            MdCigFlags::try_from_record(&rec_rev, &flags_rev)?
        ],
    )?;

    // fragment is initialized with seg_i = 0 (the forward segment)

    assert!(!fragment.requires_revcmp(0)); // Same segment -> false

    assert!(fragment.requires_revcmp(1)); // Other segment has different ori -> true

    // Move to the revcmp segment
    fragment.seg_i = 1;
    assert!(fragment.requires_revcmp(0)); // Other segment has different ori -> true
    assert!(!fragment.requires_revcmp(1)); // Same segment -> false

    Ok(())
}

#[test]
fn test_score_variants_in_window_boundaries() -> Result<()> {
    // Start the read at position 1. Spans reference [1, 11).
    let rec = create_record(b"read1", "10M", &[b'A'; 10], &[30; 10], "1", false)?;
    let flags = rec.flags();
    let p = setup_penalties();
    let fragment = Fragment::new(
        &p,
        smallvec![&rec],
        smallvec![MdCigFlags::try_from_record(&rec, &flags)?],
    )?;

    // We create three variants around the right edge of a window ending at 6:
    // v1: fully inside the window (start=3, ref_end=4)
    // v2: overlapping the right edge (start=5, ref_end=7)
    // v3: starts exactly on the right edge (start=6, ref_end=7) -> should NOT be processed
    let v1 = TestVariant::new(3, 1, 0.1);
    let v2 = TestVariant::new(5, 2, 0.1);
    let v3 = TestVariant::new(6, 1, 0.1);

    let mut dvnt: FragEvalVec =
        smallvec![smallvec![make_eval(&v1), make_eval(&v2), make_eval(&v3)]];
    let mut finished: FragEvalVec = smallvec![smallvec![]];
    let mut scratch = Scratch::new();

    // Create a WindowCtx that covers [1, 6)
    let ctx = WindowCtx::new(0, 0, 1, 6, 0.0);

    fragment.score_variants_in_window(&mut scratch, &mut dvnt, &mut finished, ctx)?;

    // v1 is fully processed (ref_end 4 <= 6). It gets moved to `finished`.
    assert_eq!(finished[0].len(), 1);
    assert_eq!(finished[0][0].start(), 3);

    // dvnt should now hold v2 and v3.
    // v2 was partially processed (start 5 < 6, but ref_end 7 > 6), so it remains in `dvnt`.
    // v3 was NOT processed because its start (6) is NOT < ctx.ref_end (6).
    assert_eq!(dvnt[0].len(), 2);
    assert_eq!(dvnt[0][0].start(), 5);
    assert_eq!(dvnt[0][1].start(), 6);

    Ok(())
}

#[test]
fn maximize_delta_exact_zero_mutant() {
    let mut dvnt: FragEvalVec = smallvec![smallvec![mk_eval(10, 1, 1, 0.0)]];
    let mut dp = SmallVec::new();

    // If it were >=, the variant with delta 0.0 would be included and processed.
    assert_eq!(maximize_delta(&mut dvnt, &mut dp), 0.0);
    assert!(
        dp.is_empty(),
        "dp should be empty because variants with <= 0 delta should be filtered out"
    );
}

#[test]
fn test_simple_rec_for_sam_record_buf() {
    let mut buf = RecordBuf::default();
    *buf.reference_sequence_id_mut() = Some(15);

    assert_eq!(buf.ref_seq_id().unwrap().unwrap(), 15);
}

#[test]
fn test_stitched_fragment_creation() -> Result<()> {
    let record1 = create_record(b"read1", "5M3S", &[b'A'; 8], &[30; 8], "5", false)?;
    let record2 = create_record(b"read1", "4M4S", &[b'A'; 8], &[30; 8], "4", false)?;
    let flags1 = record1.flags();
    let flags2 = record2.flags();

    let records: SmallVec<[&RecordBuf; READ_CT]> = smallvec![&record1, &record2];
    let p = setup_penalties();
    let mut md_cig_flags: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();
    md_cig_flags.push(MdCigFlags::try_from_record(&record1, &flags1)?);
    md_cig_flags.push(MdCigFlags::try_from_record(&record2, &flags2)?);
    let _stitched = Fragment::new(&p, records, md_cig_flags)?;
    Ok(())
}
#[test]
fn test_complement() {
    assert_eq!(complement(b'A'), b'T');
    assert_eq!(complement(b'C'), b'G');
    assert_eq!(complement(b'G'), b'C');
    assert_eq!(complement(b'T'), b'A');
    assert_eq!(complement(b'N'), b'N');
}

#[test]
fn test_q() -> Result<()> {
    let record = create_record(
        b"read1",
        "5M",
        &[b'A'; 5],
        &[30, 31, 32, 33, 34],
        "5",
        false,
    )?;
    let flags = record.flags();
    let seg = smallvec![&record];
    let md_cig_flags = smallvec![MdCigFlags::try_from_record(&record, &flags)?];
    let p = setup_penalties();
    let fragment = Fragment::new(&p, seg, md_cig_flags)?;
    assert_eq!(fragment.q(0, 0)?, 30);
    assert_eq!(fragment.q(0, 4)?, 34);
    assert!(fragment.q(0, 5).is_err());
    Ok(())
}
fn mk_eval(pos: usize, ref_len: usize, alt_len: usize, delta: f64) -> Eval<'static> {
    let v: &'static TestVariant = Box::leak(Box::new(TestVariant {
        pos,
        ref_a: vec![b'A'; ref_len],
        alt_a: vec![b'A'; alt_len],
        p_variant: 0.1,
    }));
    let mut eval = Eval::new();
    eval.set_variant(v);
    eval.update(0.0, delta); // delta() = alt_score - incurred = delta - 0 = delta
    eval
}

#[test]
fn maximize_delta_empty_returns_zero() {
    let mut dvnt: FragEvalVec = smallvec![smallvec![]];
    let mut dp = SmallVec::new();

    assert_eq!(maximize_delta(&mut dvnt, &mut dp), 0.0);
}

#[test]
fn maximize_delta_ignores_negative_deltas() {
    let mut dvnt: FragEvalVec =
        smallvec![smallvec![mk_eval(10, 1, 1, -5.0), mk_eval(20, 1, 1, -10.0),]];

    let mut dp = SmallVec::new();

    assert_eq!(maximize_delta(&mut dvnt, &mut dp), 0.0);
}

#[test]
fn maximize_delta_single_variant() {
    let mut dvnt: FragEvalVec = smallvec![smallvec![mk_eval(10, 1, 1, 7.5)]];

    let mut dp = SmallVec::new();

    assert!((maximize_delta(&mut dvnt, &mut dp) - 7.5).abs() < 1e-9);
}

#[test]
fn maximize_delta_sums_non_overlapping_variants() {
    let mut dvnt: FragEvalVec = smallvec![smallvec![
        mk_eval(10, 1, 1, 5.0),
        mk_eval(20, 1, 1, 7.0),
        mk_eval(30, 1, 1, 11.0),
    ]];

    let mut dp = SmallVec::new();

    assert!((maximize_delta(&mut dvnt, &mut dp) - 23.0).abs() < 1e-9);
}

#[test]
fn maximize_delta_prefers_larger_overlapping_variant() {
    let mut dvnt: FragEvalVec = smallvec![smallvec![
        mk_eval(10, 10, 10, 5.0), // [10,20)
        mk_eval(15, 10, 10, 8.0), // [15,25)
    ]];

    let mut dp = SmallVec::new();

    assert!((maximize_delta(&mut dvnt, &mut dp) - 8.0).abs() < 1e-9);
}

#[test]
fn maximize_delta_finds_best_chain() {
    let mut dvnt: FragEvalVec = smallvec![smallvec![
        mk_eval(10, 10, 10, 5.0),   // [10,20)
        mk_eval(15, 10, 10, 100.0), // [15,25)
        mk_eval(30, 10, 10, 6.0),   // [30,40)
    ]];

    let mut dp = SmallVec::new();

    assert!((maximize_delta(&mut dvnt, &mut dp) - 106.0).abs() < 1e-9);
}

#[test]
fn maximize_delta_handles_nested_variants() {
    let mut dvnt: FragEvalVec = smallvec![smallvec![
        mk_eval(10, 20, 20, 5.0), // [10,30)
        mk_eval(15, 2, 2, 20.0),  // [15,17)
        mk_eval(40, 1, 1, 3.0),   // [40,41)
    ]];

    let mut dp = SmallVec::new();

    assert!((maximize_delta(&mut dvnt, &mut dp) - 23.0).abs() < 1e-9);
}

#[test]
fn maximize_delta_boundary_touching_is_non_overlapping() {
    let mut dvnt: FragEvalVec = smallvec![smallvec![
        mk_eval(10, 10, 10, 5.0), // [10,20)
        mk_eval(20, 10, 10, 7.0), // [20,30)
    ]];

    let mut dp = SmallVec::new();

    assert!((maximize_delta(&mut dvnt, &mut dp) - 12.0).abs() < 1e-9);
}

#[test]
fn maximize_delta_uses_max_ref_or_alt_end_for_insertions() {
    let mut dvnt: FragEvalVec = smallvec![smallvec![
        mk_eval(10, 1, 20, 10.0), // end = 30
        mk_eval(25, 1, 1, 50.0),  // overlaps insertion span
    ]];

    let mut dp = SmallVec::new();

    assert!((maximize_delta(&mut dvnt, &mut dp) - 50.0).abs() < 1e-9);
}

#[test]
fn snp_alt_support_gives_positive_delta() -> Result<()> {
    let rec = create_record(b"read1", "5M", b"AAGAA", &[30; 5], "5", false)?;
    let flags = rec.flags();
    let p = setup_penalties();
    let mut frag = Fragment::new(
        &p,
        smallvec![&rec],
        smallvec![MdCigFlags::try_from_record(&rec, &flags)?],
    )?;

    // p_variant must exceed 0.5 for the current formula to ever favor alt
    // over the weighted-reference baseline (delta = (2p-1)*(lm-lmm)).
    let v = TestVariant::new(3, 1, 1.0);

    let mut dvnt = smallvec![smallvec![make_eval(&v)]];
    let mut scratch = Scratch::new();
    let score = frag.score(&mut scratch, &mut dvnt)?;

    assert!(
        score > 0.0,
        "Expected positive score for read supporting alt allele, got {score}"
    );
    Ok(())
}

#[test]
fn snp_ref_support_gives_no_bonus() -> Result<()> {
    let rec = create_record(b"read1", "5M", b"AAAAA", &[30; 5], "5", false)?;

    let flags = rec.flags();

    let p = setup_penalties();

    let mut frag = Fragment::new(
        &p,
        smallvec![&rec],
        smallvec![MdCigFlags::try_from_record(&rec, &flags)?],
    )?;

    let v = TestVariant::new(3, 1, 0.1);

    let mut dvnt = smallvec![smallvec![make_eval(&v)]];
    let mut scratch = Scratch::new();

    let score = frag.score(&mut scratch, &mut dvnt)?;

    assert!(score <= 0.0);

    Ok(())
}

#[test]
fn test_maximize_delta_no_variants_is_zero() -> Result<()> {
    let mut dvnt: FragEvalVec<'_> = smallvec![SmallVec::new()];
    let mut dp = SmallVec::new();
    assert_eq!(maximize_delta(&mut dvnt, &mut dp), 0.0);
    Ok(())
}

#[test]
fn test_maximize_delta_ignores_non_positive_delta() -> Result<()> {
    let v = TestVariant::new(0, 1, 0.1);
    let mut e = Eval::new();
    e.set_variant(&v);
    e.update(5.0, 5.0); // delta == 0, discarded out
    let mut dvnt: FragEvalVec<'_> = smallvec![smallvec![e]];
    let mut dp = SmallVec::new();
    assert_eq!(maximize_delta(&mut dvnt, &mut dp), 0.0);
    Ok(())
}

#[test]
fn test_maximize_delta_sums_non_overlapping_variants() -> Result<()> {
    let v1 = TestVariant::new(0, 1, 0.1); // end = 1
    let v2 = TestVariant::new(10, 1, 0.1); // end = 11
    let mut e1 = Eval::new();
    e1.set_variant(&v1);
    e1.update(0.0, 3.0); // delta 3
    let mut e2 = Eval::new();
    e2.set_variant(&v2);
    e2.update(0.0, 2.0); // delta 2
    let mut dvnt: FragEvalVec<'_> = smallvec![smallvec![e1, e2]];
    let mut dp = SmallVec::new();
    assert_eq!(maximize_delta(&mut dvnt, &mut dp), 5.0);
    Ok(())
}

#[test]
fn test_maximize_delta_picks_best_not_sum_for_overlapping_variants() -> Result<()> {
    let v1 = TestVariant::new(0, 5, 0.1); // [0,5)
    let v2 = TestVariant::new(2, 5, 0.1); // [2,7) overlaps v1
    let mut e1 = Eval::new();
    e1.set_variant(&v1);
    e1.update(0.0, 10.0); // delta 10
    let mut e2 = Eval::new();
    e2.set_variant(&v2);
    e2.update(0.0, 3.0); // delta 3
    let mut dvnt: FragEvalVec<'_> = smallvec![smallvec![e1, e2]];
    let mut dp = SmallVec::new();
    assert_eq!(maximize_delta(&mut dvnt, &mut dp), 10.0); // not 13.0
    Ok(())
}

#[test]
fn snp_no_alt_support_gives_nonpositive_delta() -> Result<()> {
    // Read has 'A' at the variant position (matches declared ref, not alt) — no rescue expected.
    let rec = create_record(b"read1", "5M", b"AAAAA", &[30; 5], "5", false)?;
    let flags = rec.flags();
    let p = setup_penalties();
    let mut frag = Fragment::new(
        &p,
        smallvec![&rec],
        smallvec![MdCigFlags::try_from_record(&rec, &flags)?],
    )?;
    let v = TestVariant::new(3, 1, 0.1);
    let mut dvnt = smallvec![smallvec![make_eval(&v)]];
    let mut scratch = Scratch::new();
    let score = frag.score(&mut scratch, &mut dvnt)?;
    assert!(score <= 0.0);
    Ok(())
}

#[test]
fn test_test_variant_p_variant_reflects_constructed_field() {
    let v = TestVariant {
        pos: 0,
        ref_a: b"A".to_vec(),
        alt_a: b"G".to_vec(),
        p_variant: 1.0,
    };
    assert_eq!(v.p_variant(), 1.0);
    let v2 = TestVariant {
        pos: 0,
        ref_a: b"A".to_vec(),
        alt_a: b"G".to_vec(),
        p_variant: 0.25,
    };
    assert_eq!(v2.p_variant(), 0.25);
}

#[test]
fn snp_alt_support_with_low_prior_is_not_rescued() -> Result<()> {
    let rec = create_record(b"read1", "5M", b"AAGAA", &[30; 5], "5", false)?;
    let flags = rec.flags();
    let p = setup_penalties();
    let mut frag = Fragment::new(
        &p,
        smallvec![&rec],
        smallvec![MdCigFlags::try_from_record(&rec, &flags)?],
    )?;
    let v = TestVariant::new(3, 1, 0.1);
    let mut dvnt = smallvec![smallvec![make_eval(&v)]];
    let mut scratch = Scratch::new();
    let score = frag.score(&mut scratch, &mut dvnt)?;
    assert!(
        score <= 0.0,
        "expected no rescue at p_variant <= 0.5, got {score}"
    );
    Ok(())
}

#[test]
fn snp_alt_support_boundary_p_variant_half_gives_zero_delta() -> Result<()> {
    let rec = create_record(b"read1", "5M", b"AAGAA", &[30; 5], "5", false)?;
    let flags = rec.flags();
    let p = setup_penalties();
    let mut frag = Fragment::new(
        &p,
        smallvec![&rec],
        smallvec![MdCigFlags::try_from_record(&rec, &flags)?],
    )?;
    let v = TestVariant::new(3, 1, 0.5);
    let mut dvnt = smallvec![smallvec![make_eval(&v)]];
    let mut scratch = Scratch::new();
    let score = frag.score(&mut scratch, &mut dvnt)?;
    assert!(score.abs() < 1e-9);
    Ok(())
}
