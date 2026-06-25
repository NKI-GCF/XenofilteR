use super::*;
use crate::tests::create_record;

#[test]
fn test_is_perfect_single_op_no_mismatches() -> Result<()> {
    let rec = create_record(b"r", "10M", &[], &[], "10", false)?;
    let flags = rec.flags();
    assert!(MdCigFlags::try_from_record(&rec, &flags)?.is_perfect());
    Ok(())
}

#[test]
fn test_is_perfect_false_with_mismatch() -> Result<()> {
    let rec = create_record(b"r", "10M", &[], &[], "5A4", false)?;
    let flags = rec.flags();
    assert!(!MdCigFlags::try_from_record(&rec, &flags)?.is_perfect());
    Ok(())
}

#[test]
fn test_is_perfect_false_multiple_cigar_ops_even_if_md_clean() -> Result<()> {
    let rec = create_record(b"r", "5M5S", &[], &[], "5", false)?;
    let flags = rec.flags();
    assert!(!MdCigFlags::try_from_record(&rec, &flags)?.is_perfect());
    Ok(())
}

#[test]
fn test_is_reverse_complemented() -> Result<()> {
    let fwd = create_record(b"r", "5M", &[], &[], "5", false)?;
    let rev = create_record(b"r", "5M", &[], &[], "5", true)?;
    assert!(!MdCigFlags::try_from_record(&fwd, &fwd.flags())?.is_reverse_complemented());
    assert!(MdCigFlags::try_from_record(&rev, &rev.flags())?.is_reverse_complemented());
    Ok(())
}

#[test]
fn test_is_last_segment() -> Result<()> {
    let mut last = create_record(b"r", "5M", &[], &[], "5", false)?;
    *last.flags_mut() = Flags::from_bits(0x80).unwrap(); // last segment, mapped
    let mut first = create_record(b"r", "5M", &[], &[], "5", false)?;
    *first.flags_mut() = Flags::from_bits(0x40).unwrap(); // first segment, mapped

    assert!(MdCigFlags::try_from_record(&last, &last.flags())?.is_last_segment());
    assert!(!MdCigFlags::try_from_record(&first, &first.flags())?.is_last_segment());
    Ok(())
}

#[test]
fn test_try_from_record_unmapped_is_rejected() {
    let rec = create_record(b"r", "", &[], &[], "", false).unwrap(); // unmapped
    let flags = rec.flags();
    assert!(MdCigFlags::try_from_record(&rec, &flags).is_err());
}

#[test]
fn test_try_from_record_missing_md_tag_errors() -> Result<()> {
    let mut rec = create_record(b"r", "5M", &[], &[], "5", false)?;
    *rec.data_mut() = Default::default(); // strip the MD tag that create_record set
    let flags = rec.flags();
    assert!(MdCigFlags::try_from_record(&rec, &flags).is_err());
    Ok(())
}

#[test]
fn test_partial_ord_and_eq() -> Result<()> {
    // Two perfect records
    let p1 = create_record(b"r1", "10M", &[], &[], "10", false)?;
    let p2 = create_record(b"r2", "10M", &[], &[], "10", false)?;

    // Two imperfect records (containing mismatches in MD tag)
    let i1 = create_record(b"r3", "10M", &[], &[], "5A4", false)?;
    let i2 = create_record(b"r4", "10M", &[], &[], "5A4", false)?;

    let p1_flags = p1.flags();
    let md_p1 = MdCigFlags::try_from_record(&p1, &p1_flags)?;
    let p2_flags = p2.flags();
    let md_p2 = MdCigFlags::try_from_record(&p2, &p2_flags)?;
    let i1_flags = i1.flags();
    let md_i1 = MdCigFlags::try_from_record(&i1, &i1_flags)?;
    let i2_flags = i2.flags();
    let md_i2 = MdCigFlags::try_from_record(&i2, &i2_flags)?;

    // --- Testing PartialEq ---

    assert!(md_p1 == md_p2, "perfect == perfect should be true");

    assert!(md_p1 != md_i1, "perfect != imperfect should be true");

    assert!(
        !(md_i1 == md_i2),
        "imperfect == imperfect defaults to None, so eq is false"
    );

    // --- Testing PartialOrd ---

    assert_eq!(md_p1.partial_cmp(&md_p2), Some(Ordering::Equal));
    assert_eq!(md_p1.partial_cmp(&md_i1), Some(Ordering::Less));
    assert_eq!(md_i1.partial_cmp(&md_p1), Some(Ordering::Greater));
    assert_eq!(md_i1.partial_cmp(&md_i2), None);

    Ok(())
}
