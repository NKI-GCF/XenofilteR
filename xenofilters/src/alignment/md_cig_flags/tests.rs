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
