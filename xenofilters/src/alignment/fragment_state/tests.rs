use crate::alignment::fragment_state::FragmentState;
use crate::tests::create_record;
use anyhow::Result;
use std::cmp::Ordering;
use smallvec::{SmallVec, smallvec};
use noodles::sam::alignment::record::Flags;

#[test]
fn test_fragment_state_ordering() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let state1 = FragmentState::from_record(rec1, 0);
    let state2 = FragmentState::from_record(rec2, 1);
    assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Equal));
    Ok(())
}

#[test]
fn test_fragment_state_first_qname() -> Result<()> {
    let qual = vec![37; 100];
    let rec = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let state = FragmentState::from_record(rec, 0);
    assert_eq!(state.first_qname(), b"read1");
    Ok(())
}

#[test]
fn test_fragment_state_get_nr() -> Result<()> {
    let qual = vec![37; 100];
    let rec = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let state = FragmentState::from_record(rec, 42);
    assert_eq!(state.get_nr(), 42);
    Ok(())
}

#[test]
fn test_fragment_state_equality() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let state1 = FragmentState::from_record(rec1, 0);
    let state2 = FragmentState::from_record(rec2, 0);
    assert_eq!(state1, state2);
    Ok(())
}

#[test]
fn test_fragment_state_inequality() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read2", "100M", &[], &qual, "100", false)?;
    let state1 = FragmentState::from_record(rec1, 0);
    let state2 = FragmentState::from_record(rec2, 0);
    assert_ne!(state1, state2);
    Ok(())
}

#[test]
fn test_fragment_state_multiple_records() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let mut state = FragmentState::from_record(rec1, 0);
    state.records.push(rec2);
    assert_eq!(state.records.len(), 2);
    assert_eq!(state.first_qname(), b"read1");
    Ok(())
}

#[test]
fn test_fragment_state_order_mates_multiple_records() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "100", true)?;
    let mut state = FragmentState::from_record(rec1, 0);
    state.records.push(rec2);
    let order = state.order_mates();
    let expected : SmallVec<[usize; 2]> = smallvec![0, 1];
    assert_eq!(order, expected); // Forward read should come before reverse read
    Ok(())
}

#[test]
fn test_fragment_state_partial_ord_with_unmapped() -> Result<()> {
    let qual = vec![37; 100];
    let seq = vec![b'A'; 100];
    let mut rec1 = create_record(b"read1", "", &seq, &qual, "", false)?;
    let rec2 = create_record(b"read2", "100M", &seq, &qual, "", false)?;
    rec1.flags_mut().toggle(Flags::from_bits(0x4).unwrap()); // Set unmapped flag
    let state1 = FragmentState::from_record(rec1, 0);
    let state2 = FragmentState::from_record(rec2, 0);
    assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Less));
    Ok(())
}

#[test]
fn test_fragment_state_partial_ord_both_unmapped() -> Result<()> {
    let qual = vec![37; 100];
    let seq = vec![b'A'; 100];
    let mut rec1 = create_record(b"read1", "", &seq, &qual, "", false)?;
    let mut rec2 = create_record(b"read2", "", &seq, &qual, "", false)?;
    rec1.flags_mut().toggle(Flags::from_bits(0x4).unwrap()); // Set unmapped flag
    rec2.flags_mut().toggle(Flags::from_bits(0x4).unwrap()); // Set unmapped flag
    let state1 = FragmentState::from_record(rec1, 0);
    let state2 = FragmentState::from_record(rec2, 0);
    assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Equal));
    Ok(())
}

#[test]
fn test_fragment_state_partial_ord_no_quick_balance() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "80T20", false)?;
    let state1 = FragmentState::from_record(rec1, 0);
    let state2 = FragmentState::from_record(rec2, 0);
    assert_eq!(state1.partial_cmp(&state2), None); // No quick balance
    Ok(())
}

#[test]
fn test_fragment_state_partial_ord_perfect_vs_imperfect() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let state1 = FragmentState::from_record(rec1, 0);
    let state2 = FragmentState::from_record(rec2, 0);
    assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Greater)); // Perfect match is better
    Ok(())
}

#[test]
fn test_fragment_state_partial_ord_imperfect_vs_perfect() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let state1 = FragmentState::from_record(rec1, 0);
    let state2 = FragmentState::from_record(rec2, 0);
    assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Less)); // Perfect match is better
    Ok(())
}

#[test]
fn test_fragment_state_partial_ord_equal_imperfects() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let state1 = FragmentState::from_record(rec1, 0);
    let state2 = FragmentState::from_record(rec2, 0);
    assert_eq!(state1.partial_cmp(&state2), None); // Same imperfect matches
    Ok(())
}

#[test]
fn test_fragment_state_partial_ord_multiple_records() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let mut state1 = FragmentState::from_record(rec1, 0);
    state1
        .records
        .push(create_record(b"read1", "100M", &[], &qual, "100", false)?);
    let mut state2 = FragmentState::from_record(rec2, 0);
    state2
        .records
        .push(create_record(b"read1", "100M", &[], &qual, "90A10", false)?);
    assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Greater)); // Perfect matches are better
    Ok(())
}

#[test]
fn test_fragment_state_partial_ord_multiple_records_no_quick_balance() -> Result<()> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "80T20", false)?;
    let mut state1 = FragmentState::from_record(rec1, 0);
    state1
        .records
        .push(create_record(b"read1", "100M", &[], &qual, "85G15", false)?);
    let mut state2 = FragmentState::from_record(rec2, 0);
    state2
        .records
        .push(create_record(b"read1", "100M", &[], &qual, "80T20", false)?);
    assert_eq!(state1.partial_cmp(&state2), None); // No quick balance
    Ok(())
}
