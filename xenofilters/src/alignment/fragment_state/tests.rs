use crate::alignment::fragment_state::FragmentState;
use crate::tests::create_record;
use noodles::core::Position;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{smallvec, SmallVec};
use std::cmp::Ordering;
use crate::Error;

fn segment(qname: &[u8], flag_bits: u16, tid: usize, start: usize) -> Result<RecordBuf, Error> {
    let mut rec = create_record(qname, "5M", &[], &[], "5", false)?;
    *rec.flags_mut() = Flags::from_bits(flag_bits).unwrap();
    *rec.reference_sequence_id_mut() = Some(tid);
    *rec.alignment_start_mut() = Some(Position::new(start).expect("nonzero position"));
    Ok(rec)
}

// Tests ok
#[test]
fn test_fragment_state_multiple_records() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let mut state = FragmentState::from_record(rec1, 0)?;
    state.add_record(rec2)?;
    assert_eq!(state.records.len(), 2);
    assert_eq!(state.first_qname(), b"read1");
    Ok(())
}
#[test]
fn test_fragment_state_first_qname() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let state = FragmentState::from_record(rec, 0)?;
    assert_eq!(state.first_qname(), b"read1");
    Ok(())
}
#[test]
fn test_fragment_state_equality() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let state1 = FragmentState::from_record(rec1, 0)?;
    let state2 = FragmentState::from_record(rec2, 0)?;
    assert_eq!(state1, state2);
    Ok(())
}
#[test]
fn test_fragment_state_partial_ord_equal_imperfects() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let state1 = FragmentState::from_record(rec1, 0)?;
    let state2 = FragmentState::from_record(rec2, 0)?;
    assert_eq!(state1.partial_cmp(&state2), None); // Same imperfect matches
    Ok(())
}
#[test]
fn test_fragment_state_get_nr() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let state = FragmentState::from_record(rec, 42)?;
    assert_eq!(state.get_nr(), 42);
    Ok(())
}
#[test]
fn test_fragment_state_inequality() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read2", "100M", &[], &qual, "100", false)?;
    let state1 = FragmentState::from_record(rec1, 0)?;
    let state2 = FragmentState::from_record(rec2, 0)?;
    assert_ne!(state1, state2);
    Ok(())
}
#[test]
fn test_fragment_state_partial_ord_multiple_records_no_quick_balance() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "80T20", false)?;
    let mut state1 = FragmentState::from_record(rec1, 0)?;
    state1.add_record(create_record(b"read1", "100M", &[], &qual, "85G15", false)?)?;
    let mut state2 = FragmentState::from_record(rec2, 0)?;
    state2.add_record(create_record(b"read1", "100M", &[], &qual, "80T20", false)?)?;
    assert_eq!(state1.partial_cmp(&state2), None); // No quick balance
    Ok(())
}
#[test]
fn test_fragment_state_partial_ord_no_quick_balance() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "80T20", false)?;
    let state1 = FragmentState::from_record(rec1, 0)?;
    let state2 = FragmentState::from_record(rec2, 0)?;
    assert_eq!(state1.partial_cmp(&state2), None); // No quick balance
    Ok(())
}

//tests fail
#[test]
fn test_fragment_state_order_mates_multiple_records() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "100", true)?;
    let mut state = FragmentState::from_record(rec1, 0)?;
    state.add_record(rec2)?;
    let order = state.order_mates();
    let expected: SmallVec<[usize; 2]> = smallvec![0, 1];
    assert_eq!(order, expected); // Forward read should come before reverse read
    Ok(())
}
#[test]
fn test_fragment_state_partial_ord_multiple_records() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let mut state1 = FragmentState::from_record(rec1, 0)?;
    state1.add_record(create_record(b"read1", "100M", &[], &qual, "100", false)?)?;
    let mut state2 = FragmentState::from_record(rec2, 0)?;
    state2.add_record(create_record(b"read1", "100M", &[], &qual, "90A10", false)?)?;
    assert_eq!(state1.partial_cmp(&state2), None); // Perfect matches are better
    let mut ord: Option<Ordering> = None;
    let _ = state1.cmp_perfect(&state2, &mut ord)?;
    assert_eq!(ord, Some(Ordering::Greater)); // Perfect matches are better
    Ok(())
}
#[test]
fn test_fragment_state_ordering() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let state1 = FragmentState::from_record(rec1, 0)?;
    let state2 = FragmentState::from_record(rec2, 1)?;
    assert_eq!(state1.partial_cmp(&state2), None);

    let mut ord: Option<Ordering> = None;
    let _ = state1.cmp_perfect(&state2, &mut ord)?;
    assert_eq!(ord, Some(Ordering::Equal));
    Ok(())
}
#[test]
fn test_fragment_state_partial_ord_with_unmapped() -> Result<(), Error> {
    let qual = vec![37; 100];
    let seq = vec![b'A'; 100];
    let rec1 = create_record(b"read1", "", &seq, &qual, "", false)?;
    let rec2 = create_record(b"read2", "100M", &seq, &qual, "", false)?;
    let state1 = FragmentState::from_record(rec1, 0)?;
    let state2 = FragmentState::from_record(rec2, 0)?;
    assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Less));
    Ok(())
}
#[test]
fn test_fragment_state_partial_ord_both_unmapped() -> Result<(), Error> {
    let qual = vec![37; 100];
    let seq = vec![b'A'; 100];
    let rec1 = create_record(b"read1", "", &seq, &qual, "", false)?;
    let rec2 = create_record(b"read2", "", &seq, &qual, "", false)?;
    let state1: FragmentState<RecordBuf> = FragmentState::from_record(rec1, 0)?;
    let state2: FragmentState<RecordBuf> = FragmentState::from_record(rec2, 0)?;
    assert_eq!(state1.partial_cmp(&state2), Some(Ordering::Equal));
    Ok(())
}
#[test]
fn test_fragment_state_partial_ord_perfect_vs_imperfect() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let state1 = FragmentState::from_record(rec1, 0)?;
    let state2 = FragmentState::from_record(rec2, 0)?;
    assert_eq!(state1.partial_cmp(&state2), None);

    let mut ord: Option<Ordering> = None;
    let _ = state1.cmp_perfect(&state2, &mut ord)?;
    assert_eq!(ord, Some(Ordering::Greater)); // Perfect match is better
    Ok(())
}
#[test]
fn test_fragment_state_partial_ord_imperfect_vs_perfect() -> Result<(), Error> {
    let qual = vec![37; 100];
    let rec1 = create_record(b"read1", "100M", &[], &qual, "90A10", false)?;
    let rec2 = create_record(b"read1", "100M", &[], &qual, "100", false)?;
    let state1 = FragmentState::from_record(rec1, 0)?;
    let state2 = FragmentState::from_record(rec2, 0)?;
    assert_eq!(state1.partial_cmp(&state2), None); // Perfect match is better
    let mut ord: Option<Ordering> = None;
    let _ = state1.cmp_perfect(&state2, &mut ord)?;
    assert_eq!(ord, Some(Ordering::Less)); // Perfect match is better
    Ok(())
}

#[test]
fn test_order_mates_sorts_by_segment_then_secondary_then_tid_then_pos() -> Result<(), Error> {
    let primary_first = segment(b"r", 0x40, 0, 100)?; // ord 0
    let secondary_first = segment(b"r", 0x40 | 0x100, 0, 200)?; // ord 1
    let primary_last = segment(b"r", 0x80, 1, 50)?; // ord 2
    let secondary_last = segment(b"r", 0x80 | 0x100, 1, 999)?; // ord 3

    // Insert deliberately out of order to verify sorting, not insertion order.
    let mut state = FragmentState::from_record(secondary_last, 0)?; // idx 0
    state.add_record(primary_last)?; // idx 1
    state.add_record(secondary_first)?; // idx 2
    state.add_record(primary_first)?; // idx 3

    let order = state.order_mates();
    let expected: SmallVec<[usize; 2]> = smallvec![3, 2, 1, 0];
    assert_eq!(order, expected);
    Ok(())
}

#[test]
fn test_fragment_state_flags_boundary_and_value() -> Result<(), Error> {
    // Use a non-default flag (e.g., 0x40 for first segment)
    let rec = segment(b"read", 0x40, 0, 100)?;
    let state = FragmentState::from_record(rec, 0)?;

    // Kills the leak/default mutant by ensuring valid index matches exact flags
    assert_eq!(state.flags(0), Some(&Flags::from_bits(0x40).unwrap()));

    // Kills the mutant by ensuring out-of-bounds correctly evaluates to None
    assert_eq!(state.flags(1), None);
    Ok(())
}

#[test]
fn test_order_mates_ord_math_isolation() -> Result<(), Error> {
    // Record A: is_last_segment = false, is_secondary = true -> valid ord = 1
    // We give it a higher tid/pos so if ord ties under mutation, it sorts incorrectly.
    let rec_a = segment(b"r", 0x100, 1, 500)?;

    // Record B: is_last_segment = true, is_secondary = false -> valid ord = 2
    // We give it a lower tid/pos so if ord ties under mutation, it sorts incorrectly.
    let rec_b = segment(b"r", 0x80, 0, 100)?;

    // Insert B first, then A
    let mut state = FragmentState::from_record(rec_b, 0)?;
    state.add_record(rec_a)?;

    let order = state.order_mates();
    // Expected: Record A (idx 1, ord 1) MUST sort before Record B (idx 0, ord 2)
    // If operators switch (+ <-> *), ord values tie, and tid forces [0, 1] instead.
    let expected: SmallVec<[usize; 2]> = smallvec![1, 0];
    assert_eq!(order, expected);
    Ok(())
}

#[test]
fn test_order_mates_reverse_complement_pos_math() -> Result<(), Error> {
    // We need values where (start1 + len1 < start2 + len2) BUT (start1 * len1 > start2 * len2)
    // Record 1: start = 5, cigar ops len = 5 ("1M1M1M1M1M") -> true pos = 10, mutated pos = 25
    let mut rec1 = create_record(b"r", "1M1M1M1M1M", &[], &[], "5", false)?;
    *rec1.flags_mut() = Flags::from_bits(0x10).unwrap(); // Reverse complemented
    *rec1.reference_sequence_id_mut() = Some(0);
    *rec1.alignment_start_mut() = Some(Position::new(5).unwrap());

    // Record 2: start = 1, cigar ops len = 10 ("1M1M1M1M1M1M1M1M1M1M") -> true pos = 11, mutated pos = 10
    let mut rec2 = create_record(b"r", "1M1M1M1M1M1M1M1M1M1M", &[], &[], "10", false)?;
    *rec2.flags_mut() = Flags::from_bits(0x10).unwrap(); // Reverse complemented
    *rec2.reference_sequence_id_mut() = Some(0);
    *rec2.alignment_start_mut() = Some(Position::new(1).unwrap());

    // Insert rec2 (idx 0) then rec1 (idx 1)
    let mut state = FragmentState::from_record(rec2, 0)?;
    state.add_record(rec1)?;

    let order = state.order_mates();
    // Expected: rec1 (idx 1, true pos 10) comes before rec2 (idx 0, true pos 11)
    let expected: SmallVec<[usize; 2]> = smallvec![1, 0];
    assert_eq!(order, expected);
    Ok(())
}
