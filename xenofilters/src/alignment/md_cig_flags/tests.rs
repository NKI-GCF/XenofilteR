// Table-driven style: each logical original test is expressed as one or more
// table-driven cases. Failures are collected and reported together at the end
// (with a human-readable case name derived from the original test name).
//
// These tests reproduce the assertions from the original file but run them
// data-driven and accumulate all misses so a single run shows everything that fails.

use super::*;
use crate::tests::create_record;
use noodles::sam::alignment::record::Flags;
use std::cmp::Ordering;

// ---------------------------------------------------------------------
// Group 1: is_perfect permutations (single-op clean, mismatch, multi-op)
// ---------------------------------------------------------------------
#[test]
fn table_is_perfect_variants_collect_misses() {
    let mut misses: Vec<String> = Vec::new();

    struct Case {
        name: &'static str,
        cigar: &'static str,
        md: &'static str,
        expect_perfect: bool,
    }

    let cases = vec![
        Case {
            name: "test_is_perfect_single_op_no_mismatches",
            cigar: "10M",
            md: "10",
            expect_perfect: true,
        },
        Case {
            name: "test_is_perfect_false_with_mismatch",
            cigar: "10M",
            md: "5A4",
            expect_perfect: false,
        },
        Case {
            name: "test_is_perfect_false_multiple_cigar_ops_even_if_md_clean",
            cigar: "5M5S",
            md: "5",
            expect_perfect: false,
        },
    ];

    for c in cases {
        match create_record(b"r", c.cigar, &[], &[], c.md, false) {
            Ok(rec) => {
                let flags = rec.flags();
                match MdCigFlags::try_from_record(&rec, &flags, false) {
                    Ok(mcf) => {
                        let got = mcf.is_perfect();
                        if got != c.expect_perfect {
                            misses.push(format!(
                                "{}: cigar='{}' md='{}' expected is_perfect={} got={}",
                                c.name, c.cigar, c.md, c.expect_perfect, got
                            ));
                        }
                    }
                    Err(e) => {
                        misses.push(format!(
                                "{}: MdCigFlags::try_from_record failed for cigar='{}' md='{}' err={:?}",
                                c.name, c.cigar, c.md, e
                            ));
                    }
                }
            }
            Err(e) => {
                misses.push(format!(
                    "{}: create_record failed for cigar='{}' md='{}' err={:?}",
                    c.name, c.cigar, c.md, e
                ));
            }
        }
    }

    if !misses.is_empty() {
        panic!("is_perfect table failures:\n\n{}", misses.join("\n"));
    }
}

// ---------------------------------------------------------------------
// Group 2: reverse complement & last-segment flag checks
// ---------------------------------------------------------------------
#[test]
fn table_revcomp_and_last_segment_collect_misses() {
    let mut misses: Vec<String> = Vec::new();

    // Reverse complement test (two cases combined to reduce duplication)
    match create_record(b"r", "5M", &[], &[], "5", false) {
        Ok(fwd) => match create_record(b"r", "5M", &[], &[], "5", true) {
            Ok(rev) => {
                match (
                        MdCigFlags::try_from_record(&fwd, &fwd.flags(), false),
                        MdCigFlags::try_from_record(&rev, &rev.flags(), false),
                    ) {
                        (Ok(mf), Ok(mr)) => {
                            if mf.is_reverse_complemented() {
                                misses.push("test_is_reverse_complemented: forward read unexpectedly reverse_complemented".to_string());
                            }
                            if !mr.is_reverse_complemented() {
                                misses.push("test_is_reverse_complemented: reverse read not reported reverse_complemented".to_string());
                            }
                        }
                        (Err(e), _) => misses.push(format!("test_is_reverse_complemented: MdCigFlags try_from_record failed for forward: {:?}", e)),
                        (_, Err(e)) => misses.push(format!("test_is_reverse_complemented: MdCigFlags try_from_record failed for reverse: {:?}", e)),
                    }
            }
            Err(e) => misses.push(format!(
                "test_is_reverse_complemented: create_record(rev) failed: {:?}",
                e
            )),
        },
        Err(e) => misses.push(format!(
            "test_is_reverse_complemented: create_record(fwd) failed: {:?}",
            e
        )),
    }

    // last/first segment flags
    // Build two records and set their Flags bits manually as in original test.
    match create_record(b"r", "5M", &[], &[], "5", false) {
        Ok(mut last) => {
            match create_record(b"r", "5M", &[], &[], "5", false) {
                Ok(mut first) => {
                    *last.flags_mut() = Flags::from_bits(0x80).unwrap(); // last segment mapped
                    *first.flags_mut() = Flags::from_bits(0x40).unwrap(); // first segment mapped
                    match (
                        MdCigFlags::try_from_record(&last, &last.flags(), false),
                        MdCigFlags::try_from_record(&first, &first.flags(), false),
                    ) {
                        (Ok(ml), Ok(mf)) => {
                            if !ml.is_last_segment() {
                                misses.push("test_is_last_segment: last record not reported as last_segment".to_string());
                            }
                            if mf.is_last_segment() {
                                misses.push("test_is_last_segment: first record incorrectly reported as last_segment".to_string());
                            }
                        }
                        (Err(e), _) => misses.push(format!(
                            "test_is_last_segment: try_from_record failed for last: {:?}",
                            e
                        )),
                        (_, Err(e)) => misses.push(format!(
                            "test_is_last_segment: try_from_record failed for first: {:?}",
                            e
                        )),
                    }
                }
                Err(e) => misses.push(format!(
                    "test_is_last_segment: create_record(first) failed: {:?}",
                    e
                )),
            }
        }
        Err(e) => misses.push(format!(
            "test_is_last_segment: create_record(last) failed: {:?}",
            e
        )),
    }

    if !misses.is_empty() {
        panic!(
            "reverse/last-segment table failures:\n\n{}",
            misses.join("\n")
        );
    }
}

// ---------------------------------------------------------------------
// Group 3: try_from_record error conditions (unmapped, missing MD)
// ---------------------------------------------------------------------
#[test]
fn table_try_from_record_errors_collect_misses() {
    let mut misses: Vec<String> = Vec::new();

    // unmapped should be rejected
    match create_record(b"r", "", &[], &[], "", false) {
        Ok(rec_unmapped) => {
            let flags = rec_unmapped.flags();
            if MdCigFlags::try_from_record(&rec_unmapped, &flags, false).is_ok() {
                misses.push(
                    "test_try_from_record_unmapped_is_rejected: expected error but got Ok"
                        .to_string(),
                );
            }
        }
        Err(e) => misses.push(format!(
            "test_try_from_record_unmapped_is_rejected: create_record failed: {:?}",
            e
        )),
    }

    // missing MD tag should error
    match create_record(b"r", "5M", &[], &[], "5", false) {
        Ok(mut rec) => {
            *rec.data_mut() = Default::default(); // strip MD tag
            let flags = rec.flags();
            if MdCigFlags::try_from_record(&rec, &flags, false).is_ok() {
                misses.push(
                    "test_try_from_record_missing_md_tag_errors: expected error but got Ok"
                        .to_string(),
                );
            }
        }
        Err(e) => misses.push(format!(
            "test_try_from_record_missing_md_tag_errors: create_record failed: {:?}",
            e
        )),
    }

    if !misses.is_empty() {
        panic!(
            "try_from_record error cases failures:\n\n{}",
            misses.join("\n")
        );
    }
}

// ---------------------------------------------------------------------
// Group 4: PartialEq and PartialOrd behavior (perfect vs imperfect)
// This is a single multi-step table-like test reproducing the original logic
// ---------------------------------------------------------------------
#[test]
fn partial_ord_and_eq_reproduced_collect_misses() {
    let mut misses: Vec<String> = Vec::new();

    // Setup records
    let p1 = create_record(b"r1", "10M", &[], &[], "10", false);
    let p2 = create_record(b"r2", "10M", &[], &[], "10", false);
    let i1 = create_record(b"r3", "10M", &[], &[], "5A4", false);
    let i2 = create_record(b"r4", "10M", &[], &[], "5A4", false);

    // If any creation failed, record as test failure and abort checks
    if p1.is_err() || p2.is_err() || i1.is_err() || i2.is_err() {
        misses.push(format!(
            "partial_ord_and_eq: create_record failed: p1={:?} p2={:?} i1={:?} i2={:?}",
            p1.err(),
            p2.err(),
            i1.err(),
            i2.err()
        ));
    } else {
        let p1 = p1.unwrap();
        let p2 = p2.unwrap();
        let i1 = i1.unwrap();
        let i2 = i2.unwrap();

        // Build MdCigFlags
        let flags_p1 = p1.flags();
        let flags_p2 = p2.flags();
        let flags_i1 = i1.flags();
        let flags_i2 = i2.flags();
        let md_p1 = MdCigFlags::try_from_record(&p1, &flags_p1, false);
        let md_p2 = MdCigFlags::try_from_record(&p2, &flags_p2, false);
        let md_i1 = MdCigFlags::try_from_record(&i1, &flags_i1, false);
        let md_i2 = MdCigFlags::try_from_record(&i2, &flags_i2, false);
        let _ = drop(flags_p1);
        let _ = drop(flags_p2);
        let _ = drop(flags_i1);
        let _ = drop(flags_i2);

        if md_p1.is_err() || md_p2.is_err() || md_i1.is_err() || md_i2.is_err() {
            misses.push(format!(
                    "partial_ord_and_eq: MdCigFlags construction failed: p1={:?} p2={:?} i1={:?} i2={:?}",
                    md_p1.err(), md_p2.err(), md_i1.err(), md_i2.err()
                ));
        } else {
            let md_p1 = md_p1.unwrap();
            let md_p2 = md_p2.unwrap();
            let md_i1 = md_i1.unwrap();
            let md_i2 = md_i2.unwrap();

            // --- PartialEq checks
            if !(md_p1 == md_p2) {
                misses.push("partial_ord_and_eq: perfect == perfect should be true".to_string());
            }
            if md_p1 == md_i1 {
                misses.push(
                    "partial_ord_and_eq: perfect != imperfect should be true (got equal)"
                        .to_string(),
                );
            }
            if md_i1 == md_i2 {
                misses.push("partial_ord_and_eq: imperfect == imperfect should be false (expected None semantics)".to_string());
            }

            // --- PartialOrd checks
            if md_p1.partial_cmp(&md_p2) != Some(Ordering::Equal) {
                misses.push(
                    "partial_ord_and_eq: md_p1.partial_cmp(md_p2) != Some(Equal)".to_string(),
                );
            }
            if md_p1.partial_cmp(&md_i1) != Some(Ordering::Less) {
                misses
                    .push("partial_ord_and_eq: md_p1.partial_cmp(md_i1) != Some(Less)".to_string());
            }
            if md_i1.partial_cmp(&md_p1) != Some(Ordering::Greater) {
                misses.push(
                    "partial_ord_and_eq: md_i1.partial_cmp(md_p1) != Some(Greater)".to_string(),
                );
            }
            if md_i1.partial_cmp(&md_i2) != None {
                misses.push("partial_ord_and_eq: md_i1.partial_cmp(md_i2) != None".to_string());
            }
        }
    }

    if !misses.is_empty() {
        panic!("partial_ord_and_eq failures:\n\n{}", misses.join("\n"));
    }
}
