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
    let _misses: Vec<String> = Vec::new();

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
    crate::tests::common::run_collecting(
        &cases,
        |c| c.name.to_string(),
        |c| {
            let rec = create_record(b"r", c.cigar, &[], &[], c.md, false)
                .map_err(|e| format!("create_record failed: {e:?}"))?;
            let flags = rec.flags();
            let mcf = MdCigFlags::try_from_record(&rec, &flags, false)
                .map_err(|e| format!("try_from_record failed: {e:?}"))?;
            let got = mcf.is_perfect();
            (got == c.expect_perfect).then_some(()).ok_or_else(|| {
                format!(
                    "cigar='{}' md='{}' expected {} got {}",
                    c.cigar, c.md, c.expect_perfect, got
                )
            })
        },
    );
}

// ---------------------------------------------------------------------
// Group 2: reverse complement & last-segment flag checks
// ---------------------------------------------------------------------
#[test]
fn table_revcomp_and_last_segment_collect_misses() {
    struct Case {
        name: &'static str,
    }
    let cases = vec![
        Case {
            name: "test_is_reverse_complemented",
        },
        Case {
            name: "test_is_last_segment",
        },
    ];

    crate::tests::common::run_collecting(
        &cases,
        |c| c.name.to_string(),
        |c| {
            if c.name == "test_is_reverse_complemented" {
                let fwd = create_record(b"r", "5M", &[], &[], "5", false)
                    .map_err(|e| format!("create_record(fwd) failed: {:?}", e))?;
                let rev = create_record(b"r", "5M", &[], &[], "5", true)
                    .map_err(|e| format!("create_record(rev) failed: {:?}", e))?;
                let fwd_flags = fwd.flags();
                let rev_flags = rev.flags();
                let mf = MdCigFlags::try_from_record(&fwd, &fwd_flags, false).map_err(|e| {
                    format!("MdCigFlags try_from_record failed for forward: {:?}", e)
                })?;
                let mr = MdCigFlags::try_from_record(&rev, &rev_flags, false).map_err(|e| {
                    format!("MdCigFlags try_from_record failed for reverse: {:?}", e)
                })?;

                if mf.is_reverse_complemented() {
                    return Err("forward read unexpectedly reverse_complemented".to_string());
                }
                if !mr.is_reverse_complemented() {
                    return Err("reverse read not reported reverse_complemented".to_string());
                }
                Ok(())
            } else {
                let mut last = create_record(b"r", "5M", &[], &[], "5", false)
                    .map_err(|e| format!("create_record(last) failed: {:?}", e))?;
                let mut first = create_record(b"r", "5M", &[], &[], "5", false)
                    .map_err(|e| format!("create_record(first) failed: {:?}", e))?;

                *last.flags_mut() = Flags::from_bits(0x80).unwrap(); // last segment mapped
                *first.flags_mut() = Flags::from_bits(0x40).unwrap(); // first segment mapped
                let last_flags = last.flags();
                let first_flags = first.flags();
                let ml = MdCigFlags::try_from_record(&last, &last_flags, false)
                    .map_err(|e| format!("try_from_record failed for last: {:?}", e))?;
                let mf = MdCigFlags::try_from_record(&first, &first_flags, false)
                    .map_err(|e| format!("try_from_record failed for first: {:?}", e))?;

                if !ml.is_last_segment() {
                    return Err("last record not reported as last_segment".to_string());
                }
                if mf.is_last_segment() {
                    return Err("first record incorrectly reported as last_segment".to_string());
                }
                Ok(())
            }
        },
    );
}

// ---------------------------------------------------------------------
// Group 3: try_from_record error conditions (unmapped, missing MD)
// ---------------------------------------------------------------------
#[test]
fn table_try_from_record_errors_collect_misses() {
    struct Case {
        name: &'static str,
    }
    let cases = vec![
        Case {
            name: "test_try_from_record_unmapped_is_rejected",
        },
        Case {
            name: "test_try_from_record_missing_md_tag_errors",
        },
    ];

    crate::tests::common::run_collecting(
        &cases,
        |c| c.name.to_string(),
        |c| {
            if c.name == "test_try_from_record_unmapped_is_rejected" {
                let rec_unmapped = create_record(b"r", "", &[], &[], "", false)
                    .map_err(|e| format!("create_record failed: {:?}", e))?;
                let flags = rec_unmapped.flags();
                if MdCigFlags::try_from_record(&rec_unmapped, &flags, false).is_ok() {
                    return Err("expected error but got Ok".to_string());
                }
                Ok(())
            } else {
                let mut rec = create_record(b"r", "5M", &[], &[], "5", false)
                    .map_err(|e| format!("create_record failed: {:?}", e))?;
                *rec.data_mut() = Default::default(); // strip MD tag
                let flags = rec.flags();
                if MdCigFlags::try_from_record(&rec, &flags, false).is_ok() {
                    return Err("expected error but got Ok".to_string());
                }
                Ok(())
            }
        },
    );
}

// ---------------------------------------------------------------------
// Group 4: PartialEq and PartialOrd behavior (perfect vs imperfect)
// This is a single multi-step table-like test reproducing the original logic
// ---------------------------------------------------------------------
#[test]
fn partial_ord_and_eq_reproduced_collect_misses() {
    struct Case {
        name: &'static str,
        check: fn(&MdCigFlags, &MdCigFlags, &MdCigFlags, &MdCigFlags) -> Result<(), String>,
    }

    let cases = vec![
        Case {
            name: "perfect == perfect",
            check: |p1, p2, _, _| {
                (p1 == p2)
                    .then_some(())
                    .ok_or_else(|| "perfect == perfect should be true".to_string())
            },
        },
        Case {
            name: "perfect != imperfect",
            check: |p1, _, i1, _| {
                (p1 != i1)
                    .then_some(())
                    .ok_or_else(|| "perfect != imperfect should be true (got equal)".to_string())
            },
        },
        Case {
            name: "imperfect == imperfect",
            check: |_, _, i1, i2| {
                (i1 != i2).then_some(()).ok_or_else(|| {
                    "imperfect == imperfect should be false (expected None semantics)".to_string()
                })
            },
        },
        Case {
            name: "md_p1.partial_cmp(md_p2) == Some(Equal)",
            check: |p1, p2, _, _| {
                (p1.partial_cmp(p2) == Some(Ordering::Equal))
                    .then_some(())
                    .ok_or_else(|| "md_p1.partial_cmp(md_p2) != Some(Equal)".to_string())
            },
        },
        Case {
            name: "md_p1.partial_cmp(md_i1) == Some(Greater)",
            check: |p1, _, i1, _| {
                (p1.partial_cmp(i1) == Some(Ordering::Greater))
                    .then_some(())
                    .ok_or_else(|| "md_p1.partial_cmp(md_i1) != Some(Greater)".to_string())
            },
        },
        Case {
            name: "md_i1.partial_cmp(md_p1) == Some(Less)",
            check: |p1, _, i1, _| {
                (i1.partial_cmp(p1) == Some(Ordering::Less))
                    .then_some(())
                    .ok_or_else(|| "md_i1.partial_cmp(md_p1) != Some(Less)".to_string())
            },
        },
        Case {
            name: "md_i1.partial_cmp(md_i2) == None",
            check: |_, _, i1, i2| {
                i1.partial_cmp(i2)
                    .is_none()
                    .then_some(())
                    .ok_or_else(|| "md_i1.partial_cmp(md_i2) != None".to_string())
            },
        },
    ];

    let p1 = create_record(b"r1", "10M", &[], &[], "10", false).unwrap();
    let p2 = create_record(b"r2", "10M", &[], &[], "10", false).unwrap();
    let i1 = create_record(b"r3", "10M", &[], &[], "5A4", false).unwrap();
    let i2 = create_record(b"r4", "10M", &[], &[], "5A4", false).unwrap();

    let flags_p1 = p1.flags();
    let flags_p2 = p2.flags();
    let flags_i1 = i1.flags();
    let flags_i2 = i2.flags();

    let md_p1 = MdCigFlags::try_from_record(&p1, &flags_p1, false).unwrap();
    let md_p2 = MdCigFlags::try_from_record(&p2, &flags_p2, false).unwrap();
    let md_i1 = MdCigFlags::try_from_record(&i1, &flags_i1, false).unwrap();
    let md_i2 = MdCigFlags::try_from_record(&i2, &flags_i2, false).unwrap();

    crate::tests::common::run_collecting(
        &cases,
        |c| c.name.to_string(),
        |c| (c.check)(&md_p1, &md_p2, &md_i1, &md_i2),
    );
}
