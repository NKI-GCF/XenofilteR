use crate::alignment::{BaseOp, MdCigFlags, ScoreOpIter};
use crate::tests::create_record;
use crate::Error;

fn op_repr(op: &BaseOp) -> String {
    match op {
        BaseOp::Match => "M".into(),
        BaseOp::Mis => "X".into(),
        BaseOp::Del(n) => format!("D{n}"),
        BaseOp::Ins(n) => format!("I{n}"),
        BaseOp::Clip(n) => format!("C{n}"),
        BaseOp::RefSkip(n) => format!("N{n}"),
        BaseOp::BisulfiteConversion => "B".into(),
    }
}

fn ops_for(cigar: &str, md: &str) -> Result<Vec<String>, Error> {
    let rec = create_record(b"r", cigar, &[], &[], md, false)?;
    let flags = rec.flags();
    let mcf = MdCigFlags::try_from_record(&rec, &flags, false)?;
    let ops: Result<Vec<BaseOp>, _> = ScoreOpIter::new(&mcf).collect();
    //Ok(ops.iter().map(op_repr).collect())
    Ok(ops?.iter().map(op_repr).collect())
}

// Table-driven unified test that replaces multiple small test functions.
#[test]
fn table_driven_ops_tests_collect_misses() {
    enum Expectation {
        Exact(Vec<&'static str>),
        Len(usize),
        Err,
    }

    struct Case {
        msg: &'static str,
        cigar: &'static str,
        md: &'static str,
        expect: Expectation,
    }

    let cases = vec![
        Case {
            msg: "pure match",
            cigar: "5M",
            md: "5",
            expect: Expectation::Exact(vec!["M", "M", "M", "M", "M"]),
        },
        Case {
            msg: "single mismatch",
            cigar: "5M",
            md: "2A2",
            expect: Expectation::Exact(vec!["M", "M", "X", "M", "M"]),
        },
        Case {
            msg: "consecutive mismatches",
            cigar: "5M",
            md: "1AC2",
            expect: Expectation::Exact(vec!["M", "X", "X", "M", "M"]),
        },
        Case {
            msg: "insertion not consumed by md and grouped",
            cigar: "2M2I2M",
            md: "4",
            expect: Expectation::Exact(vec!["M", "M", "I2", "M", "M"]),
        },
        Case {
            msg: "deletion consumes caret block",
            cigar: "2M3D2M",
            md: "2^AAA2",
            expect: Expectation::Exact(vec!["M", "M", "D3", "M", "M"]),
        },
        Case {
            msg: "soft clip is single grouped op",
            cigar: "3S5M",
            md: "5",
            expect: Expectation::Exact(vec!["C3", "M", "M", "M", "M", "M"]),
        },
        Case {
            msg: "ref skip not consumed by md",
            cigar: "2M3N2M",
            md: "4",
            expect: Expectation::Exact(vec!["M", "M", "N3", "M", "M"]),
        },
        Case {
            msg: "hard clip and pad are invisible to iterator 1",
            cigar: "2H5M2H",
            md: "5",
            expect: Expectation::Exact(vec!["M", "M", "M", "M", "M"]),
        },
        Case {
            msg: "hard clip and pad are invisible to iterator 2",
            cigar: "2M2P2M",
            md: "4",
            expect: Expectation::Exact(vec!["M", "M", "M", "M"]),
        },
        Case {
            msg: "sequence match mismatch kinds route through md",
            cigar: "5=",
            md: "5",
            expect: Expectation::Exact(vec!["M", "M", "M", "M", "M"]),
        },
        Case {
            msg: "multidigit md run length",
            cigar: "12M",
            md: "12",
            expect: Expectation::Len(12),
        },
        Case {
            msg: "deletion length mismatch is error",
            cigar: "2M3D2M",
            md: "2^AA2",
            expect: Expectation::Err,
        },
        Case {
            msg: "invalid md character is error",
            cigar: "4M",
            md: "2x2",
            expect: Expectation::Err,
        },
    ];

    crate::tests::common::run_collecting(
        &cases,
        |c| c.msg.to_string(),
        |c| match ops_for(c.cigar, c.md) {
            Ok(got) => match &c.expect {
                Expectation::Exact(evec) => {
                    let expected: Vec<String> = evec.iter().map(|s| s.to_string()).collect();
                    if got != expected {
                        return Err(format!(
                            "input: cigar='{}' md='{}'\n  expected: {:?}\n  got     : {:?}",
                            c.cigar, c.md, expected, got
                        ));
                    }
                    Ok(())
                }
                Expectation::Len(n) => {
                    if got.len() != *n {
                        return Err(format!(
                            "input: cigar='{}' md='{}'\n  expected len: {}\n  got len     : {}",
                            c.cigar,
                            c.md,
                            n,
                            got.len()
                        ));
                    }
                    Ok(())
                }
                Expectation::Err => Err(format!(
                    "expected error but got Ok:\n  input: cigar='{}' md='{}'\n  got: {:?}",
                    c.cigar, c.md, got
                )),
            },
            Err(e) => match c.expect {
                Expectation::Err => Ok(()),
                _ => Err(format!(
                    "unexpected error:\n  input: cigar='{}' md='{}'\n  error: {:?}",
                    c.cigar, c.md, e
                )),
            },
        },
    );
}

#[test]
fn bisulfite_forward_conversion_table() {
    // Forward strand: ref=C, read=T → BisulfiteConversion; other combos → Mis
    struct Row {
        ref_b: u8,
        read_b: u8,
        is_rev: bool,
        want: BaseOp,
    }
    let cases = [
        Row {
            ref_b: b'C',
            read_b: b'T',
            is_rev: false,
            want: BaseOp::BisulfiteConversion,
        },
        Row {
            ref_b: b'C',
            read_b: b'A',
            is_rev: false,
            want: BaseOp::Mis,
        },
        Row {
            ref_b: b'G',
            read_b: b'A',
            is_rev: true,
            want: BaseOp::BisulfiteConversion,
        },
        Row {
            ref_b: b'G',
            read_b: b'T',
            is_rev: true,
            want: BaseOp::Mis,
        },
        Row {
            ref_b: b'A',
            read_b: b'T',
            is_rev: false,
            want: BaseOp::Mis,
        },
        Row {
            ref_b: b'C',
            read_b: b'T',
            is_rev: true,
            want: BaseOp::Mis,
        },
    ];
    for c in &cases {
        // Use a minimal ScoreOpIter::classify_mismatch harness.
        let seq = &[c.read_b];
        let mut it = ScoreOpIter {
            seq: Some(seq),
            is_reverse: c.is_rev,
            read_pos: 0,
            ..Default::default()
        };
        let got = it.classify_mismatch(c.ref_b);
        assert_eq!(
            got, c.want,
            "ref={} read={} is_rev={}",
            c.ref_b as char, c.read_b as char, c.is_rev
        );
    }
}

#[test]
fn error_model_quality_calibration() {
    use crate::penalty::{ErrorModel, Penalty};
    let ill = Penalty::build(6.0, 1.0, 4.0, 20, ErrorModel::Illumina);
    let ont = Penalty::build(2.0, 0.3, 2.0, 20, ErrorModel::Ont);
    // At Q30, ONT should trust the quality less → higher (less negative) log_lik_mismatch
    // i.e. the mismatch is less penalised for ONT at the same reported Q.
    assert!(
        ont.log_likelihood_mismatch[30] > ill.log_likelihood_mismatch[30],
        "ONT Q30 should be less penalised than Illumina Q30 (calibration factor 0.7)"
    );
}

#[test]
fn error_model_gap_defaults() {
    use crate::penalty::ErrorModel;
    assert!(
        ErrorModel::Ont.default_gap_open() < ErrorModel::Illumina.default_gap_open(),
        "ONT gap_open should be lower than Illumina (indel-prone chemistry)"
    );
    assert!(
        ErrorModel::HiFi.default_gap_open() < ErrorModel::Illumina.default_gap_open(),
        "HiFi gap_open should be lower (few indels in CCS)"
    );
}
