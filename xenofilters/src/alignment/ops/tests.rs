use crate::alignment::{BaseOp, MdCigFlags, ScoreOpIter};
use crate::Error;
use noodles::core::Position;
use noodles::sam::alignment::{
    record::data::field::Tag,
    record_buf::{data::field::Value, Data},
};
use noodles::sam::alignment::{
    record::{
        cigar::{op::Kind, Op},
        Flags,
    },
    record_buf::{Cigar, QualityScores, RecordBuf, Sequence},
};
use std::iter::repeat;

pub(crate) fn create_cigar(cigar: &str) -> Result<Cigar, Error> {
    let mut ops = Vec::new();
    let mut num = 0;
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as u8 - b'0') as usize;
        } else {
            let kind = match c {
                'M' => Kind::Match,
                'I' => Kind::Insertion,
                'D' => Kind::Deletion,
                'N' => Kind::Skip,
                'S' => Kind::SoftClip,
                'H' => Kind::HardClip,
                'P' => Kind::Pad,
                '=' => Kind::SequenceMatch,
                'X' => Kind::SequenceMismatch,
                _ => return Err(Error::UnknownCigarOp(c as u32)),
            };
            ops.push(Op::new(kind, num));
            num = 0;
        }
    }
    Ok(Cigar::from(ops))
}

pub(crate) fn read_len_from_cigar(cigar: &str) -> usize {
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

pub(crate) fn create_record(
    qname: &[u8],
    cig_str: &str,
    seq: &[u8],
    qual: &[u8],
    md: &str,
    is_rev: bool,
) -> Result<RecordBuf, Error> {
    let mut record = RecordBuf::default();
    let is_unmapped = cig_str.is_empty() || cig_str == "*";

    *record.name_mut() = Some(qname.into());
    let read_len = if is_unmapped {
        *record.flags_mut() = Flags::from_bits(4).unwrap(); // unmapped
        0
    } else {
        if is_rev {
            *record.flags_mut() = Flags::from_bits(16).unwrap();
        } else {
            *record.flags_mut() = Flags::empty();
        }
        let data: Data = [(Tag::MISMATCHED_POSITIONS, Value::from(md))]
            .into_iter()
            .collect();
        *record.data_mut() = data;
        *record.cigar_mut() = create_cigar(cig_str)?;
        *record.alignment_start_mut() = Some(Position::MIN);
        *record.reference_sequence_id_mut() = Some(0);
        read_len_from_cigar(cig_str)
    };
    *record.sequence_mut() = if seq.is_empty() {
        Sequence::from(repeat(b'A').take(read_len).collect::<Vec<u8>>())
    } else {
        Sequence::from(seq)
    };

    *record.quality_scores_mut() = if qual.is_empty() {
        QualityScores::from_iter(repeat(30u8).take(read_len).collect::<Vec<u8>>())
    } else {
        QualityScores::from_iter(qual.iter().cloned())
    };

    Ok(record)
}

fn op_repr(op: &BaseOp) -> String {
    match op {
        BaseOp::Match => "M".into(),
        BaseOp::Mis => "X".into(),
        BaseOp::Del(n) => format!("D{n}"),
        BaseOp::Ins(n) => format!("I{n}"),
        BaseOp::Clip(n) => format!("C{n}"),
        BaseOp::RefSkip(n) => format!("N{n}"),
    }
}

fn ops_for(cigar: &str, md: &str) -> Result<Vec<String>, Error> {
    let rec = create_record(b"r", cigar, &[], &[], md, false)?;
    let flags = rec.flags();
    let mcf = MdCigFlags::try_from_record(&rec, &flags)?;
    let ops: Result<Vec<BaseOp>, _> = ScoreOpIter::new(&mcf).collect();
    Ok(ops
        .map_err(|e| anyhow::anyhow!("{e}"))?
        .iter()
        .map(op_repr)
        .collect())
}

// Table-driven unified test that replaces multiple small test functions.
#[test]
fn table_driven_ops_tests_collect_misses() {
    use std::fmt::Write as _;

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

    let mut misses: Vec<String> = Vec::new();

    for c in cases {
        match ops_for(c.cigar, c.md) {
            Ok(got) => match c.expect {
                Expectation::Exact(ref evec) => {
                    let expected: Vec<String> = evec.iter().map(|s| s.to_string()).collect();
                    if got != expected {
                        let mut s = String::new();
                        writeln!(&mut s, "{} failed:", c.msg).ok();
                        writeln!(&mut s, "  input: cigar='{}' md='{}'", c.cigar, c.md).ok();
                        writeln!(&mut s, "  expected: {:?}", expected).ok();
                        writeln!(&mut s, "  got     : {:?}", got).ok();
                        misses.push(s);
                    }
                }
                Expectation::Len(n) => {
                    if got.len() != n {
                        let mut s = String::new();
                        writeln!(&mut s, "{} failed (len):", c.msg).ok();
                        writeln!(&mut s, "  input: cigar='{}' md='{}'", c.cigar, c.md).ok();
                        writeln!(&mut s, "  expected len: {}", n).ok();
                        writeln!(&mut s, "  got len     : {}", got.len()).ok();
                        misses.push(s);
                    }
                }
                Expectation::Err => {
                    let mut s = String::new();
                    writeln!(&mut s, "{} failed (expected error but got Ok):", c.msg).ok();
                    writeln!(&mut s, "  input: cigar='{}' md='{}'", c.cigar, c.md).ok();
                    writeln!(&mut s, "  got: {:?}", got).ok();
                    misses.push(s);
                }
            },
            Err(e) => {
                match c.expect {
                    Expectation::Err => { /* expected error -> OK */ }
                    _ => {
                        let mut s = String::new();
                        writeln!(&mut s, "{} failed (unexpected error):", c.msg).ok();
                        writeln!(&mut s, "  input: cigar='{}' md='{}'", c.cigar, c.md).ok();
                        writeln!(&mut s, "  error: {:?}", e).ok();
                        misses.push(s);
                    }
                }
            }
        }
    }

    if !misses.is_empty() {
        panic!("Table-driven ops tests failures:\n\n{}", misses.join("\n"));
    }
}
