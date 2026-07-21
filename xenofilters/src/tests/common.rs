//! Parametric test infrastructure shared across all test modules.
//! Import with `use crate::tests::common::*;`

use crate::config::args::ScoringArgs;
use crate::config::run_config::RunConfig;
use crate::tests::create_record;
use noodles::sam::alignment::record_buf::RecordBuf;

pub(crate) fn cfg() -> RunConfig {
    RunConfig {
        scoring: ScoringArgs {
            gap_open: 6.0,
            gap_extend: 1.0,
            mismatch_penalty: 4.0,
            ..ScoringArgs::default()
        },
        ..RunConfig::default()
    }
}

pub(crate) fn r(name: &[u8], cigar: &str, md: &str) -> RecordBuf {
    create_record(name, cigar, &[], &[30u8; 20], md, false).unwrap()
}
pub(crate) fn u(name: &[u8]) -> RecordBuf {
    create_record(name, "", &vec![b'A'; 10], &[30u8; 10], "", false).unwrap()
}

/// A single row in a table-driven test.
pub(crate) struct Case<I, O> {
    pub(crate) label: &'static str,
    pub(crate) input: I,
    pub(crate) want: O,
}

/// Run every case; panics with label on mismatch.
pub(crate) fn run<I, O, F>(cases: &[Case<I, O>], f: F)
where
    F: Fn(&I) -> O,
    O: PartialEq + std::fmt::Debug,
{
    for c in cases {
        let got = f(&c.input);
        assert_eq!(got, c.want, "[{}]", c.label);
    }
}

/// Convenience: build a mapped record.
pub(crate) fn mapped(cigar: &str, md: &str, qual: u8, len: usize) -> RecordBuf {
    let q = vec![qual; len];
    create_record(b"r", cigar, &[], &q, md, false).unwrap()
}

/// Convenience: unmapped record (empty cigar).
pub(crate) fn unmapped(len: usize) -> RecordBuf {
    let seq = vec![b'A'; len];
    let q = vec![30u8; len];
    create_record(b"r", "", &seq, &q, "", false).unwrap()
}

/// Table-driven test runner: collects all failing cases and panics once
/// with every mismatch, instead of the manual `Vec<String>` + panic
/// boilerplate duplicated across table_*_collect_misses tests.
///
/// `check` returns `Err(reason)` on mismatch, `Ok(())` on pass.
pub(crate) fn run_collecting<C>(
    cases: &[C],
    label: impl Fn(&C) -> String,
    check: impl Fn(&C) -> Result<(), String>,
) {
    let misses: Vec<String> = cases
        .iter()
        .filter_map(|c| check(c).err().map(|e| format!("[{}] {e}", label(c))))
        .collect();
    if !misses.is_empty() {
        panic!(
            "{} of {} cases failed:\n\n{}",
            misses.len(),
            cases.len(),
            misses.join("\n")
        );
    }
}
