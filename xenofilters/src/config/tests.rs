use super::*;

fn base_config() -> Config {
    Config {
        alignment: vec!["a.bam".into(), "b.bam".into()],
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Default::default()
    }
}

#[test]
fn test_validate_ok_minimal() {
    assert!(base_config().validate_and_init().is_ok());
}

#[test]
fn test_validate_rejects_single_alignment() {
    let mut c = base_config();
    c.alignment = vec!["a.bam".into()];
    assert!(c.validate_and_init().is_err());
}

#[test]
fn test_validate_rejects_too_many_outputs() {
    let mut c = base_config();
    c.output = vec!["o1".into(), "o2".into(), "o3".into()];
    assert!(c.validate_and_init().is_err());
}

#[test]
fn test_validate_rejects_too_many_filtered_outputs() {
    let mut c = base_config();
    c.filtered_output = vec!["f1".into(), "f2".into(), "f3".into()];
    assert!(c.validate_and_init().is_err());
}

#[test]
fn test_validate_rejects_too_many_ambiguous_outputs() {
    let mut c = base_config();
    c.ambiguous_output = vec!["a1".into(), "a2".into(), "a3".into()];
    assert!(c.validate_and_init().is_err());
}

#[test]
fn test_validate_flips_positive_gap_penalties_to_negative() {
    let mut c = base_config();
    c.validate_and_init().unwrap();
    assert_eq!(c.gap_open, -6.0);
    assert_eq!(c.gap_extend, -1.0);
}

#[test]
fn test_validate_leaves_negative_gap_penalties_unchanged() {
    let mut c = base_config();
    c.gap_open = -6.0;
    c.gap_extend = -1.0;
    c.validate_and_init().unwrap();
    assert_eq!(c.gap_open, -6.0);
    assert_eq!(c.gap_extend, -1.0);
}

#[test]
fn test_validate_rejects_zero_gap_open() {
    let mut c = base_config();
    c.gap_open = 0.0;
    assert!(c.validate_and_init().is_err());
}

#[test]
fn test_validate_rejects_nonpositive_mismatch_penalty() {
    let mut c = base_config();
    c.mismatch_penalty = 0.0;
    assert!(c.validate_and_init().is_err());
    c.mismatch_penalty = -1.0;
    assert!(c.validate_and_init().is_err());
}

#[test]
fn test_validate_rejects_plain_default_config() {
    // Default::default() gives gap_open = gap_extend = mismatch_penalty = 0.0
    // (clap's `default_value` only applies via Parser::parse(), not Default).
    let mut c = Config {
        alignment: vec!["a".into(), "b".into()],
        ..Default::default()
    };
    assert!(c.validate_and_init().is_err());
}

/// BUG (documented, not just a coverage gap): `read_from_stdin` requires
/// `alignment.len() >= 2` (one ensure!) AND `alignment.len() == 1` (the stdin
/// ensure!) simultaneously. These can never both hold, so `--read-from-stdin`
/// is currently rejected unconditionally by validate_and_init, regardless of
/// alignment count.
#[test]
fn test_read_from_stdin_is_currently_always_rejected() {
    let mut c = base_config(); // 2 alignments
    c.read_from_stdin = true;
    assert!(c.validate_and_init().is_err());

    let mut c2 = base_config();
    c2.alignment = vec!["a.bam".into()]; // 1 alignment
    c2.read_from_stdin = true;
    assert!(c2.validate_and_init().is_err()); // fails the ">=2" check instead
}

#[test]
fn test_to_penalties_q0_is_certain_error() {
    let p = base_config().to_penalties();
    assert_eq!(p.gap_open, 6.0); // to_penalties doesn't flip sign itself
    assert!(p.log_likelihood_match[0].is_infinite()); // error_prob=1.0 -> log10(0)
}

#[test]
fn test_to_penalties_mismatch_scaling_factor() {
    let mut c = base_config();
    c.mismatch_penalty = 8.0; // scaling_factor = 8/4 = 2.0
    let p = c.to_penalties();
    assert!((p.log_likelihood_mismatch[10] - (-2.0)).abs() < 1e-9);

    c.mismatch_penalty = 4.0; // scaling_factor = 1.0
    let p = c.to_penalties();
    assert!((p.log_likelihood_mismatch[10] - (-1.0)).abs() < 1e-9);
}

#[test]
fn test_to_penalties_match_likelihood_improves_with_quality() {
    let p = base_config().to_penalties();
    assert!(p.log_likelihood_match[40] > p.log_likelihood_match[10]);
}
