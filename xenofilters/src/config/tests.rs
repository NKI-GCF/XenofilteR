use super::*;
use std::path::PathBuf;

fn base_config() -> Config {
    Config {
        alignment: vec!["a.bam".into(), "b.bam".into()],
        output: vec!["out1.bam".into(), "out2.bam".into()],
        ambiguous_output: vec![],
        sample_variants: vec![],
        population_variants: vec![],
        single_alignment_mode: false,
        read_from_stdin: false,
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

/// Fixes the failing test. `read_from_stdin` with 2 alignments is now perfectly valid
/// (returns Ok), but reading from stdin with 1 alignment under single-alignment mode is rejected.
#[test]
fn test_read_from_stdin_constraints() {
    // Multi-alignment standard streaming is fine
    let mut c = base_config();
    c.read_from_stdin = true;
    assert!(
        c.validate_and_init().is_ok(),
        "Standard multi-stream from stdin should pass validation"
    );

    // Single-alignment mode from stdin is strictly banned (cannot read twice)
    let mut c2 = base_config();
    c2.alignment = vec!["a.bam".into()];
    c2.output = vec!["out1.bam".into(), "out2.bam".into()];
    c2.single_alignment_mode = true;
    c2.read_from_stdin = true;
    c2.sample_variants = vec!["0:var1.vcf".into()];
    c2.population_variants = vec!["1:var2.vcf".into()];

    let res = c2.validate_and_init();
    assert!(res.is_err(), "Single alignment mode must reject stdin");
}

#[test]
fn test_single_alignment_mode_requires_flag() {
    let mut c = base_config();
    c.alignment = vec!["single_strain.bam".into()];
    c.output = vec!["out1.bam".into()];

    let res = c.validate_and_init();
    assert!(res.is_err());
    assert!(
        res.unwrap_err()
            .to_string()
            .contains("--single-alignment-mode"),
        "Should fail with explicit flag reminder string"
    );
}
#[test]
fn test_single_alignment_mode_variant_assertions() {
    // Case A: Missing variants entirely (0 variants total)
    let mut c = base_config();
    c.alignment = vec!["single_strain.bam".into()];
    c.single_alignment_mode = true;

    let res = c.validate_and_init();
    assert!(res.is_err());

    // Case B: Only 1 variant profile given (insufficient to differentiate strains)
    let mut c = base_config();
    c.alignment = vec!["single_strain.bam".into()];
    c.single_alignment_mode = true;
    c.sample_variants = vec!["0:strain_a.vcf".into()];

    assert!(c.validate_and_init().is_err());

    // Case C: Mix match of 1 sample and 1 population variant profile (Total = 2, Valid!)
    let mut c = base_config();
    c.alignment = vec!["single_strain.bam".into()];
    c.single_alignment_mode = true;
    c.sample_variants = vec!["0:strain_a.vcf".into()];
    c.population_variants = vec!["1:strain_b.vcf".into()];

    assert!(
        c.validate_and_init().is_ok(),
        "1 sample + 1 population profile is fully valid"
    );
}

#[test]
fn test_single_alignment_ambiguous_output_cap() {
    let mut c = base_config();
    c.alignment = vec!["single_strain.bam".into()];
    c.single_alignment_mode = true;
    c.sample_variants = vec!["0:a.vcf".into(), "1:b.vcf".into()];

    // Attaching 2 ambiguous outputs for single-stream virtual splits is forbidden
    c.ambiguous_output = vec!["ambig1.bam".into(), "ambig2.bam".into()];

    let res = c.validate_and_init();
    assert!(res.is_err());
}

#[test]
fn test_flag_mismatch_on_multi_alignment() {
    let mut c = base_config(); // Has 2 alignments
    c.single_alignment_mode = true; // User passed the flag accidentally

    let res = c.validate_and_init();
    assert!(res.is_err());
}

#[test]
fn test_invalid_variant_index_grouping_on_single_alignment() {
    let mut c = base_config();
    c.alignment = vec!["single_strain.bam".into()];
    c.single_alignment_mode = true;

    // Error case: both variations are targeting strain index 0. Strain index 1 remains empty.
    c.sample_variants = vec!["0:file1.vcf".into()];
    c.population_variants = vec!["0:pop_snps.vcf".into()];

    let res = c.validate_and_init();
    assert!(res.is_err());
}

#[test]
fn test_variant_array_padding_normalization() {
    let mut c = base_config();
    c.alignment = vec!["single_strain.bam".into()];
    c.single_alignment_mode = true;

    // Provide sample variant only for index 0 and population variant only for index 1
    c.sample_variants = vec!["0:varA.vcf".into()];
    c.population_variants = vec!["1:varB.vcf".into()];

    assert!(c.validate_and_init().is_ok());

    // Assert structural backfilling:
    // Slot 0 sample variant exists, slot 1 sample variant must be empty string
    assert_eq!(c.sample_variants[0], "varA.vcf");
    assert_eq!(c.sample_variants[1], "");

    // Slot 0 pop variant must be empty string, slot 1 pop variant exists
    assert_eq!(c.population_variants[0], "");
    assert_eq!(c.population_variants[1], "varB.vcf");
}

fn valid_base_config() -> Config {
    Config {
        alignment: vec!["aln1.bam".to_string(), "aln2.bam".to_string()],
        gap_open: 6.0,
        mismatch_penalty: 4.0,
        ..Default::default()
    }
}

#[test]
fn test_validate_non_namesorted_allows_only_one_score_thread() {
    let mut c = base_config();
    // Entering the non-namesorted branch
    c.matching_algorithm = MatchingAlgorithm::Collated;
    c.score_threads = 1;

    assert!(
        c.validate_and_init().is_ok(),
        "Collated matching algorithm with exactly 1 score thread must be valid"
    );

    // Ensure > 1 threads correctly trigger the validation error
    c.score_threads = 2;
    assert!(
        c.validate_and_init().is_err(),
        "Collated matching algorithm with > 1 score threads must be invalid"
    );
}

#[test]
fn test_zero_gap_penalties_do_not_flip_sign() {
    let mut c = base_config();
    c.gap_open = 0.0;

    // Validation will error out later in the function because gap_open == 0.0 is invalid,
    // but the mutant's sign-flipping side effect would have already executed.
    let _ = c.validate_and_init();
    assert!(
        !c.gap_open.is_sign_negative(),
        "Mutant killed: gap_open > 0.0 changed to >= 0.0 (flipped 0.0 to -0.0)"
    );

    let mut c2 = base_config();
    c2.gap_extend = 0.0; // A gap_extend of 0.0 is perfectly valid

    assert!(
        c2.validate_and_init().is_ok(),
        "gap_extend = 0.0 is valid and should pass initialization"
    );
    assert!(
        !c2.gap_extend.is_sign_negative(),
        "Mutant killed: gap_extend > 0.0 changed to >= 0.0 (flipped 0.0 to -0.0)"
    );
}
