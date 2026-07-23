use crate::config::args::{IoArgs, RelatedArgs, ScoringArgs};
use crate::file_spec::FileSpec;
use crate::filter_algorithm::strain::StrainArgs;

#[test]
fn test_scoring_validate_ok() {
    let mut s = ScoringArgs {
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Default::default()
    };
    assert!(s.validate().is_ok());
}

#[test]
fn test_scoring_validate_rejects_nonpositive_mismatch_penalty() {
    let mut s = ScoringArgs {
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 0.0,
        ..Default::default()
    };
    assert!(s.validate().is_err());
    s.mismatch_penalty = -1.0;
    assert!(s.validate().is_err());
}

#[test]
fn test_scoring_validate_rejects_zero_gap_open() {
    let mut s = ScoringArgs {
        gap_open: 0.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Default::default()
    };
    assert!(s.validate().is_err());
}

#[test]
fn test_scoring_validate_rejects_invalid_warn_ambig_fraction() {
    let mut s = ScoringArgs {
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        warn_ambig_fraction: 1.5,
        ..Default::default()
    };
    assert!(s.validate().is_err());
}

#[test]
fn test_to_penalty_q0_is_certain_error() {
    let s = ScoringArgs {
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Default::default()
    };
    let p = s.to_penalty();
    assert_eq!(p.gap_open, 6.0);
    assert!(p.log_likelihood_match[0].is_infinite());
}

#[test]
fn test_to_penalty_mismatch_scaling_factor() {
    let mut s = ScoringArgs {
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 8.0,
        ..Default::default()
    };
    assert!((s.to_penalty().log_likelihood_mismatch[10] - (-2.0)).abs() < 1e-9);
    s.mismatch_penalty = 4.0;
    assert!((s.to_penalty().log_likelihood_mismatch[10] - (-1.0)).abs() < 1e-9);
}

#[test]
fn test_to_penalty_match_likelihood_improves_with_quality() {
    let s = ScoringArgs {
        gap_open: 6.0,
        gap_extend: 1.0,
        mismatch_penalty: 4.0,
        ..Default::default()
    };
    let p = s.to_penalty();
    assert!(p.log_likelihood_match[40] > p.log_likelihood_match[10]);
}

#[test]
fn test_io_validate_ok_minimal() {
    let io = IoArgs {
        alignment: vec!["a.bam".into(), "b.bam".into()],
        output: vec!["out1.bam".into(), "out2.bam".into()],
        ..Default::default()
    };
    assert!(io.validate(1, 2..=2).is_ok());
}

#[test]
fn test_io_validate_rejects_too_many_output_paths() {
    let io = IoArgs {
        alignment: vec!["a.bam".into(), "b.bam".into()],
        output: vec!["o1".into(), "o2".into(), "o3".into()],
        ..Default::default()
    };
    assert!(io.validate(1, 2..=2).is_err());
}

#[test]
fn test_io_validate_rejects_invalid_stream_count() {
    let io = IoArgs {
        alignment: vec!["a.bam".into()],
        ..Default::default()
    };
    assert!(io.validate(1, 2..=2).is_err());
}

#[test]
fn test_strain_stream_count() {
    let mut args = StrainArgs {
        io: IoArgs {
            alignment: vec!["a.bam".into(), "b.bam".into()],
            ..Default::default()
        },
        scoring: ScoringArgs::default(),
        variants: RelatedArgs::default(),
        parallel: Default::default(),
    };
    assert!(args.validate_and_init().is_err());
}

#[test]
fn test_strain_requires_variant_profiles_for_both_strains() {
    let mut args = StrainArgs {
        io: IoArgs {
            alignment: vec!["single_strain.bam".into()],
            ..Default::default()
        },
        scoring: ScoringArgs {
            gap_open: 6.0,
            gap_extend: 1.0,
            mismatch_penalty: 4.0,
            ..Default::default()
        },
        variants: RelatedArgs::default(),
        parallel: Default::default(),
    };
    assert!(args.validate_and_init().is_err());

    args.variants.sample_variants = vec!["0:a.vcf".parse::<FileSpec>().unwrap()];
    assert!(args.validate_and_init().is_err());

    args.variants.population_variants = vec!["1:b.vcf".parse::<FileSpec>().unwrap()];
    assert!(args.validate_and_init().is_ok());
}
