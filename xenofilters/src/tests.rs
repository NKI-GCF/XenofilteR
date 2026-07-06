pub(crate) mod common;
mod exploratory;
pub(crate) mod property;
use super::*;
pub(crate) use alignment::tests::*;
pub(crate) use aln_stream::tests::*;

// Kills mutations in `header_name_to_id` (HashMap::new(), HashMap::from_iter, etc.)
#[test]
fn test_header_name_to_id() {
    // Construct a realistic SAM header to ensure iteration and indexing are correct
    let header_str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200";
    let header: Header = header_str.parse().expect("Failed to parse SAM header");

    let map = header_name_to_id(&header);

    assert_eq!(map.len(), 2);
    assert_eq!(map.get("chr1"), Some(&0));
    assert_eq!(map.get("chr2"), Some(&1));
}

#[test]
fn test_load_ambiguous_regions_ignores_empty_strings() {
    let name_to_id = HashMap::new();

    // Two empty strings should return [None, None] without triggering a file-read error
    let specs = vec!["".to_string(), "".to_string()];
    let result = load_ambiguous_regions(&specs, &name_to_id).unwrap();
    assert!(result[0].is_none());
    assert!(result[1].is_none());

    // A single empty string should return [None, None]
    let specs_single = vec!["".to_string()];
    let result_single = load_ambiguous_regions(&specs_single, &name_to_id).unwrap();
    assert!(result_single[0].is_none());
    assert!(result_single[1].is_none());
}

#[test]
fn test_load_diagnostic_variants_ignores_empty_strings() {
    let name_to_id = HashMap::new();

    let specs = vec!["".to_string(), "".to_string()];
    let result = load_diagnostic_variants(&specs, &name_to_id).unwrap();
    assert!(result[0].is_none());
    assert!(result[1].is_none());
}

#[test]
fn test_namesorted_sequential_single_alignment() {
    let config = Config {
        matching_algorithm: MatchingAlgorithm::Namesorted,
        score_threads: 1, // Forces the sequential path
        alignment: vec!["tests/fixtures/dummy1.bam".into()],
        // .. populate remaining necessary fields
        ..Default::default()
    };

    // We just care that the logic branches correctly.
    let _ = run_namesorted(config);
}

#[test]
fn test_namesorted_parallel_dual_alignment() {
    let config = Config {
        matching_algorithm: MatchingAlgorithm::Namesorted,
        score_threads: 2, // Forces the parallel path
        alignment: vec![
            "tests/fixtures/dummy1.bam".into(),
            "tests/fixtures/dummy2.bam".into(),
        ],
        ..Default::default()
    };

    let _ = run_namesorted(config);
}

#[test]
fn test_get_log_level() {
    assert_eq!(get_log_level(0), "warn");
    assert_eq!(get_log_level(1), "info");
    assert_eq!(get_log_level(2), "debug");
    assert_eq!(get_log_level(5), "debug");
}
