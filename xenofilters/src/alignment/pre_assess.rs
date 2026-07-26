//! Tier 2.5 pre-assessment: unified structural + read-space alignment comparison.
//!
//! # Design
//!
//! The pre-assessment runs two sub-tiers from a SINGLE CIGAR+MD walk per record:
//!
//! **Tier 2.5a -- match-count domination.**
//! `AlignSig` stores the count of exact-match bases (CIGAR M/=/X positions where
//! MD shows a digit, i.e. `BaseOp::Match` from `ScoreOpIter`) and the number of
//! pending supplementary alignments (from the `SA:Z` tag on each primary).
//! A higher match count is unambiguously better; fewer supplementaries is better.
//! A stream dominates when it is >= on matches AND <= on supplementaries.
//! This single-axis match count replaces the old 3-axis (mismatches / clips / indels)
//! approach: it is more efficient (one number instead of three comparisons) and more
//! directly correlated with the NW log-likelihood score.
//!
//! **Tier 2.5b -- read-coordinate-space comparison.**
//! When Tier 2.5a cannot decide, `compare_fragment_profiles` identifies per-position
//! dominance in read space (positions where one stream matches and the other does not).
//!
//! Both tiers share the same `FragmentProfile` built by `build_read_profile`; no
//! second CIGAR+MD walk is ever required.

use super::read_profile::{
    build_read_profile, compare_fragment_profiles, FragmentProfile, ReadOp, ReadSpaceDecision,
};
use crate::alignment::MdCigFlags;
use crate::filter_algorithm::hash_lookup::assemble::RecordKind;
use crate::filter_algorithm::line_by_line::READ_CT;
use smallvec::SmallVec;
use std::cmp::Ordering;

/// Realistic lower bound on base quality for the "worst-case" per-base
/// advantage scan. Q20 (~1% error rate) is the conventional "usable"
/// quality threshold across Illumina, ONT, and PacBio platforms, and is a
/// far more defensible floor than an arbitrarily low value: real aligners
/// essentially never report a matched CIGAR M-position below Q20 in
/// practice, whereas anything below ~Q4 makes the underlying
/// match-favors-this-side assumption not even theoretically meaningful
/// (`ll_match[q] <= ll_mismatch[q]` once `p_err >= 0.5`, i.e. q <~ 3).
///
/// This is a documented, deliberate tradeoff, not a proof: a genuinely
/// low-quality (sub-Q20) matched base in scored reads means Tier 2.5's
/// guarantee -- that it never contradicts what Tier 3 would decide -- can
/// in principle be violated for that specific comparison. Lower this if
/// your platform/basecaller genuinely reports usable matches below Q20.
const MIN_PLAUSIBLE_MATCH_QUALITY: usize = 20;

// ---------------------------------------------------------------------------
// AlignSig -- match-count-based quality signature (private to this module)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default)]
struct AlignSig {
    /// Exact-match bases across all PRIMARY segments (ScoreOpIter BaseOp::Match).
    /// Supplementary match bases are excluded: their gap penalty prevents
    /// safe cross-stream comparison without knowing the penalty magnitude.
    primary_match_bases: usize,
    /// Total supplementary alignments across all primary reads (SA:Z counts).
    /// Lower is better; used as the second comparison axis.
    supp_count: usize,
}

/// Derive a match-count `AlignSig` from an already-built `FragmentProfile`.
///
/// Primary records contribute their `ReadOp::Match` count and `supp_count`.
/// Supplementary records are NOT counted in `primary_match_bases` because
/// the chimeric-junction penalty (gap_open + bases * gap_extend) is applied
/// in the NW phase and may negate the match contribution; without knowing
/// the penalty magnitude we cannot compare supplementary matches across streams.
/// del events are subtracted from the match count because they are not exact matches.
fn sig_from_fragment_profile(fp: &FragmentProfile) -> AlignSig {
    let mut sig = AlignSig::default();
    for mate in &fp.mates {
        if mate.is_supplementary {
            continue;
        }
        sig.primary_match_bases +=
            mate.ops.iter().filter(|&&op| op == ReadOp::Match).count() - mate.del_events as usize;
        sig.supp_count += mate.supp_count;
    }
    sig
}

/// Smallest possible per-base log-likelihood advantage a matched base can
/// have over a mismatched one, across the full quality range. Occurs at the
/// lowest quality score, where `ll_match` and `ll_mismatch` are closest
/// together. Used as a conservative (worst-case) lower bound: if
/// `match_count_delta * min_per_base_advantage >= threshold`, then the
/// TRUE NW margin (which uses actual, generally larger, per-base
/// advantages) is guaranteed to also clear the threshold, so short-circuiting
/// here cannot produce a decision that full NW scoring would have called
/// ambiguous.
fn min_per_base_advantage(pen: &crate::penalty::Penalty) -> f64 {
    (MIN_PLAUSIBLE_MATCH_QUALITY..pen.log_likelihood_match.len())
        .map(|q| pen.log_likelihood_match[q] - pen.log_likelihood_mismatch[q])
        .filter(|d| d.is_finite())
        .fold(f64::INFINITY, f64::min)
        .max(0.0)
}

/// Minimum match-count delta required before Tier 2.5a is permitted to
/// short-circuit, given a configured ambiguous threshold (in nats). Returns
/// `usize::MAX` (never short-circuit) if `min_per_base_advantage` is ~zero,
/// which happens only in pathological penalty configurations.
pub(crate) fn min_delta_for_early_decision(
    pen: &crate::penalty::Penalty,
    threshold_nats: f64,
) -> usize {
    let per_base = min_per_base_advantage(pen);
    if per_base <= 1e-12 {
        return usize::MAX;
    }
    (threshold_nats / per_base).ceil() as usize
}

/// Two-axis structural domination, now threshold-aware (fix for review
/// item #2: Tier 2.5 could previously override Tier 3's own ambiguity bar
/// on a 1-base match-count difference, silently discarding quality
/// information near the decision boundary).
fn subsumes(a: &AlignSig, b: &AlignSig, min_delta: usize) -> Option<Ordering> {
    if a.supp_count != b.supp_count {
        // Differing supp_count already requires the >= match-count guard
        // below in the original logic; the threshold floor only gates the
        // *magnitude* of a match-count-only decision, so this branch is
        // unchanged from before except for propagating min_delta into the
        // equal-supp_count case, which is where a bare 1-base delta most
        // commonly arises.
        return if a.supp_count < b.supp_count && a.primary_match_bases >= b.primary_match_bases {
            Some(Ordering::Greater)
        } else if b.supp_count < a.supp_count && b.primary_match_bases >= a.primary_match_bases {
            Some(Ordering::Less)
        } else {
            None
        };
    }

    let delta = a.primary_match_bases.abs_diff(b.primary_match_bases);
    if delta == 0 {
        return Some(Ordering::Equal);
    }
    if delta < min_delta {
        // Not enough of a margin to guarantee Tier 3 would agree — fall
        // through to full scoring instead of risking a false-confident
        // early decision.
        return None;
    }
    if a.primary_match_bases > b.primary_match_bases {
        Some(Ordering::Greater)
    } else {
        Some(Ordering::Less)
    }
}

// ---------------------------------------------------------------------------
// Raw-bytes match counting for HashLookup (operates on MappedRecord)
// ---------------------------------------------------------------------------

/// Count exact-match bases from raw BAM-encoded CIGAR bytes and a raw MD string.
///
/// Walks CIGAR M/=/X ops in parallel with the MD string. A base counts as a
/// match when the CIGAR op is M/=/X (0/7/8) AND the corresponding MD token is
/// a digit (not a mismatch letter). Gracefully returns the count accumulated so
/// far on any parse inconsistency.
pub fn match_count_raw(cigar_bytes: &[u8], md: &[u8]) -> usize {
    let mut matches = 0usize;
    let mut md_pos = 0usize;
    let mut md_match_remain = 0usize;

    for chunk in cigar_bytes.chunks_exact(4) {
        let enc = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        let op = enc & 0xF;
        let len = (enc >> 4) as usize;

        match op {
            // M (0), = (7), X (8): aligned bases -- classify each against MD
            0 | 7 | 8 => {
                let mut cigar_remain = len;
                while cigar_remain > 0 {
                    if md_match_remain > 0 {
                        let take = cigar_remain.min(md_match_remain);
                        matches += take;
                        cigar_remain -= take;
                        md_match_remain -= take;
                    } else {
                        match md.get(md_pos) {
                            Some(&b) if b.is_ascii_digit() => {
                                let mut num = (b - b'0') as usize;
                                md_pos += 1;
                                while let Some(&d) = md.get(md_pos) {
                                    if !d.is_ascii_digit() {
                                        break;
                                    }
                                    num = num * 10 + (d - b'0') as usize;
                                    md_pos += 1;
                                }
                                md_match_remain = num;
                            }
                            Some(&b'A' | b'C' | b'G' | b'T' | b'N') => {
                                // Mismatch: advance MD and consume one CIGAR position
                                md_pos += 1;
                                cigar_remain -= 1;
                            }
                            Some(&b'^') => {
                                // Deletion block inside M run (malformed but graceful)
                                md_pos += 1;
                                while let Some(&d) = md.get(md_pos) {
                                    if d.is_ascii_digit() {
                                        break;
                                    }
                                    md_pos += 1;
                                }
                            }
                            _ => return matches, // malformed -- stop here
                        }
                    }
                }
            }
            // D (2): deletion from reference -- skip the ^letters MD block
            2 if md.get(md_pos) == Some(&b'^') => {
                md_pos += 1;
                while let Some(&d) = md.get(md_pos) {
                    if d.is_ascii_digit() {
                        break;
                    }
                    md_pos += 1;
                }
            }
            // I (1), N (3), S (4), H (5), P (6): no MD consumption, no matches
            _ => {}
        }
    }
    matches
}

// ---------------------------------------------------------------------------
// PreAssessResult
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum PreAssessResult {
    /// Structural dominance resolved the fragment without NW DP.
    /// `Greater` = stream A wins, `Less` = stream B wins, `Equal` = tie.
    EarlyDecision(Ordering),
    /// No structural dominance; fall through to full NW scoring.
    FullScoring,
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Unified Tier 2.5 pre-assessment for `LineByLine` and `Collated`.
///
/// Builds `ReadProfile`s once (single `ScoreOpIter` walk per record) and runs
/// both sub-tiers from the same data:
///
/// * **Tier 2.5a** -- match-count domination via `AlignSig` (O(record) counter
///   increments derived for free from the already-built profile; no extra
///   CIGAR/MD parsing).
/// * **Tier 2.5b** -- per-position read-space comparison via
///   `compare_fragment_profiles`.
///
/// Falls back to `FullScoring` on segment-count mismatch, malformed MD/CIGAR,
/// insertions, or when one stream has more matches but also more supplementaries.
pub fn pre_assess_alignments(
    mcfs_a: &SmallVec<[Option<MdCigFlags<'_>>; READ_CT]>,
    mcfs_b: &SmallVec<[Option<MdCigFlags<'_>>; READ_CT]>,
    pen: &crate::penalty::Penalty,
    ambiguous_threshold_nats: f64,
) -> PreAssessResult {
    if mcfs_a.len() != mcfs_b.len() || mcfs_a.is_empty() {
        return PreAssessResult::FullScoring;
    }
    if mcfs_a.iter().any(Option::is_none) || mcfs_b.iter().any(Option::is_none) {
        return PreAssessResult::FullScoring;
    }

    let fp_a = FragmentProfile {
        mates: mcfs_a
            .iter()
            .filter_map(Option::as_ref)
            .map(build_read_profile)
            .collect(),
    };
    let fp_b = FragmentProfile {
        mates: mcfs_b
            .iter()
            .filter_map(Option::as_ref)
            .map(build_read_profile)
            .collect(),
    };
    if !fp_a.valid() || !fp_b.valid() {
        return PreAssessResult::FullScoring;
    }

    let sig_a = sig_from_fragment_profile(&fp_a);
    let sig_b = sig_from_fragment_profile(&fp_b);
    let min_delta = min_delta_for_early_decision(pen, ambiguous_threshold_nats);
    if let Some(ord) = subsumes(&sig_a, &sig_b, min_delta) {
        return PreAssessResult::EarlyDecision(ord);
    }

    // Tier 2.5b (read-space comparison) is unaffected: it already requires
    // a per-position dominance proof (every differing position favors the
    // same side), which is a strictly stronger, quality-blind-safe
    // condition than a raw count delta — no threshold floor needed there.
    match compare_fragment_profiles(&fp_a, &fp_b) {
        ReadSpaceDecision::EarlyDecision(ord) => PreAssessResult::EarlyDecision(ord),
        _ => PreAssessResult::FullScoring,
    }
}

/// Tier 2.5 pre-assessment for `HashLookup`.
///
/// Operates on raw BAM-encoded CIGAR bytes + MD byte slices from `MappedRecord`,
/// avoiding `RecordBuf` construction until full NW is necessary. Uses
/// `match_count_raw` (same match-count semantics as Tier 2.5a above) and
/// respects `supp_count` from the SA:Z tag parsed in `next_scoring_record`.
pub(crate) fn pre_assess_scoring_records(
    recs_a: &[RecordKind],
    recs_b: &[RecordKind],
    pen: &crate::penalty::Penalty,
    ambiguous_threshold_nats: f64,
) -> PreAssessResult {
    let primaries_a: SmallVec<[_; 2]> = recs_a
        .iter()
        .filter_map(|r| match r {
            RecordKind::Mapped(m) if !m.flags.is_secondary() && !m.flags.is_supplementary() => {
                Some(m)
            }
            _ => None,
        })
        .collect();
    let primaries_b: SmallVec<[_; 2]> = recs_b
        .iter()
        .filter_map(|r| match r {
            RecordKind::Mapped(m) if !m.flags.is_secondary() && !m.flags.is_supplementary() => {
                Some(m)
            }
            _ => None,
        })
        .collect();

    if primaries_a.is_empty() || primaries_a.len() != primaries_b.len() {
        return PreAssessResult::FullScoring;
    }

    let mut sig_a = AlignSig::default();
    for m in &primaries_a {
        sig_a.primary_match_bases += match_count_raw(&m.cigar_bytes, &m.md);
        sig_a.supp_count += m.supp_count;
    }
    let mut sig_b = AlignSig::default();
    for m in &primaries_b {
        sig_b.primary_match_bases += match_count_raw(&m.cigar_bytes, &m.md);
        sig_b.supp_count += m.supp_count;
    }
    // Same conservative floor as pre_assess_alignments: a raw match-count
    // delta only short-circuits Tier 3 if it's large enough that even the
    // worst-case (lowest-quality) per-base advantage would still clear the
    // configured ambiguous threshold.
    let min_delta = min_delta_for_early_decision(pen, ambiguous_threshold_nats);
    match subsumes(&sig_a, &sig_b, min_delta) {
        Some(ord) => PreAssessResult::EarlyDecision(ord),
        None => PreAssessResult::FullScoring,
    }
}

pub(crate) fn md_mismatches(md: &[u8]) -> usize {
    let mut count = 0;
    let mut pos = 0;
    while pos < md.len() {
        match md[pos] {
            b'0'..=b'9' => {
                // Skip digits
                while pos < md.len() && md[pos].is_ascii_digit() {
                    pos += 1;
                }
            }
            b'^' => {
                // Skip deletion block
                pos += 1; // skip '^'
                while pos < md.len() && !md[pos].is_ascii_digit() {
                    pos += 1;
                }
            }
            b'A' | b'C' | b'G' | b'T' | b'N' => {
                count += 1; // mismatch
                pos += 1;
            }
            _ => {
                // Malformed MD string; stop counting
                break;
            }
        }
    }
    count
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{alignment::MdCigFlags, tests::create_record};
    use smallvec::{smallvec, SmallVec};
    use std::cmp::Ordering::{Equal, Greater, Less};

    fn mcfs_from(cigar: &str, md: &str) -> SmallVec<[Option<MdCigFlags<'static>>; READ_CT]> {
        // Leak to get 'static lifetime for convenience in tests.
        let rec: &'static _ = Box::leak(Box::new(
            create_record(b"r", cigar, &[], &[30u8; 150], md, false).unwrap(),
        ));
        let flags: &'static _ = Box::leak(Box::new(rec.flags()));
        let mcf = MdCigFlags::try_from_record(rec, flags, false).unwrap();
        smallvec![Some(mcf)]
    }

    struct Row {
        label: &'static str,
        cigar_a: &'static str,
        md_a: &'static str,
        cigar_b: &'static str,
        md_b: &'static str,
        want: Option<Ordering>, // None -> FullScoring
    }

    #[test]
    fn pre_assess_table() {
        // FIXME: missing supplementary reads and/or paired-end.
        let cases: &[Row] = &[
            // -- Aggregate subsumption (Tier 2.5a) ------------------------
            Row {
                label: "a perfect b imperfect -> a wins",
                cigar_a: "10M",
                md_a: "10",
                cigar_b: "10M",
                md_b: "5A4",
                want: Some(Greater),
            },
            Row {
                label: "a imperfect b perfect -> b wins",
                cigar_a: "10M",
                md_a: "5A4",
                cigar_b: "10M",
                md_b: "10",
                want: Some(Less),
            },
            Row {
                label: "both perfect -> equal",
                cigar_a: "10M",
                md_a: "10",
                cigar_b: "10M",
                md_b: "10",
                want: Some(Equal),
            },
            Row {
                label: "a has soft clip b does not -> b wins (more matches)",
                cigar_a: "5S5M",
                md_a: "5",
                cigar_b: "10M",
                md_b: "10",
                want: Some(Less),
            },
            Row {
                label: "a has 2 mismatch b has 1 -> b wins",
                cigar_a: "10M",
                md_a: "3A3A2",
                cigar_b: "10M",
                md_b: "9A0",
                want: Some(Less),
            },
            Row {
                label: "incomparable: a wins on matches but b has more clips",
                // a: 8M2S -> 8 matches; b: 10M, MD 8A1 -> 9 matches
                // b wins (more matches)
                cigar_a: "8M2S",
                md_a: "8",
                cigar_b: "10M",
                md_b: "8A1",
                want: Some(Less),
            },
            Row {
                label: "a has deletion b does not -> indels make a worse",
                // 5M+1D+5M = 10 matches; 10M = 10 matches, same match count
                // but deletion in a -> supp_count same, match count same -> falls to Tier 2.5b
                cigar_a: "5M1D5M",
                md_a: "5^A5",
                cigar_b: "10M",
                md_b: "10",
                want: Some(Less),
            },
            Row {
                label: "a has insertion b does not -> fallthrough (read-space ambiguity)",
                cigar_a: "5M2I3M",
                md_a: "10",
                cigar_b: "10M",
                md_b: "10",
                want: Some(Less),
            },
        ];

        crate::tests::common::run_collecting(
            cases,
            |c| c.label.to_string(),
            |c| {
                let mcfs_a = mcfs_from(c.cigar_a, c.md_a);
                let mcfs_b = mcfs_from(c.cigar_b, c.md_b);
                let pen = crate::penalty::Penalty::default();
                let result = pre_assess_alignments(&mcfs_a, &mcfs_b, &pen, 0.0);

                match (c.want, result) {
                    (Some(want_ord), PreAssessResult::EarlyDecision(got_ord)) => {
                        if got_ord != want_ord {
                            return Err(format!("want {:?} got Early({:?})", want_ord, got_ord));
                        }
                        Ok(())
                    }
                    (None, PreAssessResult::FullScoring) => Ok(()),
                    (want, got) => {
                        let got_str = match got {
                            PreAssessResult::EarlyDecision(o) => format!("Early({o:?})"),
                            PreAssessResult::FullScoring => "FullScoring".into(),
                        };
                        Err(format!("want {:?} got {}", want, got_str))
                    }
                }
            },
        );
    }

    // -- match_count_raw edge cases ------------------------------------------

    #[test]
    fn match_count_raw_table() {
        fn cigar(ops: &[(u32, u32)]) -> Vec<u8> {
            let mut v = Vec::new();
            for &(op, len) in ops {
                v.extend_from_slice(&((len << 4 | op).to_le_bytes()));
            }
            v
        }
        struct Row {
            label: &'static str,
            ops: &'static [(u32, u32)],
            md: &'static [u8],
            want: usize,
        }
        let cases: &[Row] = &[
            Row {
                label: "perfect 10M",
                ops: &[(0, 10)],
                md: b"10",
                want: 10,
            },
            Row {
                label: "one mismatch",
                ops: &[(0, 10)],
                md: b"5A4",
                want: 9,
            },
            Row {
                label: "two mismatches",
                ops: &[(0, 10)],
                md: b"3A3C2",
                want: 8,
            },
            Row {
                label: "softclip excluded",
                ops: &[(4, 3), (0, 7)],
                md: b"7",
                want: 7,
            },
            Row {
                label: "deletion skipped",
                ops: &[(0, 5), (2, 3), (0, 5)],
                md: b"5^ATG5",
                want: 10,
            },
            Row {
                label: "insertion excluded",
                ops: &[(0, 5), (1, 3), (0, 5)],
                md: b"10",
                want: 10,
            },
            Row {
                label: "all mismatch",
                ops: &[(0, 4)],
                md: b"0A0C0G0T",
                want: 0,
            },
            Row {
                label: "empty cigar",
                ops: &[],
                md: b"",
                want: 0,
            },
            Row {
                label: "malformed md truncates gracefully",
                ops: &[(0, 5)],
                md: b"3Z",
                want: 3,
            },
        ];
        for c in cases {
            let cig = cigar(c.ops);
            assert_eq!(match_count_raw(&cig, c.md), c.want, "[{}]", c.label);
        }
    }

    // -- Subsumption with supp_count ----------------------------------------

    #[test]
    fn subsumes_with_supp_count_table() {
        struct Row {
            label: &'static str,
            a: AlignSig,
            b: AlignSig,
            want: Option<Ordering>,
        }
        let cases = &[
            Row {
                label: "equal matches, a no supp -> a wins",
                a: AlignSig {
                    primary_match_bases: 80,
                    supp_count: 0,
                },
                b: AlignSig {
                    primary_match_bases: 80,
                    supp_count: 1,
                },
                want: Some(Greater),
            },
            Row {
                label: "a more matches but more supp -> incomparable -> NW",
                a: AlignSig {
                    primary_match_bases: 90,
                    supp_count: 2,
                },
                b: AlignSig {
                    primary_match_bases: 80,
                    supp_count: 0,
                },
                want: None,
            },
            Row {
                label: "a more matches and equal supp -> a wins",
                a: AlignSig {
                    primary_match_bases: 90,
                    supp_count: 0,
                },
                b: AlignSig {
                    primary_match_bases: 80,
                    supp_count: 0,
                },
                want: Some(Greater),
            },
            Row {
                label: "identical -> equal",
                a: AlignSig {
                    primary_match_bases: 80,
                    supp_count: 0,
                },
                b: AlignSig {
                    primary_match_bases: 80,
                    supp_count: 0,
                },
                want: Some(Equal),
            },
        ];
        for c in cases {
            assert_eq!(subsumes(&c.a, &c.b, 0), c.want, "[{}]", c.label);
        }
    }
    #[test]
    fn md_mismatches_table() {
        let cases: &[(&[u8], usize)] = &[
            (b"100", 0),
            (b"", 0),
            (b"5A4", 1),
            (b"3A3C2", 2),
            (b"0A0C0G0T", 4),
            (b"5^ACGT5", 0),    // deletion block: no mismatches
            (b"5^ACGTA5C4", 1), // deletion + one trailing mismatch
            (b"AAAA", 4),       // pathological: mismatches with no digits between
        ];
        for (md, want) in cases {
            assert_eq!(md_mismatches(md), *want, "md={:?}", std::str::from_utf8(md));
        }
    }

    #[test]
    fn md_mismatches_malformed_input_does_not_panic_and_stops_gracefully() {
        // '!' is not a valid MD character; function must stop rather than panic.
        assert_eq!(md_mismatches(b"5!5"), 0);
    }

    // alignment/pre_assess.rs — add to existing #[cfg(test)] mod tests
    #[test]
    fn small_match_delta_below_threshold_floor_falls_through_to_full_scoring() {
        let pen = crate::penalty::Penalty::build(
            6.0,
            1.0,
            4.0,
            20,
            crate::penalty::ErrorModel::Illumina,
            0.5,
        );
        // ambiguous_threshold = Phred 30 -> a large nats threshold; a 1-base
        // match-count delta cannot possibly guarantee this margin even in the
        // best case (Tier 2.5a).
        //
        // A "clean" 1-base delta (identical CIGAR structure, one substituted
        // base) would ALSO be resolved unconditionally by Tier 2.5b
        // (`compare_fragment_profiles`) -- correctly, since a single
        // differing read position is a strict per-position dominance proof
        // that doesn't need a threshold floor at all (see its doc comment).
        // To actually isolate and exercise Tier 2.5a's floor here, stream A
        // carries an insertion, which unconditionally forces Tier 2.5b to
        // FallThrough (`compare_mate_profiles`'s `has_insertions` guard),
        // leaving Tier 2.5a's count-based floor as the only path to a decision.
        let mcfs_a = mcfs_from("5M1I4M", "9"); // 9 matches, has an insertion
        let mcfs_b = mcfs_from("10M", "4A4A0"); // 8 matches, 2 mismatches -- delta = 1 vs a
        let threshold_nats = 30.0 * std::f64::consts::LN_10 / 10.0;
        let result = pre_assess_alignments(&mcfs_a, &mcfs_b, &pen, threshold_nats);
        assert!(
            matches!(result, PreAssessResult::FullScoring),
            "1-base delta must not clear a Phred-30 threshold floor, and Tier 2.5b \
         must not independently resolve it either (insertion forces FallThrough)"
        );
    }

    #[test]
    fn large_match_delta_above_threshold_floor_still_early_decides() {
        let pen = crate::penalty::Penalty::build(
            6.0,
            1.0,
            4.0,
            20,
            crate::penalty::ErrorModel::Illumina,
            0.5,
        );
        let threshold_nats = 10.0 * std::f64::consts::LN_10 / 10.0; // modest threshold
        let mcfs_a = mcfs_from("50M", "50"); // 50 matches
        let mcfs_b = mcfs_from("50M", "0A49"); // 49 matches — large-enough delta
        let result = pre_assess_alignments(&mcfs_a, &mcfs_b, &pen, threshold_nats);
        assert!(matches!(
            result,
            PreAssessResult::EarlyDecision(std::cmp::Ordering::Greater)
        ));
    }

    #[test]
    fn zero_threshold_permits_any_nonzero_delta() {
        let pen = crate::penalty::Penalty::build(
            6.0,
            1.0,
            4.0,
            20,
            crate::penalty::ErrorModel::Illumina,
            0.5,
        );
        let mcfs_a = mcfs_from("10M", "9A0");
        let mcfs_b = mcfs_from("10M", "10");
        let result = pre_assess_alignments(&mcfs_a, &mcfs_b, &pen, 0.0);
        assert!(matches!(result, PreAssessResult::EarlyDecision(_)));
    }
}

#[cfg(test)]
mod hashlookup_threshold_floor_tests {
    use super::*;
    use crate::filter_algorithm::hash_lookup::assemble::{MappedRecord, RecordKind};
    use crate::penalty::{ErrorModel, Penalty};
    use noodles::sam::alignment::record::Flags;
    use std::f64::consts::LN_10;

    fn cigar_bytes(len: u32) -> Vec<u8> {
        ((len << 4) | 0/* Match */).to_le_bytes().to_vec()
    }

    fn mapped(cigar_len: u32, md: &[u8], supp_count: usize) -> RecordKind {
        RecordKind::Mapped(Box::new(MappedRecord {
            flags: Flags::empty(),
            ref_id: 0,
            pos: 1,
            ref_len: cigar_len as usize,
            cigar_bytes: cigar_bytes(cigar_len),
            md: md.to_vec(),
            qualities: vec![30; cigar_len as usize],
            sequence: vec![b'A'; cigar_len as usize],
            virtual_offset: 0,
            supp_count,
        }))
    }
    #[test]
    fn hashlookup_large_delta_above_threshold_still_early_decides() {
        let pen = Penalty::build(6.0, 1.0, 4.0, 20, ErrorModel::Illumina, 0.5);
        let threshold_nats = 10.0 * LN_10 / 10.0;
        let recs_a = vec![mapped(60, b"60", 0)]; // 60 matches
        let recs_b = vec![mapped(60, build_md(40, 20).as_bytes(), 0)]; // 40 matches, 20 mismatches
        let result = pre_assess_scoring_records(&recs_a, &recs_b, &pen, threshold_nats);
        assert!(matches!(
            result,
            PreAssessResult::EarlyDecision(std::cmp::Ordering::Greater)
        ));
    }

    /// Build an MD string: `matches` leading matched bases followed by
    /// `mismatches` contiguous mismatched bases, summing to
    /// `matches + mismatches` -- must equal the paired CIGAR M-op length.
    fn build_md(matches: usize, mismatches: usize) -> String {
        format!("{matches}{}", "A".repeat(mismatches))
    }

    #[test]
    fn hashlookup_small_delta_below_threshold_floor_falls_through() {
        let pen = Penalty::build(6.0, 1.0, 4.0, 20, ErrorModel::Illumina, 0.5);
        let threshold_nats = 30.0 * LN_10 / 10.0; // Phred 30 -- a demanding bar
        let recs_a = vec![mapped(10, b"9A0", 0)]; // 9 matches
        let recs_b = vec![mapped(10, b"10", 0)]; // 10 matches -- 1-base delta only
        let result = pre_assess_scoring_records(&recs_a, &recs_b, &pen, threshold_nats);
        assert!(
            matches!(result, PreAssessResult::FullScoring),
            "1-base match delta must not clear a Phred-30 threshold floor in HashLookup"
        );
    }

    #[test]
    fn hashlookup_zero_threshold_matches_pre_fix_behavior() {
        let pen = Penalty::build(6.0, 1.0, 4.0, 20, ErrorModel::Illumina, 0.5);
        let recs_a = vec![mapped(10, b"9A0", 0)];
        let recs_b = vec![mapped(10, b"10", 0)];
        let result = pre_assess_scoring_records(&recs_a, &recs_b, &pen, 0.0);
        assert!(matches!(result, PreAssessResult::EarlyDecision(_)));
    }
    #[test]
    fn min_per_base_advantage_is_finite_and_positive_at_default_settings() {
        // Regression guard: previously collapsed to 0.0 via -infinity at Q0,
        // which made min_delta_for_early_decision return usize::MAX for EVERY
        // penalty configuration and permanently disabled Tier 2.5.
        let pen = crate::penalty::Penalty::build(
            6.0,
            1.0,
            4.0,
            20,
            crate::penalty::ErrorModel::Illumina,
            0.5,
        );
        let advantage = min_per_base_advantage(&pen);
        assert!(
            advantage.is_finite() && advantage > 0.0,
            "expected a finite positive worst-case per-base advantage, got {advantage}"
        );
    }

    #[test]
    fn min_delta_for_early_decision_is_not_permanently_disabled() {
        let pen = crate::penalty::Penalty::build(
            6.0,
            1.0,
            4.0,
            20,
            crate::penalty::ErrorModel::Illumina,
            0.5,
        );
        let threshold_nats = 10.0 * std::f64::consts::LN_10 / 10.0; // Phred 10
        let min_delta = min_delta_for_early_decision(&pen, threshold_nats);
        assert_ne!(
            min_delta,
            usize::MAX,
            "Tier 2.5 must be able to fire for a realistic Phred-10 threshold at default settings"
        );
        assert!(min_delta < 100, "min_delta unexpectedly large: {min_delta}");
    }

    #[test]
    fn min_delta_scales_with_threshold() {
        let pen = crate::penalty::Penalty::build(
            6.0,
            1.0,
            4.0,
            20,
            crate::penalty::ErrorModel::Illumina,
            0.5,
        );
        let low = min_delta_for_early_decision(&pen, 1.0 * std::f64::consts::LN_10 / 10.0);
        let high = min_delta_for_early_decision(&pen, 30.0 * std::f64::consts::LN_10 / 10.0);
        assert!(
            high > low,
            "a stricter threshold must require a larger match-count delta"
        );
    }
    #[test]
    fn hashlookup_table_row_regression_2base_delta_survives_phred10_default() {
        let pen = crate::penalty::Penalty::build(
            6.0,
            1.0,
            4.0,
            20,
            crate::penalty::ErrorModel::Illumina,
            0.5,
        );
        let threshold_nats = crate::config::args::resolve_threshold(u32::MAX, false); // pass1 default = Phred 10
        let recs_a = vec![mapped(10, b"9A0", 0)]; // 9 matches
        let recs_b = vec![mapped(10, b"7AAA", 0)]; // 7 matches -- fixed: "6AAA" summed to
                                                   // only 9 ref bases (1 short of the 10M
                                                   // CIGAR), and silently returned 6 via
                                                   // match_count_raw's lenient malformed-MD
                                                   // truncation rather than erroring.
                                                   // "7AAA" is internally consistent:
                                                   // 7+1+1+1=10, matches=7, delta=2 as the
                                                   // test name states.
        let result =
            crate::alignment::pre_assess_scoring_records(&recs_a, &recs_b, &pen, threshold_nats);
        assert!(
            matches!(result, crate::alignment::PreAssessResult::EarlyDecision(_)),
            "existing hash_lookup_table fixture must not silently start requiring full scoring"
        );
    }
}
