//! Tier 2.5 pre-assessment: unified structural + read-space alignment comparison.
//!
//! # Design
//!
//! The pre-assessment runs two sub-tiers from a SINGLE CIGAR+MD walk per record:
//!
//! **Tier 2.5a — match-count domination.**
//! `AlignSig` stores the count of exact-match bases (CIGAR M/=/X positions where
//! MD shows a digit, i.e. `BaseOp::Match` from `ScoreOpIter`) and the number of
//! pending supplementary alignments (from the `SA:Z` tag on each primary).
//! A higher match count is unambiguously better; fewer supplementaries is better.
//! A stream dominates when it is ≥ on matches AND ≤ on supplementaries.
//! This single-axis match count replaces the old 3-axis (mismatches / clips / indels)
//! approach: it is more efficient (one number instead of three comparisons) and more
//! directly correlated with the NW log-likelihood score.
//!
//! **Tier 2.5b — read-coordinate-space comparison.**
//! When Tier 2.5a cannot decide, `compare_fragment_profiles` identifies per-position
//! dominance in read space (positions where one stream matches and the other does not).
//!
//! Both tiers share the same `FragmentProfile` built by `build_read_profile`; no
//! second CIGAR+MD walk is ever required.

use super::read_profile::{
    build_read_profile, compare_fragment_profiles, FragmentProfile, ReadOp, ReadSpaceDecision,
};
use crate::alignment::MdCigFlags;
use crate::filter_algorithm::line_by_line::READ_CT;
use smallvec::SmallVec;
use std::cmp::Ordering;

// ---------------------------------------------------------------------------
// AlignSig — match-count-based quality signature (private to this module)
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
/// the chimeric-junction penalty (gap_open + bases × gap_extend) is applied
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

/// Two-axis structural domination:
/// - more `primary_match_bases` is better (higher NW log-likelihood contribution)
/// - fewer `supp_count` is better (fewer chimeric-junction penalties)
///
/// Returns `None` when one axis favours A and the other favours B; the NW score
/// must then break the tie.
fn subsumes(a: &AlignSig, b: &AlignSig) -> Option<Ordering> {
    if a.supp_count == b.supp_count {
        if a.primary_match_bases > b.primary_match_bases {
            Some(Ordering::Greater)
        } else if a.primary_match_bases < b.primary_match_bases {
            Some(Ordering::Less)
        } else {
            Some(Ordering::Equal) // identical on both axes
        }
    } else if a.supp_count < b.supp_count {
        if a.primary_match_bases >= b.primary_match_bases {
            Some(Ordering::Greater)
        } else {
            None // a has fewer supps but fewer matches → incomparable
        }
    } else {
        if b.primary_match_bases >= a.primary_match_bases {
            Some(Ordering::Less)
        } else {
            None // b has fewer supps but fewer matches → incomparable
        }
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
            // M (0), = (7), X (8): aligned bases — classify each against MD
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
                            _ => return matches, // malformed — stop here
                        }
                    }
                }
            }
            // D (2): deletion from reference — skip the ^letters MD block
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
/// * **Tier 2.5a** — match-count domination via `AlignSig` (O(record) counter
///   increments derived for free from the already-built profile; no extra
///   CIGAR/MD parsing).
/// * **Tier 2.5b** — per-position read-space comparison via
///   `compare_fragment_profiles`.
///
/// Falls back to `FullScoring` on segment-count mismatch, malformed MD/CIGAR,
/// insertions, or when one stream has more matches but also more supplementaries.
pub fn pre_assess_alignments(
    mcfs_a: &SmallVec<[MdCigFlags<'_>; READ_CT]>,
    mcfs_b: &SmallVec<[MdCigFlags<'_>; READ_CT]>,
) -> PreAssessResult {
    if mcfs_a.len() != mcfs_b.len() || mcfs_a.is_empty() {
        #[cfg(test)]
        eprintln!(
            "pre_assess_alignments: segment count mismatch ({} vs {})",
            mcfs_a.len(),
            mcfs_b.len()
        );
        return PreAssessResult::FullScoring;
    }

    // Single CIGAR+MD walk per record; both tiers read from these profiles.
    let fp_a = FragmentProfile {
        mates: mcfs_a.iter().map(build_read_profile).collect(),
    };
    let fp_b = FragmentProfile {
        mates: mcfs_b.iter().map(build_read_profile).collect(),
    };

    if !fp_a.valid() || !fp_b.valid() {
        #[cfg(test)]
        eprintln!("pre_assess_alignments: invalid fragment profile (malformed CIGAR/MD)");
        return PreAssessResult::FullScoring;
    }

    // Tier 2.5a: match-count domination — derived from profiles at zero extra cost.
    let sig_a = sig_from_fragment_profile(&fp_a);
    let sig_b = sig_from_fragment_profile(&fp_b);
    if let Some(ord) = subsumes(&sig_a, &sig_b) {
        #[cfg(test)]
        eprintln!(
            "pre_assess_alignments: Tier 2.5a early decision1: {:?} vs {:?} => {:?}",
            sig_a, sig_b, ord
        );
        return PreAssessResult::EarlyDecision(ord);
    }

    // Tier 2.5b: per-position read-space comparison — reuses the same profiles.
    match compare_fragment_profiles(&fp_a, &fp_b) {
        ReadSpaceDecision::EarlyDecision(ord) => {
            #[cfg(test)]
            eprintln!(
                "pre_assess_alignments: Tier 2.5b early decision2: {:?} vs {:?} => {:?}",
                sig_a, sig_b, ord
            );
            PreAssessResult::EarlyDecision(ord)
        }
        ReadSpaceDecision::PartialScoring { .. } | ReadSpaceDecision::FallThrough => {
            #[cfg(test)]
            eprintln!(
                "pre_assess_alignments: Tier 2.5b fallthrough: {:?} vs {:?} => FullScoring",
                sig_a, sig_b
            );
            PreAssessResult::FullScoring
        }
    }
}

/// Tier 2.5 pre-assessment for `HashLookup`.
///
/// Operates on raw BAM-encoded CIGAR bytes + MD byte slices from `MappedRecord`,
/// avoiding `RecordBuf` construction until full NW is necessary. Uses
/// `match_count_raw` (same match-count semantics as Tier 2.5a above) and
/// respects `supp_count` from the SA:Z tag parsed in `next_scoring_record`.
pub(crate) fn pre_assess_scoring_records(
    recs_a: &[crate::filter_algorithm::hash_lookup::assemble::RecordKind],
    recs_b: &[crate::filter_algorithm::hash_lookup::assemble::RecordKind],
) -> PreAssessResult {
    use crate::filter_algorithm::hash_lookup::assemble::RecordKind;

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

    match subsumes(&sig_a, &sig_b) {
        Some(ord) => PreAssessResult::EarlyDecision(ord),
        None => PreAssessResult::FullScoring,
    }
}

#[cfg(test)]
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

    fn mcfs_from(cigar: &str, md: &str) -> SmallVec<[MdCigFlags<'static>; READ_CT]> {
        // Leak to get 'static lifetime for convenience in tests.
        let rec: &'static _ = Box::leak(Box::new(
            create_record(b"r", cigar, &[], &vec![30u8; 150], md, false).unwrap(),
        ));
        let flags: &'static _ = Box::leak(Box::new(rec.flags()));
        let mcf = MdCigFlags::try_from_record(rec, flags, false).unwrap();
        smallvec![mcf]
    }

    struct Row {
        label: &'static str,
        cigar_a: &'static str,
        md_a: &'static str,
        cigar_b: &'static str,
        md_b: &'static str,
        want: Option<Ordering>, // None → FullScoring
    }

    #[test]
    fn pre_assess_table() {
        // FIXME: missing supplementary reads and/or paired-end.
        let cases: &[Row] = &[
            // -- Aggregate subsumption (Tier 2.5a) ------------------------
            Row {
                label: "a perfect b imperfect → a wins",
                cigar_a: "10M",
                md_a: "10",
                cigar_b: "10M",
                md_b: "5A4",
                want: Some(Greater),
            },
            Row {
                label: "a imperfect b perfect → b wins",
                cigar_a: "10M",
                md_a: "5A4",
                cigar_b: "10M",
                md_b: "10",
                want: Some(Less),
            },
            Row {
                label: "both perfect → equal",
                cigar_a: "10M",
                md_a: "10",
                cigar_b: "10M",
                md_b: "10",
                want: Some(Equal),
            },
            Row {
                label: "a has soft clip b does not → b wins (more matches)",
                cigar_a: "5S5M",
                md_a: "5",
                cigar_b: "10M",
                md_b: "10",
                want: Some(Less),
            },
            Row {
                label: "a has 2 mismatch b has 1 → b wins",
                cigar_a: "10M",
                md_a: "3A3A2",
                cigar_b: "10M",
                md_b: "9A0",
                want: Some(Less),
            },
            Row {
                label: "incomparable: a wins on matches but b has more clips",
                // a: 8M2S → 8 matches; b: 10M, MD 8A1 → 9 matches
                // b wins (more matches)
                cigar_a: "8M2S",
                md_a: "8",
                cigar_b: "10M",
                md_b: "8A1",
                want: Some(Less),
            },
            Row {
                label: "a has deletion b does not → indels make a worse",
                // 5M+1D+5M = 10 matches; 10M = 10 matches, same match count
                // but deletion in a → supp_count same, match count same → falls to Tier 2.5b
                cigar_a: "5M1D5M",
                md_a: "5^A5",
                cigar_b: "10M",
                md_b: "10",
                want: Some(Less),
            },
            Row {
                label: "a has insertion b does not → fallthrough (read-space ambiguity)",
                cigar_a: "5M2I3M",
                md_a: "10",
                cigar_b: "10M",
                md_b: "10",
                want: Some(Less),
            },
        ];
        let mut misses = Vec::new();

        for c in cases {
            let mcfs_a = mcfs_from(c.cigar_a, c.md_a);
            let mcfs_b = mcfs_from(c.cigar_b, c.md_b);
            let result = pre_assess_alignments(&mcfs_a, &mcfs_b);
            match (c.want, result) {
                (Some(want_ord), PreAssessResult::EarlyDecision(got_ord)) => {
                    assert_eq!(got_ord, want_ord, "[{}]", c.label);
                }
                (None, PreAssessResult::FullScoring) => {}
                (want, got) => {
                    let got = match got {
                        PreAssessResult::EarlyDecision(o) => format!("Early({o:?})"),
                        PreAssessResult::FullScoring => "FullScoring".into(),
                    };
                    misses.push(format!("[{}] want {:?} got {:?}", c.label, want, got));
                }
            }
        }
        if !misses.is_empty() {
            panic!("{} test cases failed:\n{}", misses.len(), misses.join("\n"));
        }
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
                label: "equal matches, a no supp → a wins",
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
                label: "a more matches but more supp → incomparable → NW",
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
                label: "a more matches and equal supp → a wins",
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
                label: "identical → equal",
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
            assert_eq!(subsumes(&c.a, &c.b), c.want, "[{}]", c.label);
        }
    }
}
