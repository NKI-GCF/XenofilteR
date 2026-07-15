//! Chimeric (cross-species) fragment detection and routing.
//!
//! # Two detection phases
//!
//! ## Phase 1 — Mate-split (paired-end only)
//!
//! Different mates map to different streams with disjoint segment identifiers:
//!
//! ```text
//!   Human stream  : read1 mapped  (0x40)
//!   HPV stream    : read2 mapped  (0x80)
//!   → disjoint  →  chimeric
//! ```
//!
//! ## Phase 2 — Read-split (single-end or paired-end)
//!
//! The **same** mate appears as a primary alignment in **both** streams, but
//! with complementary soft-clip patterns that together cover the full read:
//!
//! ```text
//!   Human stream  : read1 primary  25M25S   mapped read range [0,  25)
//!   HPV stream    : read1 primary  25S25M   mapped read range [25, 50)
//!   → non-overlapping union ≈ full read  →  split-read chimeric
//! ```
//!
//! # False-positive rejection
//!
//! Aligners produce a **supplementary** alignment of the HPV portion of read1
//! on the human reference.  This supplementary is a false positive: HPV sequence
//! aligns poorly to human and its MD mismatch count will be high.  We compare:
//!
//! - `mismatches(stream_A_supplementary)` — HPV bases on human reference
//! - `mismatches(stream_B_primary)`        — HPV bases on HPV reference
//!
//! If the supplementary is *better* (fewer mismatches) than stream B's primary,
//! the split is likely a repetitive-sequence artefact and we do not call chimeric.
//!
//! # Three-stream example
//!
//! Human (0) + HPV (1) + Mouse (2), `--chimeric-pairs 0:1`:
//!
//! ```text
//! Fragment with HPV-integration breakpoint:
//!   read1 — 5′ 25 bp → human primary  +  HPV supplementary (false positive)
//!         — 3′ 25 bp → HPV primary    (after read-split detection)
//!   read2 — entirely → HPV primary
//!
//! Outcome:
//!   human output  : read1 primary  tagged XC:Z:hpv
//!                   read1 supp     tagged XC:Z:hpv  (false positive; flag 0x800)
//!   hpv output    : read1 primary  tagged XC:Z:human
//!                   read2 primary  tagged XC:Z:human
//!   mouse output  : all records    → discarded output (normal tournament)
//! ```

use super::core::{LineByLine, COUNTER_STRIDE};
use crate::alignment::{FragmentState, SimpleRec};
use crate::filter_algorithm::line_by_line::core::FragmentBuffer;
use crate::Error;
use noodles::sam::alignment::record::{
    cigar::op::Kind,
    data::field::{Tag, Value},
    Flags,
};
use noodles::sam::alignment::record_buf::data::field::Value as RecordBufValue;
use smallvec::SmallVec;

/// Fast flag-only check for paired-end chimeric complement.
///
/// Returns `Some((stream_a, stream_b))` when, for a configured chimeric pair,
/// one stream has mate0=unmapped + mate1=mapped and the other has the inverse.
/// Runs before `detect_chimeric_event` (which also scans CIGAR/MD for read-splits);
/// this path uses only flag checks.
pub(crate) fn detect_chimeric_mate_complement<R: SimpleRec>(
    best: &FragmentBuffer<R>,
    chimeric_pairs: &[[usize; 2]],
) -> Option<(usize, usize)> {
    use crate::alignment::mate_kind::{mate_slot, segment_id};

    /// Per-slot mapping: Some(true)=mapped-primary-present, Some(false)=unmapped-primary-present.
    fn flag_mate_map<R: SimpleRec>(state: &FragmentState<R>) -> [Option<bool>; 2] {
        let mut m = [None::<bool>; 2];
        for r in state.get_records() {
            let Ok(f) = r.flags() else { continue };
            if f.is_secondary() || f.is_supplementary() {
                continue;
            }
            let slot = mate_slot(segment_id(&f));
            // Later records for the same slot overwrite; last primary wins.
            m[slot] = Some(!f.is_unmapped());
        }
        m
    }

    for &[a, b] in chimeric_pairs {
        let sa = best.iter().find(|s| s.get_nr() == a)?;
        let sb = best.iter().find(|s| s.get_nr() == b)?;

        // Both streams must have paired-end records.
        let paired = sa
            .get_records()
            .iter()
            .chain(sb.get_records().iter())
            .any(|r| r.flags().is_ok_and(|f| f.is_segmented()));
        if !paired {
            continue;
        }

        let mk_a = flag_mate_map(sa);
        let mk_b = flag_mate_map(sb);

        // Complement: [unmapped, mapped] ↔ [mapped, unmapped]
        let complement = |x: [Option<bool>; 2], y: [Option<bool>; 2]| {
            x == [Some(false), Some(true)] && y == [Some(true), Some(false)]
        };
        if complement(mk_a, mk_b) || complement(mk_b, mk_a) {
            return Some((a, b));
        }
    }
    None
}

// ---------------------------------------------------------------------------
// ChimericKind
// ---------------------------------------------------------------------------

/// How the chimeric event was detected.
#[derive(Debug, Clone, Copy)]
pub(crate) enum ChimericKind {
    /// Different mates map to different streams; segment identifiers are disjoint.
    MateSplit,

    /// The same read spans a breakpoint: its aligned positions are split across
    /// two streams with complementary soft-clip patterns.
    ReadSplit {
        /// SAM flags segment identifier of the breakpoint-spanning read:
        /// `0x40` = read 1, `0x80` = read 2, `0x00` = single-end.
        read_seg_id: u8,
    },
}

/// Outcome of `detect_chimeric_event`.
#[derive(Debug)]
pub(crate) enum ChimericDecision {
    /// A chimeric event was detected across the two nominated streams.
    Chimeric {
        stream_a: usize,
        stream_b: usize,
        kind: ChimericKind,
    },
    /// No chimeric event; proceed with the normal tournament cascade.
    Normal,
}

// ---------------------------------------------------------------------------
// Segment-identity helpers
// ---------------------------------------------------------------------------

/// Encode the segment role of a record as a single byte:
/// `0x40` = first segment (read 1), `0x80` = last segment (read 2),
/// `0x00` = single-end or unknown.
fn segment_id(flags: &Flags) -> u8 {
    match (flags.is_first_segment(), flags.is_last_segment()) {
        (true, _) => 0x40,
        (false, true) => 0x80,
        (false, false) => 0x00,
    }
}

/// Collect segment identifiers for which a **primary, mapped** record exists.
fn mapped_segment_ids<R: SimpleRec>(state: &FragmentState<R>) -> SmallVec<[u8; 2]> {
    state
        .get_records()
        .iter()
        .filter_map(|r| {
            let f = r.flags().ok()?;
            if f.is_unmapped() || f.is_secondary() || f.is_supplementary() {
                return None;
            }
            Some(segment_id(&f))
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Read-coordinate geometry
// ---------------------------------------------------------------------------

/// Return the half-open range `[start, end)` of read positions that are
/// **aligned** (not soft- or hard-clipped) in the record's primary CIGAR.
///
/// Returns `(0, 0)` for unmapped records or on CIGAR parse failure.
fn mapped_read_range<R: SimpleRec>(rec: &R) -> (usize, usize) {
    if rec.flags().map_or(true, |f| f.is_unmapped()) {
        return (0, 0);
    }

    let mut read_pos = 0usize;
    let mut map_start: Option<usize> = None;
    let mut map_end = 0usize;

    for op_result in rec.cigar().as_ref().iter() {
        let op = match op_result {
            Ok(o) => o,
            Err(_) => break,
        };
        let len = op.len();
        match op.kind() {
            Kind::SoftClip => {
                read_pos += len;
            }
            Kind::HardClip | Kind::Pad => {}
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Insertion => {
                if map_start.is_none() {
                    map_start = Some(read_pos);
                }
                read_pos += len;
                map_end = read_pos;
            }
            Kind::Deletion | Kind::Skip => {}
        }
    }
    (map_start.unwrap_or(0), map_end)
}

/// True when `range_a` and `range_b` are largely non-overlapping and
/// together cover most of `read_len` bases.
///
/// Thresholds:
/// - Each arm must contribute ≥ 15 aligned bases.
/// - Overlap < 20 % of `read_len` (tolerance for aligner boundary fuzz).
/// - Union ≥ 80 % of `read_len` (together they explain the full read).
fn is_complementary(range_a: (usize, usize), range_b: (usize, usize), read_len: usize) -> bool {
    let (a0, a1) = range_a;
    let (b0, b1) = range_b;

    if a1 <= a0 || b1 <= b0 || read_len == 0 {
        return false;
    }
    let mapped_a = a1 - a0;
    let mapped_b = b1 - b0;
    if mapped_a < 15 || mapped_b < 15 {
        return false;
    }

    // Overlap between mapped ranges.
    let overlap = a1.min(b1).saturating_sub(a0.max(b0));
    // Union of mapped ranges.
    let union = a1.max(b1) - a0.min(b0);

    // overlap < 20 % of read_len  AND  union ≥ 80 % of read_len
    overlap * 5 < read_len && union * 10 >= read_len * 8
}

// ---------------------------------------------------------------------------
// Mismatch counting from MD tag
// ---------------------------------------------------------------------------

/// Count mismatch bases encoded in a raw MD tag byte string.
///
/// Digit runs (matches) and `^`-prefixed deletion blocks are skipped;
/// SNV letters (A/C/G/T/N) are counted as mismatches.  Returns 0 on
/// empty or malformed input and never panics.
fn md_mismatches_from_record<R: SimpleRec>(rec: &R) -> usize {
    match rec
        .data()
        .get(&Tag::MISMATCHED_POSITIONS)
        .and_then(|v| v.ok())
    {
        Some(Value::String(s)) => crate::alignment::pre_assess::md_mismatches(s.as_ref()),
        _ => usize::MAX, // MD absent → treat as maximally mismatched
    }
}

// ---------------------------------------------------------------------------
// Phase 2 — read-split detection
// ---------------------------------------------------------------------------

/// Attempt to detect a read-split chimeric event between `state_a` and `state_b`.
///
/// For each possible mate (read 1, read 2, single-end), checks whether:
///
/// 1. Both streams have a **primary** alignment for that mate.
/// 2. The two primaries' aligned read-position ranges are **complementary**.
/// 3. If stream A has a **supplementary** alignment for the same mate,
///    its mismatch count must be ≥ stream B's primary mismatch count
///    (otherwise stream A's supplementary is a better explanation for
///    the complementary region and the split is likely a false positive).
///
/// Returns `Some(seg_id)` on success, `None` if no read-split is found.
fn detect_split_read<R: SimpleRec>(
    state_a: &FragmentState<R>,
    state_b: &FragmentState<R>,
) -> Option<u8> {
    // Candidate segment identifiers (read1, read2, single-end).
    for &seg_id in &[0x40u8, 0x80u8, 0x00u8] {
        // Primary in stream A for this segment.
        let primary_a = state_a.get_records().iter().find(|r| {
            let Ok(f) = r.flags() else { return false };
            !f.is_secondary()
                && !f.is_supplementary()
                && !f.is_unmapped()
                && segment_id(&f) == seg_id
        });
        // Primary in stream B for this segment.
        let primary_b = state_b.get_records().iter().find(|r| {
            let Ok(f) = r.flags() else { return false };
            !f.is_secondary()
                && !f.is_supplementary()
                && !f.is_unmapped()
                && segment_id(&f) == seg_id
        });

        let (rec_a, rec_b) = match (primary_a, primary_b) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };

        // Read length: use whichever record has quality scores available.
        let read_len = rec_a
            .quality_scores()
            .as_ref()
            .len()
            .max(rec_b.quality_scores().as_ref().len());
        if read_len == 0 {
            continue;
        }

        let range_a = mapped_read_range(rec_a);
        let range_b = mapped_read_range(rec_b);

        if !is_complementary(range_a, range_b, read_len) {
            continue;
        }

        // ── False-positive rejection ──────────────────────────────────────
        // Stream A's supplementary for this segment maps the *same read region*
        // as stream B's primary but on the wrong reference genome.  If the
        // supplementary has *fewer* mismatches than stream B's primary, the
        // alignment in the other stream is not actually an improvement and the
        // complementary clip pattern is likely a repetitive-sequence artefact.
        let supp_a = state_a.get_records().iter().find(|r| {
            let Ok(f) = r.flags() else { return false };
            f.is_supplementary() && segment_id(&f) == seg_id
        });

        if let Some(sa) = supp_a {
            let mis_supp = md_mismatches_from_record(sa);
            let mis_b = md_mismatches_from_record(rec_b);

            if mis_supp < mis_b {
                // Stream A's supplementary alignment of the complementary region
                // is better than stream B's primary → false positive; skip.
                tracing::debug!(
                    seg_id,
                    mismatches_supp_a = mis_supp,
                    mismatches_primary_b = mis_b,
                    "Read-split candidate rejected: stream A supplementary \
                     explains the complementary region better than stream B"
                );
                continue;
            }
        }
        // No supplementary, or supplementary is worse → read-split confirmed.
        return Some(seg_id);
    }
    None
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Inspect `best` for chimeric events against every configured pair.
///
/// Runs two detection phases in order:
///
/// 1. **Mate-split** (paired-end only): segment-ID disjointness.
/// 2. **Read-split** (any library): complementary soft-clip analysis with
///    false-positive rejection via supplementary mismatch comparison.
///
/// Returns the first matching `Chimeric` decision, or `Normal`.
pub(crate) fn detect_chimeric_event<R: SimpleRec>(
    best: &FragmentBuffer<R>,
    chimeric_pairs: &[[usize; 2]],
) -> ChimericDecision {
    for &[a, b] in chimeric_pairs {
        let state_a = best.iter().find(|s| s.get_nr() == a);
        let state_b = best.iter().find(|s| s.get_nr() == b);

        let (sa, sb) = match (state_a, state_b) {
            (Some(x), Some(y)) => (x, y),
            _ => continue,
        };

        // ── Phase 1: mate-split (requires paired-end data) ───────────────
        let has_paired = sa
            .get_records()
            .iter()
            .chain(sb.get_records().iter())
            .any(|r| r.flags().is_ok_and(|f| f.is_segmented()));

        if has_paired {
            let ids_a = mapped_segment_ids(sa);
            let ids_b = mapped_segment_ids(sb);

            if !ids_a.is_empty() && !ids_b.is_empty() {
                let disjoint = ids_a.iter().all(|id| !ids_b.contains(id));
                if disjoint {
                    return ChimericDecision::Chimeric {
                        stream_a: a,
                        stream_b: b,
                        kind: ChimericKind::MateSplit,
                    };
                }
            }
        }

        // ── Phase 2: read-split (single-end or paired-end) ───────────────
        if let Some(seg_id) = detect_split_read(sa, sb) {
            return ChimericDecision::Chimeric {
                stream_a: a,
                stream_b: b,
                kind: ChimericKind::ReadSplit {
                    read_seg_id: seg_id,
                },
            };
        }
    }
    ChimericDecision::Normal
}

impl<R: SimpleRec> LineByLine<R> {
    pub(super) fn chimeric_label(&self, stream: usize) -> String {
        self.stream_labels
            .get(stream)
            .cloned()
            .unwrap_or_else(|| format!("stream_{stream}"))
    }

    /// Write all records from a chimeric fragment.
    ///
    /// For streams in the chimeric pair: records go to assigned output with
    /// `XC:Z:<other_stream_label>` tag.  For every other stream present in
    /// `best`: records go to filtered (discarded) output.
    ///
    /// The `XC:Z:` tag value is the human-readable label of the *other* stream
    /// (configured via `--stream-labels`), making it easy to filter chimeric
    /// reads with `samtools view -d XC:hpv`.
    pub(super) fn emit_chimeric(
        &mut self,
        best: &mut super::core::FragmentBuffer<R>,
        decision: ChimericDecision,
    ) -> Result<(), Error> {
        let (chimeric_a, chimeric_b, kind) = match decision {
            ChimericDecision::Chimeric {
                stream_a,
                stream_b,
                kind,
            } => (stream_a, stream_b, kind),
            ChimericDecision::Normal => unreachable!("emit_chimeric on non-chimeric"),
        };
        let label_a = self.chimeric_label(chimeric_a);
        let label_b = self.chimeric_label(chimeric_b);
        let tag_xc = Tag::new(b'X', b'C');

        tracing::debug!(kind = ?kind, stream_a = chimeric_a, stream_b = chimeric_b, "Emitting chimeric fragment");

        best.drain(..)
            .try_for_each(|mut state| -> Result<(), Error> {
                let nr = state.get_nr();
                let is_chimeric = nr == chimeric_a || nr == chimeric_b;
                if is_chimeric {
                    let other_label = if nr == chimeric_a { &label_b } else { &label_a };
                    let xc_value = RecordBufValue::String(other_label.as_bytes().into());

                    // For ReadSplit chimerism, stream A may contain a supplementary
                    // alignment that is a false-positive mapping of the split read's
                    // complementary portion to the wrong reference.  It is written
                    // here with the XC:Z tag so that its provenance is clear;
                    // downstream tools can identify it via the supplementary flag
                    // (SAM 0x800 / `samtools view -F 2048`).
                    //
                    // A future `--chimeric-suppress-supplementary` option could
                    // silently drop these records.
                    state
                        .drain_records()
                        .try_for_each(|r| -> Result<(), Error> {
                            let header = self.aln[nr].header();
                            let mut rb = r.as_record_buf(header)?;
                            rb.data_mut().insert(tag_xc, xc_value.clone());
                            self.routing_counters[nr * COUNTER_STRIDE + 3] += 1;
                            self.aln[nr].write_record(rb, Some(true))
                        })
                } else {
                    // Streams outside the chimeric pair are discarded normally.
                    state.drain_records().try_for_each(|r| {
                        let header = self.aln[nr].header();
                        let rb = r.as_record_buf(header)?;
                        self.routing_counters[nr * COUNTER_STRIDE] += 1;
                        self.aln[nr].write_record(rb, Some(false))
                    })
                }
            })
    }
}
