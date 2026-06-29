//! Chimeric (cross-species) fragment detection and routing.
//!
//! # What is a chimeric fragment?
//!
//! In paired-end sequencing of tissues harbouring viral integration (e.g. HPV),
//! a read pair that spans an integration breakpoint has one mate aligning to the
//! host genome and the other aligning to the viral genome.  Neither mate
//! maps well to the other species' reference, so the normal scoring cascade would
//! call the fragment ambiguous or discard one stream.
//!
//! Xenofilter detects this as a **chimeric fragment** when:
//!  - The two streams form a configured `--chimeric-pairs A:B` pair.
//!  - Stream A has at least one mapped primary segment.
//!  - Stream B has at least one mapped primary segment.
//!  - The *sets of mapped segment identifiers* (read 1 / read 2) are disjoint:
//!    no mate maps well in *both* streams simultaneously.
//!
//! # What happens to chimeric fragments?
//!
//! Both streams are treated as "winners" for their respective mapped mates:
//!  - Stream A's records are written to stream A's assigned output.
//!  - Stream B's records are written to stream B's assigned output.
//!  - An `XC:Z:<other_stream_label>` aux tag is added to every written record.
//!  - Streams not involved in the chimeric pair compete in the normal tournament
//!    and are discarded when they lose.
//!
//! # Single-end data
//!
//! Chimeric detection requires paired-end reads; `is_segmented()` must be true
//! for at least one record in each `FragmentState`.  Single-end reads that span
//! an integration breakpoint appear as supplementary alignments and are handled
//! by the supplementary structural penalty already applied in Tier 3.
//!
//! # Three-stream example (HPV + human tissue xenografted in mouse)
//!
//! ```text
//! --chimeric-pairs 0:1          human ↔ HPV  may be chimeric
//! stream 2 (mouse)              competes normally for non-chimeric fragments
//! ```
//!
//! If a fragment has read1 → human (0) and read2 → HPV (1):
//!   - Chimeric event detected for pair [0,1].
//!   - Stream 0 records → human output with `XC:Z:hpv`.
//!   - Stream 1 records → HPV output with `XC:Z:human`.
//!   - Stream 2 (mouse) records → mouse discarded output.

use crate::alignment::{FragmentState, SimpleRec};
use crate::filter_algorithm::line_by_line::core::FragmentBuffer;
use smallvec::SmallVec;

// ---------------------------------------------------------------------------
// Segment-identity helpers
// ---------------------------------------------------------------------------

/// Returns the set of *segment identifiers* for which a primary, mapped record
/// exists in `state`.
///
/// Segment identifier encoding:
///  - `0x40`  = first segment in template (read 1)
///  - `0x80`  = last  segment in template (read 2)
///  - `0x00`  = single-end (neither flag set)
///
/// Secondary and supplementary records are excluded; they do not represent
/// independent mate alignments.
fn mapped_segment_ids<R: SimpleRec>(state: &FragmentState<R>) -> SmallVec<[u8; 2]> {
    state
        .get_records()
        .iter()
        .filter_map(|r| {
            let f = r.flags().ok()?;
            if f.is_unmapped() || f.is_secondary() || f.is_supplementary() {
                return None;
            }
            // Encode which segment this primary record represents.
            let seg_id: u8 = match (f.is_first_segment(), f.is_last_segment()) {
                (true,  _)     => 0x40,
                (false, true)  => 0x80,
                (false, false) => 0x00,  // single-end
            };
            Some(seg_id)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// ChimericDecision
// ---------------------------------------------------------------------------

/// Outcome of `detect_chimeric_event`.
#[derive(Debug)]
pub(crate) enum ChimericDecision {
    /// The fragment spans a species boundary: mates split across these two streams.
    ///
    /// Both streams' records should be written as winners with an `XC:Z:` tag.
    /// Streams not in the pair compete normally and are discarded when they lose.
    Chimeric {
        /// First (lower-index) stream in the detected pair.
        stream_a: usize,
        /// Second (higher-index) stream in the detected pair.
        stream_b: usize,
    },
    /// No chimeric event detected; proceed with the normal tournament cascade.
    Normal,
}

/// Inspect `best` for cross-species chimeric events against the configured pairs.
///
/// Returns the first matching `ChimericDecision::Chimeric` found, or
/// `ChimericDecision::Normal` if none of the configured pairs match.
///
/// **Detection conditions (both must hold for a pair [A, B]):**
/// 1. Both stream A and stream B have at least one mapped primary segment.
/// 2. The sets of mapped segment identifiers are *disjoint*: no mate maps
///    as a primary alignment in both streams simultaneously.
///
/// Condition 2 ensures we do not classify genuinely ambiguous fragments
/// (where a mate aligns well to both species) as chimeric.
pub(crate) fn detect_chimeric_event<R: SimpleRec>(
    best: &FragmentBuffer<R>,
    chimeric_pairs: &[[usize; 2]],
) -> ChimericDecision {
    for &[a, b] in chimeric_pairs {
        let state_a = best.iter().find(|s| s.get_nr() == a);
        let state_b = best.iter().find(|s| s.get_nr() == b);

        let (sa, sb) = match (state_a, state_b) {
            (Some(x), Some(y)) => (x, y),
            _ => continue,  // one or both streams absent for this fragment
        };

        // Require paired-end data: at least one stream must be segmented.
        let is_paired = sa.get_records().iter().chain(sb.get_records().iter())
            .any(|r| r.flags().map_or(false, |f| f.is_segmented()));
        if !is_paired {
            continue;
        }

        let ids_a = mapped_segment_ids(sa);
        let ids_b = mapped_segment_ids(sb);

        // Both streams must contribute at least one mapped segment.
        if ids_a.is_empty() || ids_b.is_empty() {
            continue;
        }

        // Disjointness check: no segment identifier appears in both sets.
        let disjoint = ids_a.iter().all(|id| !ids_b.contains(id));
        if disjoint {
            return ChimericDecision::Chimeric { stream_a: a, stream_b: b };
        }
    }
    ChimericDecision::Normal
}
