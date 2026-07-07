//! Fragment-level decision logic and the parallel scoring pipeline.
//!
//! # Architecture
//!
//! Reading and writing are inherently sequential (name-sorted BAM streams must
//! advance in lockstep), but **scoring is embarrassingly parallel**: once a
//! fragment's records have been collected from every stream, scoring it is
//! independent of every other fragment.
//!
//! ```text
//!  IO thread (main)
//!  ------------------------------------------------------------------
//!  ingest_record() per stream → assemble FragmentBundle → work_tx ->.
//!                                                                   | N workers
//!  result_rx <-- ScoredFragment <---------------------------------- +
//!  |                                                                |
//!  +-> write output (writers stay on IO thread; no Mutex needed)    |
//!                                                                   |
//!  Worker threads  (one per --score-threads)                        |
//!  ---------------------------------------------                    |
//!  <----------------------------------------------------------------'
//!  own Scratch (DP tables)
//!  read Arc<Penalty> + Arc<dyn StoreTrait>  (shared, immutable)
//!  compute compute_score_delta + decision
//!  send ScoredFragment back
//! ```
//!
//! Writers stay on the IO thread so `BamOutput` needs no `Mutex`.
//! Output order is nondeterministic across fragments when
//! `score_threads > 1` (pass `--score-threads 1` to restore order).
//!
//! # No-variant fast path
//!
//! When no stream has a variant store (`stores` are all `None`), the inner
//! per-segment `overlapping_multi` calls are skipped entirely.  The flag is
//! computed once per bundle in `score_candidate_owned` before the mate loop,
//! so the common no-variant case pays only one `Option::is_some()` check per
//! candidate rather than one per read segment.

use super::core::{FragmentBuffer, LineByLine, Scratch};
use crate::alignment::{
    mate_slot, segment_id, BaseOp, FragmentState, MateClassifiable, MateKind, MdCigFlags,
    ScoreOpIter, SimpleRec,
};
use crate::filter_algorithm::line_by_line::{
    chimeric::{detect_chimeric_mate_complement, ChimericKind},
    detect_chimeric_event, ChimericDecision, MAX_STREAMS, READ_CT,
};
use crate::variant::StoreTrait;
use crate::Error;
use noodles::sam::alignment::RecordBuf;
use smallvec::{smallvec, SmallVec};
use std::f64::consts::LN_10;
use std::sync::Arc;
use crate::config::Config;

pub(crate) fn compute_cancel_slot<R: SimpleRec>(best: &FragmentBuffer<R>) -> [bool; 2] {
    // -- Pass 0: per-mate classification (cheap, borrow-only) -----------
    // Must precede scoring so cancel_slot is known before any NW call.
    let mut mate_kinds_arr = [[None, None]; MAX_STREAMS];
    for s in best.iter() {
        mate_kinds_arr[s.get_nr()] = s.mate_kinds();
    }
    let mut cancel_kind: [Option<MateKind>; 2] = [None, None];
    let mut cancel_ok = [true, true];
    for s in best.iter() {
        let mk = mate_kinds_arr[s.get_nr()];
        for slot in 0..2 {
            match (mk[slot], cancel_kind[slot]) {
                (Some(MateKind::Other), _) | (_, _)
                    if matches!(mk[slot], Some(MateKind::Other)) =>
                {
                    cancel_ok[slot] = false;
                }
                (Some(k), None) => cancel_kind[slot] = Some(k),
                (Some(k), Some(prev)) if prev != k => cancel_ok[slot] = false,
                _ => {}
            }
        }
    }
    [
        cancel_ok[0] && cancel_kind[0].is_some(),
        cancel_ok[1] && cancel_kind[1].is_some(),
    ]
}

/// Per-stream metrics produced by a single CIGAR+MD scan.
/// Indexed by `FragmentState::get_nr()`, not by position in `best`.
struct TournamentMetrics {
    is_perfect: [bool; MAX_STREAMS],
    match_bases: [usize; MAX_STREAMS],
    supp_counts: [usize; MAX_STREAMS],
    nw_scores: [f64; MAX_STREAMS],
    vdeltas: [f64; MAX_STREAMS],
    any_perfect: bool,
    any_imperfect: bool,
}

/// Build `TournamentMetrics` with one CIGAR+MD walk per stream.
///
/// `score_fn(nr, state, mcfs, cancel_slot) -> Result<(nw_score, variant_delta)>`
/// is called only for imperfect streams; for perfect ones it is not invoked.
fn compute_metrics<R: SimpleRec, F>(
    best: &FragmentBuffer<R>,
    cancel_slot: [bool; 2],
    mut score_fn: F,
) -> Result<TournamentMetrics, Error>
where
    F: FnMut(
        usize,
        &FragmentState<R>,
        SmallVec<[MdCigFlags<'_>; READ_CT]>,
        [bool; 2],
    ) -> Result<(f64, f64), Error>,
{
    let mut m = TournamentMetrics {
        is_perfect: [false; MAX_STREAMS],
        match_bases: [0; MAX_STREAMS],
        supp_counts: [0; MAX_STREAMS],
        nw_scores: [f64::NEG_INFINITY; MAX_STREAMS],
        vdeltas: [0.0; MAX_STREAMS],
        any_perfect: false,
        any_imperfect: false,
    };

    for i in 0..best.len() {
        let nr = best[i].get_nr();
        let mcfs = best[i].build_mcfs()?;

        let perf = mcfs
            .iter()
            .filter(|m| !m.is_supplementary())
            .all(|m| m.is_perfect());
        m.is_perfect[nr] = perf;
        if perf {
            m.any_perfect = true;
        } else {
            m.any_imperfect = true;
        }

        for mcf in &mcfs {
            if mcf.is_supplementary() {
                continue;
            }
            for op in ScoreOpIter::new(mcf) {
                if matches!(op, Ok(BaseOp::Match)) {
                    m.match_bases[nr] += 1;
                }
            }
            m.supp_counts[nr] += mcf.supp_count();
        }

        if !perf {
            let (nw, vd) = score_fn(nr, &best[i], mcfs, cancel_slot)?;
            m.nw_scores[nr] = nw;
            m.vdeltas[nr] = vd;
        }
    }
    Ok(m)
}

/// Apply tier 2 / 2.5 / 3 elimination sweeps.
///
/// `discard` is called for each loser index (backward sweep order guaranteed).
/// Returns early `Some(Decision)` when resolved; `None` means the discard fn
/// handled all elimination and callers should emit whatever remains in `best`.
fn apply_tiers<R: SimpleRec, D>(
    best: &mut FragmentBuffer<R>,
    m: &TournamentMetrics,
    threshold: f64,
    add_tag: bool,
    mut discard: D,
) -> Result<Option<Decision>, Error>
where
    D: FnMut(&mut FragmentBuffer<R>, usize) -> Result<(), Error>,
{
    // Tier 2
    if m.any_perfect && m.any_imperfect {
        for i in (0..best.len()).rev() {
            if !m.is_perfect[best[i].get_nr()] {
                discard(best, i)?;
            }
        }
        return Ok(add_tag.then_some(if best.len() == 1 {
            Decision::First
        } else {
            Decision::Ambiguous
        }));
    }
    if !m.any_imperfect {
        return Ok(add_tag.then_some(Decision::Ambiguous));
    }

    // Tier 2.5
    let max_mb = best
        .iter()
        .map(|s| m.match_bases[s.get_nr()])
        .max()
        .unwrap_or(0);
    let min_sc = best
        .iter()
        .filter(|s| m.match_bases[s.get_nr()] == max_mb)
        .map(|s| m.supp_counts[s.get_nr()])
        .min()
        .unwrap_or(0);
    if best.iter().any(|s| {
        let nr = s.get_nr();
        m.match_bases[nr] < max_mb || m.supp_counts[nr] > min_sc
    }) {
        for i in (0..best.len()).rev() {
            let nr = best[i].get_nr();
            if m.match_bases[nr] < max_mb || m.supp_counts[nr] > min_sc {
                discard(best, i)?;
            }
        }
        if best.len() == 1 {
            return Ok(add_tag.then_some(Decision::First));
        }
    }

    // Tier 3
    let max_nw = best
        .iter()
        .map(|s| m.nw_scores[s.get_nr()])
        .fold(f64::NEG_INFINITY, f64::max);
    for i in (0..best.len()).rev() {
        if m.nw_scores[best[i].get_nr()] < max_nw - threshold {
            discard(best, i)?;
        }
    }

    if !add_tag {
        return Ok(None);
    }
    if best.len() != 1 {
        return Ok(Some(Decision::Ambiguous));
    }

    let winner_nr = best[0].get_nr();
    let second = m
        .nw_scores
        .iter()
        .enumerate()
        .filter(|&(nr, &s)| nr != winner_nr && s > f64::NEG_INFINITY)
        .map(|(_, &s)| s)
        .fold(f64::NEG_INFINITY, f64::max);
    let margin = if second.is_finite() {
        max_nw - second
    } else {
        f64::MAX
    };

    Ok(Some(if margin > threshold {
        let phred = ((10.0 * margin / LN_10) as u32).min(255) as u8;
        if m.vdeltas[winner_nr] > 0.0 {
            Decision::VariantRescued(phred)
        } else {
            Decision::PhredConfidence(phred)
        }
    } else {
        Decision::Ambiguous
    }))
}

// -- Public decision type ------------------------------------------------------

#[derive(PartialEq, Debug)]
pub(crate) enum Decision {
    First,
    Last,
    Ambiguous,
    PhredConfidence(u8),
    VariantRescued(u8),
}

// -- Channel message types -----------------------------------------------------

/// Immutable scoring parameters shared across worker threads via `Arc`.
#[derive(Clone)]
pub(super) struct ScoringContext {
    pub(super) penalties: Arc<crate::penalty::Penalty>,
    pub(super) ambiguous_log_threshold: f64,
    pub(super) add_decision_tag: bool,
}

// -- Worker --------------------------------------------------------------------

/// Score a fragment bundle and return the routing decision.
///
/// Mirrors `LineByLine::resolve` + `decide_from_delta` + `apply_ordering` but
/// operates on owned data with no `&mut self`.
pub(super) fn score_bundle(
    best: &mut FragmentBuffer<RecordBuf>,
    stores: &SmallVec<[Option<Arc<dyn StoreTrait>>; 2]>,
    ctx: &ScoringContext,
    scratch: &mut Scratch,
) -> Option<Decision> {
    if best.len() <= 1 {
        return None;
    }

    // Tier 1
    let has_mapped = best.iter().any(|s| !s.is_all_unmapped());
    if has_mapped {
        for i in (0..best.len()).rev() {
            if best[i].is_all_unmapped() {
                best.remove(i);
            }
        }
    }
    if best.len() == 1 {
        return None;
    }
    if !has_mapped {
        return ctx.add_decision_tag.then_some(Decision::Ambiguous);
    }

    let cancel_slot = compute_cancel_slot(best);
    let metrics = compute_metrics(best, cancel_slot, |nr, state, mcfs, cs| {
        let store = stores.get(nr).and_then(|s| s.as_deref());
        let score = score_candidate_owned(state, mcfs, store, ctx, scratch, cs)?;
        Ok((score, scratch.last_variant_delta))
    })
    .ok()?; // wrap_err converts Option<T> closure to Result<T> inside compute_metrics

    apply_tiers(
        best,
        &metrics,
        ctx.ambiguous_log_threshold,
        ctx.add_decision_tag,
        |b, i| {
            b.remove(i);
            Ok(())
        },
    )
    .ok()?
}

/// Score one candidate using its owned [`RecordBuf`]s.
///
/// Mirrors `LineByLine::score_candidate` in `score.rs` but takes owned/Arc
/// data instead of borrowing from `&mut self`.
///
/// # No-variant fast path
///
/// When `store` is `None` the entire `overlapping_multi` call chain is
/// skipped.  The check is a single `Option::is_some()` per candidate — not
/// per segment — so the no-variant case is as lean as possible.
fn score_candidate_owned(
    state: &FragmentState<RecordBuf>,
    mcfs: SmallVec<[MdCigFlags<'_>; READ_CT]>,
    store: Option<&dyn StoreTrait>,
    ctx: &ScoringContext,
    scratch: &mut Scratch,
    cancel_slot: [bool; 2],
) -> Result<f64, Error> {
    use crate::alignment::{stringify_record, Fragment};
    use crate::variant::FragEvalVec;

    let mut segment: SmallVec<[&RecordBuf; READ_CT]> = SmallVec::new();
    let mut seg_mcfs: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();
    let mut mcfs_opt: SmallVec<[Option<MdCigFlags>; READ_CT]> =
        mcfs.into_iter().map(Some).collect();

    // Hoist the variant check out of the per-segment loop.
    // `has_variants` is false for the vast majority of runs.
    let has_variants = store.is_some();

    let mut dvnt: FragEvalVec<'_> = SmallVec::new();
    let mut supplementary_penalty = 0.0;

    for idx in state.order_mates() {
        let flags = state
            .flags(idx)
            .ok_or(Error::NoFlagsForRecordIndex { idx })?;

        if flags.is_secondary() {
            continue;
        }

        // Skip every record (primary or supplementary) belonging to a
        // unanimously-cancelled mate slot: its contribution is provably
        // identical across all competing streams, so it is excluded from
        // the NW segment entirely rather than scored and subtracted out.
        let slot = mate_slot(segment_id(flags));
        if cancel_slot[slot] {
            mcfs_opt[idx] = None; // release the borrow; never consumed
            continue;
        }

        let rec = &state.get_records()[idx];

        // Supplementary alignments contribute BOTH a chimeric-junction
        // penalty (gap_open + chimeric_junction_bases × gap_extend) AND per-base NW
        // scoring.
        if flags.is_supplementary() {
            supplementary_penalty += ctx.penalties.chimeric_junction_penalty;
        }

        if flags.is_unmapped() || !has_variants {
            dvnt.push(SmallVec::new());
        } else {
            let tid = rec
                .ref_seq_id()
                .transpose()?
                .ok_or(Error::MappedRecordNoReferenceSequenceId)?;
            let start = rec
                .alignment_start()
                .ok_or(Error::MappedRecordNoAlignmentStart)?
                .get();
            let cig_len = mcfs_opt[idx]
                .as_ref()
                .ok_or(Error::MdCigFlagsMissing { idx })?
                .get_cigar()
                .len();
            // store is Some here (has_variants == true).
            dvnt.push(
                store
                    .unwrap()
                    .overlapping_multi(tid, start, start + cig_len),
            );
        }

        segment.push(rec);
        seg_mcfs.push(
            mcfs_opt[idx]
                .take()
                .ok_or(Error::MdCigFlagsAlreadyConsumed { idx })?,
        );
    }

    if segment.is_empty() {
        return Ok(supplementary_penalty);
    }

    let base_score = Fragment::new(&ctx.penalties, segment, seg_mcfs)?
        .score(scratch, &mut dvnt)
        .map_err(|e| Error::ScoringError {
            message: e.to_string(),
            state: state
                .get_records()
                .iter()
                .map(stringify_record)
                .collect::<Vec<_>>()
                .join("\n"),
        })?;
    Ok(base_score + supplementary_penalty)
}

// -- LineByLine::process — dispatcher -----------------------------------------

impl<R: SimpleRec> LineByLine<R> {
    /// Original single-threaded process loop.
    pub(crate) fn process_sequential(&mut self, config: &Config) -> Result<(), Error> {
        let mut best: FragmentBuffer<R> = smallvec![];
        let mut i = 0;
        let aln_len = self.aln.len();

        while i != aln_len {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.ingest_record(i, rec, &mut best)? {
                    break;
                }
            }
            let mut decision = None;
            i += 1;
            if i == aln_len {
                if best.is_empty() {
                    break;
                }

                // -- Chimeric pre-check (before normal tournament) ---------
                // If configured chimeric pairs exist, inspect the fragment for
                // complementary mate mapping across species.  When detected,
                // both streams' records are emitted with XC:Z: tags; the normal
                // scoring cascade is skipped for those streams.
                if !self.chimeric_pairs.is_empty() {
                    // Fast path: flag-only per-mate complement check (no CIGAR/MD scanning).
                    if let Some((ca, cb)) =
                        detect_chimeric_mate_complement(&best, &self.chimeric_pairs)
                    {
                        self.emit_chimeric(&mut best, ChimericDecision::Chimeric {
                            stream_a: ca,
                            stream_b: cb,
                            kind: ChimericKind::MateSplit,
                        })?;
                        i = 0;
                        continue;
                    }
                    // Full detection: read-split + false-positive rejection via supplementary MD.
                    let cd = detect_chimeric_event(&best, &self.chimeric_pairs);
                    if matches!(cd, ChimericDecision::Chimeric { .. }) {
                        tracing::debug!(
                            fragment = ?best.first().map(|s| s.first_qname()),
                            "Chimeric event detected"
                        );
                        self.emit_chimeric(&mut best, cd)?;
                        // best is now empty; reset and read next fragment.
                        i = 0;
                        continue;
                    }
                }

                // -- Tournament cascade (N-way round-robin) ----------------
                //
                // Each tier runs in O(N) over all streams simultaneously.
                // Streams are eliminated in backward-sweep discard passes after
                // each tier; emit_winners handles whatever remains in `best`.
                //
                // See `run_tournament` for the full tier description.
                if best.len() > 1 {
                    decision = self.run_tournament(&mut best)?;
                }
                assert!(!best.is_empty());
                self.emit_winners(&mut best, decision)?;
                if let Some(ref mut p) = self.progress { p.tick(); }
                i = 0;
            }
        }
        for i in 0..aln_len {
            if self.aln[i].next_rec()?.is_some() {
                return Err(Error::AlignmentStillHasReads { i });
            }
        }
        if let Some(p) = self.progress.as_ref() { p.finish() }
        config.print_routing_counters(&self.routing_counters, "namesorted");
        Ok(())
    }

    /// N-way round-robin tournament.
    ///
    /// Runs all scoring tiers in a single scan of `best`, storing per-stream
    /// scalar results in fixed-size STACK arrays (N ≤ MAX_STREAMS = 32).
    /// Discards happen after the scan in O(N) backward sweeps; no additional
    /// heap containers are allocated.
    ///
    /// # Tier progression
    ///
    /// ```text
    /// Tier 1   — unmapped filter      O(N) flag checks
    /// Tier 2   — perfect/imperfect    O(N) MCF builds  (shared with tiers below)
    /// Tier 2.5 — match-count max      O(N) ScoreOpIter walks per stream
    /// Tier 3   — NW argmax            O(N) NW score evaluations
    /// ```
    ///
    /// Each tier exits early when a winner is resolved; later tiers are skipped.
    fn run_tournament(&mut self, best: &mut FragmentBuffer<R>) -> Result<Option<Decision>, Error> {
        if best.len() <= 1 {
            return Ok(None);
        }

        // Tier 1
        let has_mapped = best.iter().any(|s| !s.is_all_unmapped());
        if has_mapped {
            for i in (0..best.len()).rev() {
                if best[i].is_all_unmapped() {
                    self.discard_at(best, i)?;
                }
            }
        }
        if best.len() == 1 {
            return Ok(self.add_decision_tag.then_some(Decision::First));
        }
        if !has_mapped {
            return Ok(self.add_decision_tag.then_some(Decision::Ambiguous));
        }

        let cancel_slot = compute_cancel_slot(best);
        let metrics = compute_metrics(best, cancel_slot, |nr, state, mcfs, cs| {
            let store = self
                .aln
                .get(nr)
                .ok_or(Error::NoAlignmentForIndex { aln_idx: nr })?
                .variant_store();
            let score = state.score_state_nw(
                mcfs,
                store.as_deref(),
                &self.penalties,
                &mut self.scratch,
                cs,
            )?;
            Ok((score, self.scratch.last_variant_delta))
        })?;

        apply_tiers(
            best,
            &metrics,
            self.ambiguous_log_threshold,
            self.add_decision_tag,
            |b, i| self.discard_at(b, i),
        )
    }

    /// Emit stream at position `idx` in `best` as filtered output and remove it.
    /// O(N) shift; N ≤ MAX_STREAMS = 32 so this is acceptable.
    fn discard_at(&mut self, best: &mut FragmentBuffer<R>, idx: usize) -> Result<(), Error> {
        let mut loser = best.remove(idx);
        let nr = loser.get_nr();
        loser.drain_records().try_for_each(|r| {
            self.write_record(nr, r.as_record_buf(self.aln[nr].header())?, Some(false))
        })
    }

    pub(crate) fn ingest_record(
        &mut self,
        i: usize,
        rec: R,
        best: &mut FragmentBuffer<R>,
    ) -> Result<bool, Error> {
        if !(self.is_secondary_skipped)(&rec)? {
            let name = rec.name().ok_or(Error::RecordHasNoReadName)?;
            if let Some(new_readname) = (self.is_new_qname)(best, name.as_ref()) {
                if new_readname {
                    #[cfg(test)]
                    if self.aln.is_empty() {
                        return Ok(true);
                    }
                    self.aln[i].un_next(rec)?;
                    return Ok(true);
                }
                // XXX: If always in order, testing last() may be sufficient.
                for state in best.iter_mut().rev() {
                    if state.get_nr() == i {
                        state.add_record(rec)?;
                        return Ok(false);
                    }
                }
            }
            best.push(FragmentState::from_record(rec, i)?);
        }
        Ok(false)
    }

    fn emit_winners(
        &mut self,
        best: &mut FragmentBuffer<R>,
        decision: Option<Decision>,
    ) -> Result<(), Error> {
        let best_state = (best.len() == 1).then_some(true);
        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.drain_records().try_for_each(|r| {
                if best_state.is_none() && (self.is_unmapped_skipped)(&r)? {
                    self.routing_counters[3 + (4 * nr)] += 1;
                    return Ok(());
                }
                let header = self
                    .aln
                    .get(nr)
                    .map(|a| a.header())
                    .ok_or(Error::NoAlignmentForIndex { aln_idx: nr })?;
                let mut rb = RecordBuf::try_from_alignment_record(header, &r)?;
                match decision {
                    Some(Decision::PhredConfidence(d)) => self.add_aux_tags(&mut rb, b"XF", d)?,
                    Some(Decision::VariantRescued(d)) => self.add_aux_tags(&mut rb, b"XR", d)?,
                    _ => {}
                }
                self.write_record(nr, rb, best_state)
            })
        })
    }
}

#[cfg(test)]
mod tests;
