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
use crate::Error;
use crate::alignment::{BaseOp, ScoreOpIter};
use crate::alignment::{FragmentState, MdCigFlags, SimpleRec};
use crate::filter_algorithm::line_by_line::{
    ChimericDecision, MAX_STREAMS, READ_CT, detect_chimeric_event,
};
use crate::variant::StoreTrait;
use crossbeam_channel::{Receiver, Sender, bounded};
use noodles::sam::alignment::RecordBuf;
use smallvec::{SmallVec, smallvec};
use std::f64::consts::LN_10;
use std::sync::Arc;
use std::thread;

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

/// Everything a worker needs to score one fragment.
///
/// All fields are owned/`Arc`-wrapped so the bundle moves across thread
/// boundaries without any lifetime entanglement.
struct FragmentBundle {
    /// Collected alignments, one [`FragmentState`] per stream. Length ≥ 2.
    best: FragmentBuffer<RecordBuf>,
    /// Per-stream variant stores.  Each `Arc` is a cheap clone (atomic
    /// increment) of the one held by [`AlnStream`].  `None` for streams with
    /// no variant file.
    stores: SmallVec<[Option<Arc<dyn StoreTrait>>; 2]>,
    /// Scoring parameters cloned from [`LineByLine`] once at pipeline startup.
    ctx: ScoringContext,
    /// Chimeric pairs to check for complementary mate mapping across species.
    chimeric_pairs: Vec<[usize; 2]>,
}

/// A scored fragment, ready to be written by the IO thread.
struct ScoredFragment {
    /// Only the winning stream(s) remain; losing entries are already stripped
    /// by the worker so the IO thread only needs to route and write.
    best: FragmentBuffer<RecordBuf>,
    decision: Option<Decision>,
}

/// Immutable scoring parameters shared across worker threads via `Arc`.
#[derive(Clone)]
struct ScoringContext {
    penalties: Arc<crate::penalty::Penalty>,
    ambiguous_log_threshold: f64,
    add_decision_tag: bool,
}

// -- Worker --------------------------------------------------------------------

/// Entry point for each scoring worker thread.
///
/// Loops until the IO thread drops the work sender, which closes `rx` and
/// causes `rx.recv()` to return `Err`, exiting cleanly.
fn worker_loop(rx: Receiver<FragmentBundle>, tx: Sender<ScoredFragment>) {
    // Each worker owns its own Scratch (NW DP tables); no contention.
    let mut scratch = Scratch::new();

    while let Ok(bundle) = rx.recv() {
        let FragmentBundle {
            mut best,
            stores,
            ctx,
            chimeric_pairs,
        } = bundle;
        let decision = score_bundle(&mut best, &stores, &ctx, &mut scratch);
        // Ignore send errors: the IO thread may have exited after an error.
        let _ = tx.send(ScoredFragment { best, decision });
    }
}

/// Score a fragment bundle and return the routing decision.
///
/// Mirrors `LineByLine::resolve` + `decide_from_delta` + `apply_ordering` but
/// operates on owned data with no `&mut self`.
fn score_bundle(
    best: &mut FragmentBuffer<RecordBuf>,
    stores: &SmallVec<[Option<Arc<dyn StoreTrait>>; 2]>,
    ctx: &ScoringContext,
    scratch: &mut Scratch,
) -> Option<Decision> {
    if best.len() <= 1 {
        return None;
    }

    // -- Tier 1 ------------------------------------------------------------
    {
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
    }

    // -- Single scan pass (mirrors run_tournament) -------------------------
    let mut is_perfect = [false; MAX_STREAMS];
    let mut match_bases = [0usize; MAX_STREAMS];
    let mut supp_count_arr = [0usize; MAX_STREAMS];
    let mut nw_scores = [f64::NEG_INFINITY; MAX_STREAMS];
    let mut vdeltas = [0.0f64; MAX_STREAMS];
    let mut any_perfect = false;
    let mut any_imperfect = false;

    for i in 0..best.len() {
        let nr = best[i].get_nr();
        let mcfs = best[i].build_mcfs().ok()?;

        let perf = mcfs
            .iter()
            .filter(|m| !m.is_supplementary())
            .all(|m| m.is_perfect());
        is_perfect[nr] = perf;
        if perf {
            any_perfect = true;
        } else {
            any_imperfect = true;
        }

        let mut pmb = 0usize;
        let mut sc = 0usize;
        for mcf in &mcfs {
            if mcf.is_supplementary() {
                continue;
            }
            for op in ScoreOpIter::new(mcf) {
                if matches!(op, Ok(BaseOp::Match)) {
                    pmb += 1;
                }
            }
            sc += mcf.supp_count();
        }
        match_bases[nr] = pmb;
        supp_count_arr[nr] = sc;

        if !perf {
            let store = stores.get(nr).and_then(|s| s.as_deref());
            nw_scores[nr] = score_candidate_owned(&best[i], mcfs, store, ctx, scratch).ok()?;
            vdeltas[nr] = scratch.last_variant_delta;
        }
    }

    // -- Tier 2 ------------------------------------------------------------
    if any_perfect && any_imperfect {
        for i in (0..best.len()).rev() {
            if !is_perfect[best[i].get_nr()] {
                best.remove(i);
            }
        }
        return if best.len() == 1 {
            ctx.add_decision_tag.then_some(Decision::First)
        } else {
            ctx.add_decision_tag.then_some(Decision::Ambiguous)
        };
    }
    if !any_imperfect {
        return ctx.add_decision_tag.then_some(Decision::Ambiguous);
    }

    // -- Tier 2.5 ----------------------------------------------------------
    {
        let max_m = best.iter().map(|s| match_bases[s.get_nr()]).max()?;
        let min_s = best
            .iter()
            .filter(|s| match_bases[s.get_nr()] == max_m)
            .map(|s| supp_count_arr[s.get_nr()])
            .min()?;
        if best.iter().any(|s| {
            let nr = s.get_nr();
            match_bases[nr] < max_m || supp_count_arr[nr] > min_s
        }) {
            for i in (0..best.len()).rev() {
                let nr = best[i].get_nr();
                if match_bases[nr] < max_m || supp_count_arr[nr] > min_s {
                    best.remove(i);
                }
            }
            if best.len() == 1 {
                return ctx.add_decision_tag.then_some(Decision::First);
            }
        }
    }

    // -- Tier 3 ------------------------------------------------------------
    let max_score = best
        .iter()
        .map(|s| nw_scores[s.get_nr()])
        .fold(f64::NEG_INFINITY, f64::max);

    for i in (0..best.len()).rev() {
        if nw_scores[best[i].get_nr()] < max_score - ctx.ambiguous_log_threshold {
            best.remove(i);
        }
    }

    if best.len() != 1 {
        return ctx.add_decision_tag.then_some(Decision::Ambiguous);
    }
    if !ctx.add_decision_tag {
        return None;
    }

    let winner_nr = best[0].get_nr();
    let second_best = nw_scores
        .iter()
        .enumerate()
        .filter(|&(nr, &s)| nr != winner_nr && s > f64::NEG_INFINITY)
        .map(|(_, &s)| s)
        .fold(f64::NEG_INFINITY, f64::max);
    let margin = if second_best.is_finite() {
        max_score - second_best
    } else {
        f64::MAX
    };

    if margin > ctx.ambiguous_log_threshold {
        let phred = if margin.is_finite() {
            ((10.0 * margin / LN_10) as u32).min(255) as u8
        } else {
            255
        };
        Some(if vdeltas[winner_nr] > 0.0 {
            Decision::VariantRescued(phred)
        } else {
            Decision::PhredConfidence(phred)
        })
    } else {
        Some(Decision::Ambiguous)
    }
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
) -> Result<f64, Error> {
    use crate::alignment::{Fragment, stringify_record};
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
    pub(crate) fn process_sequential(&mut self) -> Result<(), Error> {
        let mut best: FragmentBuffer<R> = smallvec![];
        let mut i = 0;

        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.ingest_record(i, rec, &mut best)? {
                    break;
                }
            }
            let mut decision = None;
            i += 1;
            if i == self.aln.len() {
                if best.is_empty() {
                    break;
                }

                // -- Chimeric pre-check (before normal tournament) ---------
                // If configured chimeric pairs exist, inspect the fragment for
                // complementary mate mapping across species.  When detected,
                // both streams' records are emitted with XC:Z: tags; the normal
                // scoring cascade is skipped for those streams.
                if !self.chimeric_pairs.is_empty() {
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
                i = 0;
            }
        }
        while i > 0 {
            i -= 1;
            self.print_counters(i);
            if self.aln[i].next_rec()?.is_some() {
                return Err(Error::AlignmentStillHasReads { i });
            }
        }
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

        // -- Tier 1: unmapped fast-path ------------------------------------
        // Scan once to check whether ANY stream is mapped.
        // If so, unmapped streams are discarded (backward sweep preserves indices).
        {
            let has_mapped = best.iter().any(|s| !s.is_all_unmapped());
            if has_mapped {
                for i in (0..best.len()).rev() {
                    if best[i].is_all_unmapped() {
                        self.discard_at(best, i)?;
                    }
                }
            }
            // Reached by both "all unmapped" (no discard) and "some unmapped" (discarded).
            if best.len() == 1 {
                return Ok(self.add_decision_tag.then_some(Decision::First));
            }
            if !has_mapped {
                // Entire fragment is unmapped in every stream.
                return Ok(self.add_decision_tag.then_some(Decision::Ambiguous));
            }
        }

        // -- Single scan pass: MCF build + all tier metrics -----------------
        //
        // Per-stream results are stored indexed by `get_nr()` (the alignment stream
        // number, 0..MAX_STREAMS).  This survives the backward-discard sweeps that
        // follow: after removing elements from `best`, surviving fragments are
        // looked up by `best[i].get_nr()`, not by their position in `best`.
        //
        // Stack arrays: N ≤ 32, total ≤ ~800 bytes.
        let mut is_perfect = [false; MAX_STREAMS];
        let mut match_bases = [0usize; MAX_STREAMS];
        let mut supp_count_arr = [0usize; MAX_STREAMS];
        let mut nw_scores = [f64::NEG_INFINITY; MAX_STREAMS];
        let mut vdeltas = [0.0f64; MAX_STREAMS];
        let mut any_perfect = false;
        let mut any_imperfect = false;

        for i in 0..best.len() {
            let nr = best[i].get_nr();

            // `build_mcfs` borrows `best[i]` for lifetime 'r.  The borrow is
            // immutable and released either when `mcfs` is moved into
            // `score_candidate` (and eventually dropped inside it) or when
            // `mcfs` drops at the end of this loop body.  In both cases the
            // borrow is gone before the discard sweeps below mutably access
            // `best`.
            let mcfs = best[i].build_mcfs()?;

            // -- Tier 2: perfection check ---------------------------------
            let perf = mcfs
                .iter()
                .filter(|m| !m.is_supplementary())
                .all(|m| m.is_perfect());
            is_perfect[nr] = perf;
            if perf {
                any_perfect = true;
            } else {
                any_imperfect = true;
            }

            // -- Tier 2.5a: primary exact-match count ---------------------
            // Walk each primary's ScoreOpIter once.  Supplementary records
            // list their siblings via SA:Z on the primary; skip them here to
            // avoid double-counting.
            let mut pmb = 0usize; // primary_match_bases
            let mut sc = 0usize; // SA:Z-derived pending supplementary count
            for mcf in &mcfs {
                if mcf.is_supplementary() {
                    continue;
                }
                for op in ScoreOpIter::new(mcf) {
                    if matches!(op, Ok(BaseOp::Match)) {
                        pmb += 1;
                    }
                }
                sc += mcf.supp_count();
            }
            match_bases[nr] = pmb;
            supp_count_arr[nr] = sc;

            // -- Tier 3: NW score ------------------------------------------
            // Computed only for imperfect streams; for perfect streams the
            // MCFs drop here at zero extra cost (the ScoreOpIter walks above
            // were needed for the match-count check anyway).
            if !perf {
                // `mcfs` is moved into `score_candidate`.  The SmallVec and
                // the MdCigFlags inside it (which borrow `best[i]`) are alive
                // for the duration of the call, then dropped.
                nw_scores[nr] = self.score_candidate(&best[i], mcfs, nr)?;
                vdeltas[nr] = self.scratch.last_variant_delta;
            }
            // For perf==true: mcfs drops here.
        }
        // -- End scan pass.  All borrows from elements of `best` have dropped.
        //    Backward-discard sweeps below may now mutably access `best`. --

        // -- Tier 2: apply perfect / imperfect partition --------------------
        if any_perfect && any_imperfect {
            // Discard all imperfect streams; keep only perfect ones.
            for i in (0..best.len()).rev() {
                if !is_perfect[best[i].get_nr()] {
                    self.discard_at(best, i)?;
                }
            }
            // All remaining streams are perfect → indistinguishable here.
            return Ok(if best.len() == 1 {
                self.add_decision_tag.then_some(Decision::First)
            } else {
                self.add_decision_tag.then_some(Decision::Ambiguous)
            });
        }
        if !any_imperfect {
            // Every stream is perfect: genuinely ambiguous.
            return Ok(self.add_decision_tag.then_some(Decision::Ambiguous));
        }

        // -- Tier 2.5: match-count domination ------------------------------
        // N-way generalisation of pairwise AlignSig::subsumes:
        //   dominant set = { nr : match_bases[nr] == max  AND
        //                         supp_count[nr]  == min_among_max }
        // Any stream NOT in the dominant set is discarded.
        // If the dominant set is a strict subset of the competing streams,
        // dominated streams are eliminated here (saves one NW evaluation each).
        {
            let max_m = best
                .iter()
                .map(|s| match_bases[s.get_nr()])
                .max()
                .unwrap_or(0);
            let min_s = best
                .iter()
                .filter(|s| match_bases[s.get_nr()] == max_m)
                .map(|s| supp_count_arr[s.get_nr()])
                .min()
                .unwrap_or(0);

            let any_dominated = best.iter().any(|s| {
                let nr = s.get_nr();
                match_bases[nr] < max_m || supp_count_arr[nr] > min_s
            });

            if any_dominated {
                for i in (0..best.len()).rev() {
                    let nr = best[i].get_nr();
                    if match_bases[nr] < max_m || supp_count_arr[nr] > min_s {
                        self.discard_at(best, i)?;
                    }
                }
                if best.len() == 1 {
                    return Ok(self.add_decision_tag.then_some(Decision::First));
                }
                // Multiple dominant streams: continue to NW for tiebreaking.
            }
        }

        // -- Tier 3: NW argmax ----------------------------------------------
        // All remaining streams are imperfect; NW scores computed in scan pass.
        // Find the highest score; discard every stream below max − threshold.
        let max_score = best
            .iter()
            .map(|s| nw_scores[s.get_nr()])
            .fold(f64::NEG_INFINITY, f64::max);

        for i in (0..best.len()).rev() {
            if nw_scores[best[i].get_nr()] < max_score - self.ambiguous_log_threshold {
                self.discard_at(best, i)?;
            }
        }

        // -- Decision tag ---------------------------------------------------
        if !self.add_decision_tag {
            return Ok(None);
        }
        if best.len() != 1 {
            return Ok(Some(Decision::Ambiguous));
        }
        let winner_nr = best[0].get_nr();

        // Margin = winner's score minus the next-best NW score seen in this
        // tournament round (including eliminated streams).
        let second_best = nw_scores
            .iter()
            .enumerate()
            .filter(|&(nr, &s)| nr != winner_nr && s > f64::NEG_INFINITY)
            .map(|(_, &s)| s)
            .fold(f64::NEG_INFINITY, f64::max);

        let margin = if second_best.is_finite() {
            max_score - second_best
        } else {
            f64::MAX // winner was the only imperfect stream
        };

        Ok(if margin > self.ambiguous_log_threshold {
            let phred = if margin.is_finite() {
                ((10.0 * margin / std::f64::consts::LN_10) as u32).min(255) as u8
            } else {
                255
            };
            Some(if vdeltas[winner_nr] > 0.0 {
                Decision::VariantRescued(phred)
            } else {
                Decision::PhredConfidence(phred)
            })
        } else {
            Some(Decision::Ambiguous)
        })
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

// -- Parallel pipeline — only for LineByLine<RecordBuf> -----------------------

impl LineByLine<RecordBuf> {
    /// Parallel scoring pipeline.
    ///
    /// Collects `Arc<dyn StoreTrait>` clones from each alignment stream (O(1)
    /// atomic increment per stream) before spawning workers, giving every
    /// worker full variant-rescue capability with zero unsafe code and zero
    /// extra allocation.
    pub(crate) fn process_parallel(&mut self) -> Result<(), Error> {
        let n = self.score_threads;
        let cap = n * 2; // bounded channel capacity — natural backpressure

        // -- Collect stores (Arc clones, O(1) each) ------------------------
        let stores: SmallVec<[Option<Arc<dyn StoreTrait>>; 2]> = self
            .aln
            .iter()
            .map(|a| a.variant_store()) // returns Option<Arc<dyn StoreTrait>>
            .collect();

        let any_variants = stores.iter().any(|s| s.is_some());
        tracing::debug!(
            any_variants,
            "Variant stores collected for parallel pipeline"
        );

        let ctx = ScoringContext {
            penalties: Arc::new(self.penalties),
            ambiguous_log_threshold: self.ambiguous_log_threshold,
            add_decision_tag: self.add_decision_tag,
        };

        // -- Spawn workers -------------------------------------------------
        let (work_tx, work_rx) = bounded::<FragmentBundle>(cap);
        let (result_tx, result_rx) = bounded::<ScoredFragment>(cap);

        let workers: Vec<_> = (0..n)
            .map(|id| {
                let rx = work_rx.clone();
                let tx = result_tx.clone();
                thread::Builder::new()
                    .name(format!("xenofilters-scorer-{id}"))
                    .spawn(move || worker_loop(rx, tx))
                    .expect("failed to spawn scorer thread")
            })
            .collect();

        // Drop the extra result_tx held by the spawning thread so result_rx
        // closes automatically when all workers have finished.
        drop(result_tx);

        tracing::info!(
            workers = n,
            any_variants,
            "Parallel scoring pipeline started"
        );

        let io_result = self.parallel_io_loop(work_tx, result_rx, ctx, stores);

        for w in workers {
            if let Err(e) = w.join() {
                tracing::error!(?e, "Scorer worker thread panicked");
            }
        }

        io_result
    }

    /// IO loop: assembles fragments, sends them for scoring, drains results.
    ///
    /// Sending blocks naturally when workers are all busy (bounded channel),
    /// preventing the IO thread from reading too far ahead of the workers.
    /// We interleave result draining on the blocking path to keep the write
    /// side flowing while waiting for channel capacity.
    fn parallel_io_loop(
        &mut self,
        work_tx: Sender<FragmentBundle>,
        result_rx: Receiver<ScoredFragment>,
        ctx: ScoringContext,
        stores: SmallVec<[Option<Arc<dyn StoreTrait>>; 2]>,
    ) -> Result<(), Error> {
        let mut best: FragmentBuffer<RecordBuf> = smallvec![];
        let mut i = 0;

        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.ingest_record(i, rec, &mut best)? {
                    break;
                }
            }

            i += 1;
            if i == self.aln.len() {
                if best.is_empty() {
                    break;
                }

                if best.len() == 1 {
                    // Single-stream fragment — no scoring needed.
                    self.write_scored(ScoredFragment {
                        best: best.drain(..).collect(),
                        decision: None,
                    })?;
                } else {
                    // -- Chimeric pre-check (IO thread, before scoring workers) --
                    if !self.chimeric_pairs.is_empty() {
                        let cd = detect_chimeric_event(&best, &self.chimeric_pairs);
                        if matches!(cd, ChimericDecision::Chimeric { .. }) {
                            // Chimeric routing is a pure IO-thread decision
                            // (no NW DP needed); handle inline.
                            self.emit_chimeric(&mut best, cd)?;
                            i = 0;
                            continue;
                        }
                    }
                    // Clone the Arc per bundle — O(1) per store.
                    let bundle = FragmentBundle {
                        best: best.drain(..).collect(),
                        stores: stores.iter().cloned().collect(),
                        ctx: ctx.clone(),
                        chimeric_pairs: self.chimeric_pairs.clone(),
                    };

                    // Try non-blocking send first.
                    match work_tx.try_send(bundle) {
                        Ok(()) => {}
                        Err(crossbeam_channel::TrySendError::Full(bundle)) => {
                            // Workers are busy: drain one result to free a slot,
                            // then block until the channel has space.
                            match result_rx.recv() {
                                Ok(sf) => self.write_scored(sf)?,
                                Err(_) => return Err(Error::ScorerWorkerExited),
                            }
                            work_tx.send(bundle).map_err(|_| Error::WorkChannelClosed)?;
                        }
                        Err(crossbeam_channel::TrySendError::Disconnected(_)) => {
                            return Err(Error::AllScorerWorkersExited);
                        }
                    }

                    // Non-blocking drain of any already-finished results.
                    loop {
                        match result_rx.try_recv() {
                            Ok(sf) => self.write_scored(sf)?,
                            Err(crossbeam_channel::TryRecvError::Empty) => break,
                            Err(crossbeam_channel::TryRecvError::Disconnected) => {
                                return Err(Error::ScorerWorkerExited);
                            }
                        }
                    }
                }

                i = 0;
            }
        }

        // Signal workers to finish and drain remaining results.
        drop(work_tx);
        for sf in &result_rx {
            self.write_scored(sf)?;
        }

        let mut j = self.aln.len();
        while j > 0 {
            j -= 1;
            self.print_counters(j);
            if self.aln[j].next_rec()?.is_some() {
                return Err(Error::AlignmentStillHasReadsAfterParallelProcessing { j });
            }
        }
        Ok(())
    }

    /// Write a pre-scored fragment through the appropriate output stream(s).
    fn write_scored(&mut self, sf: ScoredFragment) -> Result<(), Error> {
        let ScoredFragment { mut best, decision } = sf;
        let best_state = (best.len() == 1).then_some(true);

        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.drain_records().try_for_each(|mut r| {
                if best_state.is_none() && (self.is_unmapped_skipped)(&r)? {
                    self.routing_counters[3 + (4 * nr)] += 1;
                    return Ok(());
                }
                match decision {
                    Some(Decision::PhredConfidence(d)) => self.add_aux_tags(&mut r, b"XF", d)?,
                    Some(Decision::VariantRescued(d)) => self.add_aux_tags(&mut r, b"XR", d)?,
                    _ => {}
                }
                self.write_record(nr, r, best_state)
            })
        })
    }
}

#[cfg(test)]
mod tests;
