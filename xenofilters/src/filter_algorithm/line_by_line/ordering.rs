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
use crate::alignment::{pre_assess_alignments, PreAssessResult};
use crate::alignment::{FragmentState, MdCigFlags, SimpleRec};
use crate::filter_algorithm::line_by_line::READ_CT;
use crate::variant::StoreTrait;
use anyhow::{anyhow, ensure, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use noodles::sam::alignment::RecordBuf;
use smallvec::{smallvec, SmallVec};
use std::cmp::Ordering;
use std::f64::consts::LN_10;
use std::sync::Arc;
use std::thread;

// -- Public decision type ------------------------------------------------------

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
    while best.len() > 1 {
        let len_before = best.len();

        let mut ord = best[0].partial_cmp(&best[1]);
        if ord.is_none() {
            let (mcfs1, mcfs2) = match best[0].cmp_perfect(&best[1], &mut ord) {
                Ok(pair) => pair,
                Err(_) => return None,
            };
            if ord.is_none() {
                if let PreAssessResult::EarlyDecision(pa_ord) =
                    pre_assess_alignments(&mcfs1, &mcfs2)
                {
                    drop(mcfs1);
                    drop(mcfs2);
                    let dec = apply_ordering_owned(best, pa_ord, ctx);
                    if best.len() == len_before {
                        return dec;
                    }
                    continue;
                }
                let a = &best[0];
                let b = &best[1];
                let store1 = stores.get(a.get_nr()).and_then(|s| s.as_deref());
                let store2 = stores.get(b.get_nr()).and_then(|s| s.as_deref());
                let s1 = score_candidate_owned(a, mcfs1, store1, ctx, scratch).ok()?;
                let s1_vd = scratch.last_variant_delta;
                let s2 = score_candidate_owned(b, mcfs2, store2, ctx, scratch).ok()?;
                let s2_vd = scratch.last_variant_delta;
                let dec = decide_from_delta_owned(best, s1 - s2, s1_vd, s2_vd, ctx);
                if best.len() == len_before {
                    return dec;
                }
                continue;
            }
        }
        let dec = apply_ordering_owned(best, ord.unwrap(), ctx);
        if best.len() == len_before {
            return dec;
        }
    }
    None
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
) -> Result<f64> {
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
            .ok_or_else(|| anyhow!("No flags for record index {idx}"))?;

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
                .ok_or_else(|| anyhow!("Mapped record has no reference sequence ID"))?;
            let start = rec
                .alignment_start()
                .ok_or_else(|| anyhow!("Mapped record has no alignment start"))?
                .get();
            let cig_len = mcfs_opt[idx]
                .as_ref()
                .ok_or_else(|| anyhow!("MdCigFlags missing for index {idx}"))?
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
                .ok_or_else(|| anyhow!("MdCigFlags already consumed for index {idx}"))?,
        );
    }

    if segment.is_empty() {
        return Ok(supplementary_penalty);
    }

    let base_score = Fragment::new(&ctx.penalties, segment, seg_mcfs)?
        .score(scratch, &mut dvnt)
        .map_err(|e| {
            anyhow!(
                "Scoring error: {e}\n{}",
                state
                    .get_records()
                    .iter()
                    .map(stringify_record)
                    .collect::<Vec<_>>()
                    .join("\n")
            )
        })?;
    Ok(base_score + supplementary_penalty)
}

fn decide_from_delta_owned(
    best: &mut FragmentBuffer<RecordBuf>,
    mut delta: f64,
    s1_vd: f64,
    s2_vd: f64,
    ctx: &ScoringContext,
) -> Option<Decision> {
    match delta {
        d if d > ctx.ambiguous_log_threshold => {
            best.remove(1); // best[0] wins; discard challenger
            if ctx.add_decision_tag {
                let phred = (10.0 * delta / LN_10).min(255.0) as u8;
                return Some(if s1_vd > 0.0 {
                    Decision::VariantRescued(phred)
                } else {
                    Decision::PhredConfidence(phred)
                });
            }
        }
        d if d < -ctx.ambiguous_log_threshold => {
            best.remove(0); // best[1] wins; discard leader
            delta = -delta;
            if ctx.add_decision_tag {
                let phred = (10.0 * delta / LN_10).min(255.0) as u8;
                return Some(if s2_vd > 0.0 {
                    Decision::VariantRescued(phred)
                } else {
                    Decision::PhredConfidence(phred)
                });
            }
        }
        _ => return ctx.add_decision_tag.then_some(Decision::Ambiguous),
    }
    None
}

fn apply_ordering_owned(
    best: &mut FragmentBuffer<RecordBuf>,
    ord: Ordering,
    ctx: &ScoringContext,
) -> Option<Decision> {
    match ord {
        Ordering::Greater => {
            best.remove(1);
            ctx.add_decision_tag.then_some(Decision::First)
        }
        Ordering::Less => {
            best.remove(0);
            ctx.add_decision_tag.then_some(Decision::Last)
        }
        Ordering::Equal => ctx.add_decision_tag.then_some(Decision::Ambiguous),
    }
}

// -- LineByLine::process — dispatcher -----------------------------------------

impl<R: SimpleRec> LineByLine<R> {
    /// Original single-threaded process loop.
    pub(crate) fn process_sequential(&mut self) -> Result<()> {
        let mut best: FragmentBuffer<R> = smallvec![];
        let mut i = 0;

        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.ingest_record(i, rec, &mut best)? {
                    break;
                }
            }
            let mut decision = None;
            if best.len() > 1 {
                // Champion-vs-challenger sequential tournament.
                //
                // Each round compares best[0] (current leader) vs best[1] (next challenger).
                //   • Leader wins  → discard_at(1); leader stays at position 0.
                //   • Challenger wins → discard_at(0); challenger shifts to position 0
                //                       and becomes the new leader for the next round.
                //   • Tie detected → break; ALL remaining streams go to ambiguous output.
                //
                // Correctness: visiting pairs in order means a stream that beats [0]
                // immediately becomes the new [0] and must defeat every subsequent
                // challenger before being declared the winner.  No stream is silently
                // skipped.
                //
                // Known limitation: if best[0] = best[1] (tie) but best[2] > both,
                // best[2] never participates — all three go to ambiguous. A full
                // round-robin N-way tournament is tracked in ROADMAP v0.4.
                while best.len() > 1 {
                    let len_before = best.len();
                    decision = match self.resolve(&best)? {
                        Resolution::Early(ord) => self.apply_ordering(&mut best, ord)?,
                        Resolution::NwDelta(d, s1_vd, s2_vd) => {
                            self.decide_from_delta(&mut best, d, s1_vd, s2_vd)?
                        }
                    };
                    // If best.len() did not decrease, no stream was eliminated:
                    // the result was Equal / within ambiguous_threshold. Break to
                    // avoid an infinite loop; emit_winners will route all remaining
                    // streams to ambiguous output.
                    if best.len() == len_before {
                        break;
                    }
                }
                assert!(!best.is_empty());
            }
            i += 1;
            if i == self.aln.len() {
                if best.is_empty() {
                    break;
                }
                self.emit_winners(&mut best, decision)?;
                i = 0;
            }
        }
        while i > 0 {
            i -= 1;
            self.print_counters(i);
            ensure!(
                self.aln[i].next_rec()?.is_none(),
                "alignment {i} still has reads"
            );
        }
        Ok(())
    }

    fn compute_score_delta<'b>(
        &mut self,
        state_a: &'b FragmentState<R>,
        state_b: &'b FragmentState<R>,
        mcfs1: SmallVec<[MdCigFlags<'b>; READ_CT]>,
        mcfs2: SmallVec<[MdCigFlags<'b>; READ_CT]>,
    ) -> Result<(f64, f64, f64)> {
        let s1 = self.score_candidate(state_a, mcfs1, state_a.get_nr())?;
        let s1_vd = self.scratch.last_variant_delta;
        let s2 = self.score_candidate(state_b, mcfs2, state_b.get_nr())?;
        let s2_vd = self.scratch.last_variant_delta;
        Ok((s1 - s2, s1_vd, s2_vd))
    }

    fn resolve(&mut self, best: &FragmentBuffer<R>) -> Result<Resolution> {
        debug_assert!(best.len() >= 2, "resolve() requires at least 2 candidates");
        // Always compare the current leader (best[0]) against the next
        // challenger (best[1]).  The tournament loop in process_sequential
        // ensures these are the only two streams in play per round.
        let a = &best[0];
        let b = &best[1];
        let mut ord = a.partial_cmp(b);
        if ord.is_none() {
            let (mcfs1, mcfs2) = a.cmp_perfect(b, &mut ord)?;
            if ord.is_none() {
                if let PreAssessResult::EarlyDecision(pa_ord) =
                    pre_assess_alignments(&mcfs1, &mcfs2)
                {
                    return Ok(Resolution::Early(pa_ord));
                }
                let (delta, s1_vd, s2_vd) = self.compute_score_delta(a, b, mcfs1, mcfs2)?;
                return Ok(Resolution::NwDelta(delta, s1_vd, s2_vd));
            }
        }
        Ok(Resolution::Early(ord.unwrap()))
    }

    fn decide_from_delta(
        &mut self,
        best: &mut FragmentBuffer<R>,
        mut delta: f64,
        s1_vd: f64,
        s2_vd: f64,
    ) -> Result<Option<Decision>> {
        match delta {
            d if d > self.ambiguous_log_threshold => {
                self.discard_at(best, 1)?; // best[0] wins; discard challenger
                if self.add_decision_tag {
                    let phred = (10.0 * delta / LN_10).min(255.0) as u8;
                    return Ok(Some(if s1_vd > 0.0 {
                        Decision::VariantRescued(phred)
                    } else {
                        Decision::PhredConfidence(phred)
                    }));
                }
            }
            d if d < -self.ambiguous_log_threshold => {
                self.discard_at(best, 0)?; // best[1] wins; discard leader
                delta = -delta;
                if self.add_decision_tag {
                    let phred = (10.0 * delta / LN_10).min(255.0) as u8;
                    return Ok(Some(if s2_vd > 0.0 {
                        Decision::VariantRescued(phred)
                    } else {
                        Decision::PhredConfidence(phred)
                    }));
                }
            }
            _ => return Ok(self.add_decision_tag.then_some(Decision::Ambiguous)),
        }
        Ok(None)
    }

    /// Emit stream at position `idx` in `best` as filtered output and remove it.
    /// O(N) shift; N ≤ MAX_STREAMS = 32 so this is acceptable.
    fn discard_at(&mut self, best: &mut FragmentBuffer<R>, idx: usize) -> Result<()> {
        let mut loser = best.remove(idx);
        let nr = loser.get_nr();
        loser.drain_records().try_for_each(|r| {
            self.write_record(nr, r.as_record_buf(self.aln[nr].header())?, Some(false))
        })
    }

    fn apply_ordering(
        &mut self,
        best: &mut FragmentBuffer<R>,
        ord: Ordering,
    ) -> Result<Option<Decision>> {
        match ord {
            Ordering::Greater => {
                // best[0] (current leader) beats best[1] (current challenger).
                self.discard_at(best, 1)?;
                Ok(self.add_decision_tag.then_some(Decision::First))
            }
            Ordering::Less => {
                // best[1] (challenger) beats best[0] (leader).
                // Removing position 0 shifts the challenger to position 0,
                // making it the new leader for the next tournament round.
                self.discard_at(best, 0)?;
                Ok(self.add_decision_tag.then_some(Decision::Last))
            }
            Ordering::Equal => {
                // Tie: caller breaks the tournament loop; remaining streams
                // (including any not yet compared) all route to ambiguous output.
                Ok(self.add_decision_tag.then_some(Decision::Ambiguous))
            }
        }
    }

    pub(crate) fn ingest_record(
        &mut self,
        i: usize,
        rec: R,
        best: &mut FragmentBuffer<R>,
    ) -> Result<bool> {
        if !(self.is_secondary_skipped)(&rec)? {
            let name = rec
                .name()
                .ok_or_else(|| anyhow!("Record has no read name"))?;
            if let Some(new_readname) = (self.is_new_qname)(best, name.as_ref()) {
                if new_readname {
                    #[cfg(test)]
                    if self.aln.is_empty() {
                        return Ok(true);
                    }
                    self.aln[i].un_next(rec)?;
                    return Ok(true);
                }
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
    ) -> Result<()> {
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
                    .ok_or_else(|| anyhow!("No alignment for index {nr}"))?;
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
    pub(crate) fn process_parallel(&mut self) -> Result<()> {
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
    ) -> Result<()> {
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
                    // Clone the Arc per bundle — O(1) per store.
                    let bundle = FragmentBundle {
                        best: best.drain(..).collect(),
                        stores: stores.iter().cloned().collect(),
                        ctx: ctx.clone(),
                    };

                    // Try non-blocking send first.
                    match work_tx.try_send(bundle) {
                        Ok(()) => {}
                        Err(crossbeam_channel::TrySendError::Full(bundle)) => {
                            // Workers are busy: drain one result to free a slot,
                            // then block until the channel has space.
                            match result_rx.recv() {
                                Ok(sf) => self.write_scored(sf)?,
                                Err(_) => return Err(anyhow!("Scorer worker exited unexpectedly")),
                            }
                            work_tx
                                .send(bundle)
                                .map_err(|_| anyhow!("Work channel closed unexpectedly"))?;
                        }
                        Err(crossbeam_channel::TrySendError::Disconnected(_)) => {
                            return Err(anyhow!("All scorer workers exited unexpectedly"))
                        }
                    }

                    // Non-blocking drain of any already-finished results.
                    loop {
                        match result_rx.try_recv() {
                            Ok(sf) => self.write_scored(sf)?,
                            Err(crossbeam_channel::TryRecvError::Empty) => break,
                            Err(crossbeam_channel::TryRecvError::Disconnected) => {
                                return Err(anyhow!("Scorer workers exited unexpectedly"))
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
            ensure!(
                self.aln[j].next_rec()?.is_none(),
                "alignment {j} still has reads after parallel processing"
            );
        }
        Ok(())
    }

    /// Write a pre-scored fragment through the appropriate output stream(s).
    fn write_scored(&mut self, sf: ScoredFragment) -> Result<()> {
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

/// Outcome of `resolve()`: either an early decision from Tiers 1–2.5,
/// or a full NW score triple from Tier 3.
enum Resolution {
    /// Tiers 1, 2, or 2.5 produced a definitive ordering without NW DP.
    Early(Ordering),
    /// Full NW scoring completed; carries (delta, s1_variant_delta, s2_variant_delta).
    NwDelta(f64, f64, f64),
}

#[cfg(test)]
mod tests;
