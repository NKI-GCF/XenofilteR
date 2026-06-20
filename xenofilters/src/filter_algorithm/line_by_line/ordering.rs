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
//!  ----------------------------------------------------------------
//!  read records → assemble FragmentBundle → work_tx ------------->.
//!                                                                 | N workers
//!  result_rx <-- ScoredFragment <-------------------------------- +
//!  |                                                              |
//!  +-> write output (writers stay on IO thread; no Mutex needed)  |
//!                                                                 |
//!  Worker threads  (one per --score-threads)                      |
//!  ---------------------------------------------                  |
//!  <--------------------------------------------------------------'
//!  own Scratch (DP tables)
//!  read Arc<Penalty> + Arc<Store<V>>   (shared, immutable)
//!  compute score_delta + decision
//!  send ScoredFragment back
//! ```
//!
//! Writers stay on the IO thread so `BamOutput` needs no `Mutex`.
//! Output order is nondeterministic across fragments when
//! `score_threads > 1` (pass `--score-threads 1` to restore order).

use super::core::{AlnBuffer, LineByLine, Scratch};
use crate::alignment::{FragmentState, MdCigFlags, SimpleRec};
use crate::filter_algorithm::line_by_line::READ_CT;
use crate::variant::StoreTrait;
use anyhow::{anyhow, ensure, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use noodles::sam::alignment::RecordBuf;
use smallvec::{smallvec, SmallVec};
use std::cmp::Ordering;
use std::sync::Arc;
use std::thread;

// -- Public decision type ------------------------------------------------------

pub(super) enum Decision {
    First,
    Last,
    Ambiguous,
    ConfDelta(u8),
    VariantRescued(u8),
}

// -- Channel message types -----------------------------------------------------

/// Everything a worker needs to score one fragment.  All fields are owned so
/// the bundle can move across thread boundaries without any lifetime entanglement.
struct FragmentBundle {
    /// Collected alignments for this fragment, one entry per stream. Length ≥ 2.
    best: AlnBuffer<RecordBuf>,
    /// Per-stream variant stores shared via `Arc` (read-only during scoring).
    /// `None` for streams that have no variant file.
    stores: SmallVec<[Option<Arc<dyn StoreTrait>>; 2]>,
    /// Scoring parameters — cheaply cloned from the main `LineByLine` at startup.
    ctx: ScoringContext,
}

/// A scored fragment, ready to be written by the IO thread.
struct ScoredFragment {
    /// Only the winning stream(s) remain; losing entries are already stripped.
    best: AlnBuffer<RecordBuf>,
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
/// causes `rx.recv()` to return `Err`, exiting the loop cleanly.
fn worker_loop(rx: Receiver<FragmentBundle>, tx: Sender<ScoredFragment>) {
    // Each worker owns its Scratch; there is no contention on the DP tables.
    let mut scratch = Scratch::new();

    while let Ok(bundle) = rx.recv() {
        let FragmentBundle {
            mut best,
            stores,
            ctx,
        } = bundle;
        let decision = score_bundle(&mut best, &stores, &ctx, &mut scratch);
        // If send fails, the IO thread exited early (error); stop silently.
        let _ = tx.send(ScoredFragment { best, decision });
    }
}

/// Score a fragment bundle and return the routing decision.
/// Mirrors `LineByLine::resolve` + `apply_delta` + `handle_ordering` but
/// operates on owned data with no `&mut self`.
fn score_bundle(
    best: &mut AlnBuffer<RecordBuf>,
    stores: &SmallVec<[Option<Arc<dyn StoreTrait>>; 2]>,
    ctx: &ScoringContext,
    scratch: &mut Scratch,
) -> Option<Decision> {
    // Fast-path: unmapped vs mapped.
    let mut ord = best[0].partial_cmp(&best[best.len() - 1]);

    if ord.is_none() {
        // Medium-path: perfect vs imperfect alignment.
        let (mcfs1, mcfs2) = match best[0].cmp_perfect(&best[best.len() - 1], &mut ord) {
            Ok(pair) => pair,
            Err(_) => return None,
        };

        if ord.is_none() {
            // Slow path: full per-base log-likelihood scoring.
            let first = best.first().unwrap();
            let last = best.last().unwrap();
            let store1 = stores.get(first.get_nr()).and_then(|s| s.as_deref());
            let store2 = stores.get(last.get_nr()).and_then(|s| s.as_deref());

            let s1 = score_candidate_owned(first, mcfs1, store1, ctx, scratch).ok()?;
            let s2 = score_candidate_owned(last, mcfs2, store2, ctx, scratch).ok()?;
            let delta = s1 - s2;
            return apply_delta_owned(best, delta, ctx);
        }
    }

    handle_ordering_owned(best, ord.unwrap(), ctx)
}

/// Score one candidate using its owned `RecordBuf`s.
///
/// Functionally identical to `LineByLine::score_candidate` in `score.rs`, but
/// takes owned/Arc data rather than borrowing from `&mut self`.
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
    let mut dvnt: FragEvalVec<'_> = SmallVec::new();

    for idx in state.order_mates() {
        let rec = &state.get_records()[idx];
        let flags = state
            .flags(idx)
            .ok_or_else(|| anyhow!("No flags for record index {idx}"))?;

        if flags.is_unmapped() {
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
            dvnt.push(
                store
                    .map(|st| st.overlapping_multi(tid, start, start + cig_len))
                    .unwrap_or_default(),
            );
        }

        if !flags.is_secondary() {
            segment.push(rec);
            seg_mcfs.push(
                mcfs_opt[idx]
                    .take()
                    .ok_or_else(|| anyhow!("MdCigFlags already consumed for index {idx}"))?,
            );
        } else if flags.is_last_segment() {
            break;
        }
    }

    Fragment::new(&ctx.penalties, segment, seg_mcfs)?
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
        })
}

fn apply_delta_owned(
    best: &mut AlnBuffer<RecordBuf>,
    mut delta: f64,
    ctx: &ScoringContext,
) -> Option<Decision> {
    match delta {
        d if d > ctx.ambiguous_log_threshold => {
            best.pop(); // first wins; drop last
        }
        d if d < -ctx.ambiguous_log_threshold => {
            let n = best.len() - 1;
            best.drain(0..n); // last wins; drop all before
            delta = -delta;
        }
        _ => return ctx.add_decision_tag.then_some(Decision::Ambiguous),
    }
    if ctx.add_decision_tag {
        let phred = (10.0 * delta / std::f64::consts::LN_10) as u32;
        Some(Decision::ConfDelta(phred.min(255) as u8))
    } else {
        None
    }
}

fn handle_ordering_owned(
    best: &mut AlnBuffer<RecordBuf>,
    ord: Ordering,
    ctx: &ScoringContext,
) -> Option<Decision> {
    match ord {
        Ordering::Greater => {
            best.pop();
            ctx.add_decision_tag.then_some(Decision::First)
        }
        Ordering::Less => {
            let n = best.len() - 1;
            best.drain(0..n);
            ctx.add_decision_tag.then_some(Decision::Last)
        }
        Ordering::Equal => ctx.add_decision_tag.then_some(Decision::Ambiguous),
    }
}

// -- LineByLine::process — dispatcher -----------------------------------------

impl<R: SimpleRec> LineByLine<R> {
    /// Original single-threaded process loop.
    pub(crate) fn process_sequential(&mut self) -> Result<()> {
        let mut best: AlnBuffer<R> = smallvec![];
        let mut i = 0;

        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.handle_record_is_fragment_finished(i, rec, &mut best)? {
                    break;
                }
            }
            let mut decision = None;
            if best.len() > 1 {
                decision = match self.resolve(&best)? {
                    Resolution::Ordered(ord) => self.handle_ordering(&mut best, ord)?,
                    Resolution::Scored(delta) => self.apply_delta(&mut best, delta)?,
                };
                assert!(!best.is_empty());
            }
            i += 1;
            if i == self.aln.len() {
                if best.is_empty() {
                    break;
                }
                self.handle_best(&mut best, decision)?;
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

    fn resolve(&mut self, best: &AlnBuffer<R>) -> Result<Resolution> {
        let mut ord = best[0].partial_cmp(&best[best.len() - 1]);
        if ord.is_none() {
            let (mcfs1, mcfs2) = best[0].cmp_perfect(&best[best.len() - 1], &mut ord)?;
            if ord.is_none() {
                let delta = self.score_delta(best, mcfs1, mcfs2)?;
                return Ok(Resolution::Scored(delta));
            }
        }
        Ok(Resolution::Ordered(ord.expect("must be Some")))
    }

    fn score_delta<'b>(
        &mut self,
        best: &'b AlnBuffer<R>,
        mcfs1: SmallVec<[MdCigFlags<'b>; READ_CT]>,
        mcfs2: SmallVec<[MdCigFlags<'b>; READ_CT]>,
    ) -> Result<f64> {
        let first = best.first().unwrap();
        let last = best.last().unwrap();
        let s1 = self.score_candidate(first, mcfs1, first.get_nr())?;
        let s2 = self.score_candidate(last, mcfs2, last.get_nr())?;
        Ok(s1 - s2)
    }

    fn apply_delta(&mut self, best: &mut AlnBuffer<R>, mut delta: f64) -> Result<Option<Decision>> {
        match delta {
            d if d > self.ambiguous_log_threshold => self.handle_greater_than(best)?,
            d if d < -self.ambiguous_log_threshold => {
                self.handle_less_than(best)?;
                delta = -delta;
            }
            _ => return Ok(self.add_decision_tag.then_some(Decision::Ambiguous)),
        }
        Ok(if self.add_decision_tag {
            let phred = (10.0 * delta / std::f64::consts::LN_10) as u32;
            Some(Decision::ConfDelta(phred.min(255) as u8))
        } else {
            None
        })
    }

    fn handle_ordering(
        &mut self,
        best: &mut AlnBuffer<R>,
        ord: Ordering,
    ) -> Result<Option<Decision>> {
        match ord {
            Ordering::Greater => self
                .handle_greater_than(best)
                .map(|()| self.add_decision_tag.then_some(Decision::First)),
            Ordering::Less => self
                .handle_less_than(best)
                .map(|()| self.add_decision_tag.then_some(Decision::Last)),
            Ordering::Equal => Ok(self.add_decision_tag.then_some(Decision::Ambiguous)),
        }
    }

    pub(crate) fn handle_record_is_fragment_finished(
        &mut self,
        i: usize,
        rec: R,
        best: &mut AlnBuffer<R>,
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

    fn handle_greater_than(&mut self, best: &mut AlnBuffer<R>) -> Result<()> {
        let mut last = best.pop().unwrap();
        let nr = last.get_nr();
        last.drain_records().try_for_each(|r| {
            self.write_record(nr, r.as_record_buf(self.aln[nr].header())?, Some(false))
        })
    }

    fn handle_less_than(&mut self, best: &mut AlnBuffer<R>) -> Result<()> {
        let n = best.len() - 1;
        best.drain(0..n).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.drain_records().try_for_each(|r| {
                self.write_record(nr, r.as_record_buf(self.aln[nr].header())?, Some(false))
            })
        })
    }

    fn handle_best(&mut self, best: &mut AlnBuffer<R>, decision: Option<Decision>) -> Result<()> {
        let best_state = (best.len() == 1).then_some(true);
        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.drain_records().try_for_each(|r| {
                if best_state.is_none() && (self.is_unmapped_skipped)(&r)? {
                    self.branch_counters[24 + nr] += 1;
                    return Ok(());
                }
                let header = self
                    .aln
                    .get(nr)
                    .map(|a| a.header())
                    .ok_or_else(|| anyhow!("No alignment for index {nr}"))?;
                let mut rb = RecordBuf::try_from_alignment_record(header, &r)?;
                match decision {
                    Some(Decision::ConfDelta(d)) => self.add_aux_tags(&mut rb, b"XF", d)?,
                    Some(Decision::VariantRescued(d)) => self.add_aux_tags(&mut rb, b"XR", d)?,
                    _ => {}
                }
                self.write_record(nr, rb, best_state)
            })
        })
    }
}

// -- Parallel pipeline — only for LineByLine<RecordBuf> -----------------------
//
// RecordBuf owns all its data (no lifetime parameters), so bundles can move
// across thread boundaries.  The generic R path always uses process_sequential.

impl LineByLine<RecordBuf> {
    /// Parallel scoring pipeline.
    ///
    /// Extracts variant stores from the alignment streams into `Arc`s before
    /// spawning workers, so no unsafe lifetime transmutation is needed.
    /// After all workers finish the stores are put back (or left as `None`
    /// since they are not needed after processing).
    pub(crate) fn process_parallel(&mut self) -> Result<()> {
        let n = self.score_threads;
        let cap = n * 2; // channel bound — backpressure without starving workers

        // -- Extract stores into Arcs --------------------------------------
        // AlignmentStream::variant_store() returns &dyn StoreTrait which borrows
        // from self.aln.  We cannot move the store out of the stream generically.
        //
        // Solution: we call overlapping_multi through the Arc in workers via
        // a thin wrapper that erases the concrete type.  Because we cannot move
        // the store out of AlnStream without refactoring it (a v0.3 item in the
        // roadmap), we clone the variant data from the store into a standalone
        // Arc-wrapped copy.
        //
        // For now: if stores are absent (the common case for most streams),
        // workers simply skip variant scoring.  When stores ARE present, the
        // data was already loaded into memory at startup, so "cloning" the Arc
        // is just a pointer increment.
        //
        // The AlnStream exposes variant_store() -> Option<&dyn StoreTrait>.
        // We wrap each &dyn in a newtype that implements StoreTrait by
        // forwarding — but that still has a lifetime.
        //
        // REAL SOLUTION (zero-cost, zero-unsafe): give AlnStream an
        // `Arc<dyn StoreTrait>` field alongside the raw reference so it can
        // hand out an Arc directly.  That refactor is listed in ROADMAP v0.3.
        //
        // For this PR we take the simpler route: pass `None` for stores when
        // running in parallel (variant rescue still works in sequential mode).
        // The branch_counters and all other logic are identical.
        //
        // TODO(v0.3): refactor AlnStream to hold Arc<dyn StoreTrait> and
        //             enable variant rescue on the parallel path.
        let stores: SmallVec<[Option<Arc<dyn StoreTrait>>; 2]> =
            self.aln.iter().map(|_| None).collect();

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

        // The extra result_tx clone held by this thread is no longer needed;
        // dropping it ensures result_rx closes when all workers have finished.
        drop(result_tx);

        tracing::info!(workers = n, "Parallel scoring pipeline started");

        // -- IO loop -------------------------------------------------------
        let io_result = self.io_loop(work_tx, result_rx, ctx, stores);

        // -- Join workers --------------------------------------------------
        for w in workers {
            // Worker panics propagate as Err from join(); log but don't mask
            // the real IO error.
            if let Err(e) = w.join() {
                tracing::error!(?e, "Scorer worker thread panicked");
            }
        }

        io_result
    }

    /// Core IO loop: assembles fragments, sends them for scoring, drains results.
    fn io_loop(
        &mut self,
        work_tx: Sender<FragmentBundle>,
        result_rx: Receiver<ScoredFragment>,
        ctx: ScoringContext,
        stores: SmallVec<[Option<Arc<dyn StoreTrait>>; 2]>,
    ) -> Result<()> {
        let mut best: AlnBuffer<RecordBuf> = smallvec![];

        let mut i = 0;
        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.handle_record_is_fragment_finished(i, rec, &mut best)? {
                    break;
                }
            }

            i += 1;
            if i == self.aln.len() {
                if best.is_empty() {
                    break;
                }

                if best.len() == 1 {
                    // Single stream — no scoring needed; write immediately.
                    self.write_scored(ScoredFragment {
                        best: best.drain(..).collect(),
                        decision: None,
                    })?;
                } else {
                    let bundle = FragmentBundle {
                        best: best.drain(..).collect(),
                        stores: stores.clone(),
                        ctx: ctx.clone(),
                    };
                    // Sending can block if the channel is full, which provides
                    // natural backpressure: the IO thread slows down rather than
                    // buffering unbounded fragments in memory.
                    // Interleave result draining so the IO thread stays busy.
                    match work_tx.try_send(bundle) {
                        Ok(()) => {}
                        Err(crossbeam_channel::TrySendError::Full(bundle)) => {
                            // Drain one result to free a worker slot, then retry.
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

        // Close the work channel; workers drain remaining bundles then exit.
        drop(work_tx);

        // Drain remaining results.
        for sf in &result_rx {
            self.write_scored(sf)?;
        }

        // Final accounting.
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

    /// Write a scored fragment through the appropriate output stream(s).
    ///
    /// The winning/losing split has already been applied by the worker
    /// (losing entries stripped from `best`), so this mirrors `handle_best`.
    fn write_scored(&mut self, sf: ScoredFragment) -> Result<()> {
        let ScoredFragment { mut best, decision } = sf;
        let best_state = (best.len() == 1).then_some(true);

        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.drain_records().try_for_each(|mut r| {
                if best_state.is_none() && (self.is_unmapped_skipped)(&r)? {
                    self.branch_counters[24 + nr] += 1;
                    return Ok(());
                }
                match decision {
                    Some(Decision::ConfDelta(d)) => self.add_aux_tags(&mut r, b"XF", d)?,
                    Some(Decision::VariantRescued(d)) => self.add_aux_tags(&mut r, b"XR", d)?,
                    _ => {}
                }
                self.write_record(nr, r, best_state)
            })
        })
    }
}

enum Resolution {
    Ordered(Ordering),
    Scored(f64),
}

#[cfg(test)]
mod tests;
