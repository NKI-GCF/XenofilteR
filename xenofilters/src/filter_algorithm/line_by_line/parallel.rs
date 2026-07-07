//! Parallel fragment scoring pipeline for the `namesorted` backend.
//!
//! The IO thread assembles `FragmentBundle`s and dispatches them to a bounded
//! crossbeam channel. Each of `score_threads` worker threads owns a `Scratch`
//! (NW DP tables) and calls `score_bundle`. Results drain back to the IO thread
//! which owns all writers (no Mutex required).

use super::core::{FragmentBuffer, LineByLine, Scratch};
use super::chimeric::{ChimericDecision, detect_chimeric_event, detect_chimeric_mate_complement, ChimericKind};
use super::ordering::{score_bundle, ScoringContext};
use crate::{
    Error,
    variant::StoreTrait,
};
use crossbeam_channel::{Receiver, Sender, bounded};
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{SmallVec, smallvec};
use std::sync::Arc;
use std::thread;
use crate::config::Config;

/// Everything a scoring worker needs for one fragment.
/// All data is owned or `Arc`-wrapped; no lifetime entanglement.
pub(super) struct FragmentBundle {
    /// Collected alignments, one [`FragmentState`] per stream. Length ≥ 2.
    pub(super) best: FragmentBuffer<RecordBuf>,
    /// Per-stream variant stores.  Each `Arc` is a cheap clone (atomic
    /// increment) of the one held by [`AlnStream`].  `None` for streams with
    /// no variant file.
    pub(super) stores: SmallVec<[Option<Arc<dyn StoreTrait>>; 2]>,
    /// Scoring parameters cloned from [`LineByLine`] once at pipeline startup.
    pub(super) ctx: ScoringContext,
    /// Chimeric pairs to check for complementary mate mapping across species.
    pub(super) chimeric_pairs: Vec<[usize; 2]>,
}

/// A scored fragment ready to be written by the IO thread.
pub(super) struct ScoredFragment {
    pub(super) best:     FragmentBuffer<RecordBuf>,
    pub(super) decision: Option<super::ordering::Decision>,
}

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
            chimeric_pairs: _,
        } = bundle;
        let decision = score_bundle(&mut best, &stores, &ctx, &mut scratch);
        // Ignore send errors: the IO thread may have exited after an error.
        let _ = tx.send(ScoredFragment { best, decision });
    }
}

impl LineByLine<RecordBuf> {
    /// Parallel scoring pipeline.
    ///
    /// Collects `Arc<dyn StoreTrait>` clones from each alignment stream (O(1)
    /// atomic increment per stream) before spawning workers, giving every
    /// worker full variant-rescue capability with zero unsafe code and zero
    /// extra allocation.
    pub(crate) fn process_parallel(&mut self, config: &Config) -> Result<(), Error> {
        let n   = self.score_threads;
        let cap = n * 2; // bounded channel capacity — natural backpressure

        // -- Collect stores (Arc clones, O(1) each) ------------------------
        let stores: SmallVec<[Option<Arc<dyn StoreTrait>>; 2]> =
            self.aln.iter().map(|a| a.variant_store()).collect();

        let ctx = ScoringContext {
            penalties:              Arc::new(self.penalties),
            ambiguous_log_threshold: self.ambiguous_log_threshold,
            add_decision_tag:        self.add_decision_tag,
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

        let io_result = self.parallel_io_loop(config, work_tx, result_rx, ctx, stores);

        for w in workers {
            if let Err(e) = w.join() {
                tracing::error!(?e, "scorer thread panicked");
            }
        }
        io_result?;

        Ok(())
    }

    /// IO loop: assembles fragments, sends them for scoring, drains results.
    ///
    /// Sending blocks naturally when workers are all busy (bounded channel),
    /// preventing the IO thread from reading too far ahead of the workers.
    /// We interleave result draining on the blocking path to keep the write
    /// side flowing while waiting for channel capacity.
    pub(super) fn parallel_io_loop(
        &mut self,
        config:    &Config,
        work_tx:   Sender<FragmentBundle>,
        result_rx: Receiver<ScoredFragment>,
        ctx:       ScoringContext,
        stores:    SmallVec<[Option<Arc<dyn StoreTrait>>; 2]>,
    ) -> Result<(), Error> {
        let mut best: FragmentBuffer<RecordBuf> = smallvec![];
        let mut i = 0;
        let aln_len = self.aln.len();

        while i != aln_len {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.ingest_record(i, rec, &mut best)? { break; }
            }
            i += 1;
            if i == aln_len {
                if best.is_empty() { break; }
                if best.len() == 1 {
                    // Single-stream fragment — no scoring needed.
                    self.write_scored(ScoredFragment {
                        best: best.drain(..).collect(),
                        decision: None,
                    })?;
                } else {
                    // -- Chimeric pre-check (IO thread, before scoring workers) --
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
                        Err(_) => return Err(Error::AllScorerWorkersExited),
                    }
                    // Non-blocking drain of any already-finished results.
                    while let Ok(sf) = result_rx.try_recv() { self.write_scored(sf)?; }
                }
                i = 0;
            }
        }
        // Signal workers to finish and drain remaining results.
        drop(work_tx);
        for sf in &result_rx {
            self.write_scored(sf)?;
            if let Some(ref mut p) = self.progress { p.tick(); }
        }

        config.print_routing_counters(&self.routing_counters, "namesorted-parallel");
        for i in 0..aln_len {
            if self.aln[i].next_rec()?.is_some() {
                return Err(Error::AlignmentStillHasReadsAfterParallelProcessing { j: i });
            }
        }
        if let Some(p) = self.progress.as_ref() { p.finish() }
        Ok(())
    }

    fn write_scored(&mut self, sf: ScoredFragment) -> Result<(), Error> {
        use super::ordering::Decision;
        let ScoredFragment { mut best, decision } = sf;
        let best_state = (best.len() == 1).then_some(true);
        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.drain_records().try_for_each(|mut r| {
                if best_state.is_none() && (self.is_unmapped_skipped)(&r)? { return Ok(()); }
                match decision {
                    Some(Decision::PhredConfidence(d)) => self.add_aux_tags(&mut r, b"XF", d)?,
                    Some(Decision::VariantRescued(d))  => self.add_aux_tags(&mut r, b"XR", d)?,
                    _ => {}
                }
                self.write_record(nr, r, best_state)
            })
        })
    }
}
