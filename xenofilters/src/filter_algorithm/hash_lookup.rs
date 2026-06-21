//! [`HashLookup`] — fragment-matching backend for arbitrary (non-collated) BAM
//! input.
//!
//! Records from both streams are inserted into a shared `NameTable` keyed by
//! canonical read name. Once both streams have contributed a primary alignment
//! for a given name, the fragment is scored and staged for emission in driving-
//! stream order.
//!
//! Memory is proportional to the number of in-flight (incomplete) fragments.
//! Worst case O(N). Expected case much smaller when both BAMs are sorted by
//! the same coordinates.
//!
//! Output order follows the driving stream (stream 0).
//! Single-threaded only.

pub(crate) mod assemble;
pub(crate) mod stage;
#[cfg(test)]
pub(crate) mod tests;

use crate::alignment::{stringify_record, Fragment, FragmentState, MdCigFlags, SimpleRec};
use crate::aln_stream::AlignmentStream;
use crate::config::{Config, StripReadSuffix};
use crate::filter_algorithm::line_by_line::{ordering::Decision, Scratch, READ_CT};
use crate::penalty::Penalty;
use crate::variant::FragEvalVec;
use anyhow::{anyhow, Result};
use assemble::{insert, NameTable};
use smallvec::SmallVec;
use stage::StagedOutput;
use std::cmp::Ordering;

/// A fragment that has been scored and is waiting for in-order emission.
pub(crate) struct ScoredFragment<R: SimpleRec> {
    /// The winning fragment (written to best / ambiguous output).
    pub(crate) winner: FragmentState<R>,
    /// The losing fragment (written to filtered output), if any.
    /// None for ambiguous pairs — both fragments are treated as "winner".
    pub(crate) loser: Option<FragmentState<R>>,
    pub(crate) decision: Option<Decision>,
}

pub(crate) struct HashLookup<R: SimpleRec> {
    aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    table: NameTable<R>,
    staged: StagedOutput<R>,
    seq_counter: u64,
    /// Records per stream read so far, used to interleave reads.
    lookahead_limit: usize,
    penalties: Penalty,
    scratch: Scratch,
    pub(crate) branch_counters: [u64; 32],
    add_decision_tag: bool,
    ambiguous_log_threshold: f64,
    strip: StripReadSuffix,
}

impl<R: SimpleRec> HashLookup<R> {
    pub(crate) fn new(
        config: Config,
        mut aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    ) -> Result<Self> {
        assert_eq!(
            aln.len(),
            2,
            "HashLookup requires exactly 2 alignment streams"
        );

        let ambiguous_log_threshold = match config.ambiguous_threshold {
            0 => 0.0,
            t => (t as f64) * std::f64::consts::LN_10 / 10.0,
        };

        for (i, a) in aln.iter_mut().enumerate() {
            a.init_writers(&config, i)?;
        }

        Ok(Self {
            aln,
            table: NameTable::new(),
            staged: StagedOutput::new(),
            seq_counter: 0,
            lookahead_limit: 100_000,
            penalties: config.to_penalties(),
            scratch: Scratch::new(),
            branch_counters: [0u64; 32],
            add_decision_tag: config.add_decision_tag,
            ambiguous_log_threshold,
            strip: config.strip_read_suffix,
        })
    }

    pub(crate) fn process(&mut self) -> Result<()> {
        // Interleave reads from both streams until both are exhausted.
        let mut exhausted = [false; 2];
        loop {
            let mut progress = false;
            for nr in 0..2usize {
                if exhausted[nr] {
                    continue;
                }
                match self.aln[nr].next_rec()? {
                    None => {
                        exhausted[nr] = true;
                    }
                    Some(rec) => {
                        progress = true;
                        self.ingest(rec, nr)?;
                    }
                }
            }
            if !progress {
                break;
            }
            // Flush staged output whenever possible.
            self.staged.flush(
                self.aln.as_mut_slice(),
                &mut self.branch_counters,
                self.add_decision_tag,
            )?;
        }

        // Any fragments remaining in the table are unmatched (malformed input).
        let unmatched: Vec<_> = self.table.drain().collect();
        for (_, pending) in unmatched {
            let seq_nr = pending.seq_nr;
            let frag = if pending.driving.is_ready() {
                pending.into_driving_fragment()?
            } else {
                // Only lookup side arrived — emit from stream 1.
                // Reconstruct manually.
                let (_, b) = pending
                    .into_pair()
                    .unwrap_or_else(|_| panic!("both sides empty"));
                b
            };
            let sf = ScoredFragment {
                winner: frag,
                loser: None,
                decision: None,
            };
            self.staged.push(seq_nr, sf);
        }

        self.staged.flush_all(
            self.aln.as_mut_slice(),
            &mut self.branch_counters,
            self.add_decision_tag,
        )?;

        self.print_counters();
        Ok(())
    }

    fn ingest(&mut self, rec: R, nr: usize) -> Result<()> {
        let (key, complete) = insert(&mut self.table, rec, nr, self.strip, &mut self.seq_counter);
        if complete {
            let pending = self.table.remove(&key).unwrap();
            let seq_nr = pending.seq_nr;
            let (a, b) = pending.into_pair()?;
            let sf = self.score_pair(a, b)?;
            self.staged.push(seq_nr, sf);
        }
        Ok(())
    }

    fn score_pair(
        &mut self,
        a: FragmentState<R>,
        b: FragmentState<R>,
    ) -> Result<ScoredFragment<R>> {
        let mut ord = a.partial_cmp(&b);
        let (mcfs1, mcfs2) = if ord.is_none() {
            a.cmp_perfect(&b, &mut ord)?
        } else {
            (SmallVec::new(), SmallVec::new())
        };

        enum Res {
            Ordered(Ordering),
            Scored(f64),
        }

        let res = if ord.is_none() {
            let s1 = self.score_one(&a, mcfs1, 0)?;
            let s2 = self.score_one(&b, mcfs2, 1)?;
            Res::Scored(s1 - s2)
        } else {
            Res::Ordered(ord.unwrap())
        };

        let sf = match res {
            Res::Ordered(Ordering::Greater) => {
                let dec = self.phred_decision_first();
                ScoredFragment {
                    winner: a,
                    loser: Some(b),
                    decision: dec,
                }
            }
            Res::Ordered(Ordering::Less) => {
                let dec = self.phred_decision_last();
                ScoredFragment {
                    winner: b,
                    loser: Some(a),
                    decision: dec,
                }
            }
            Res::Ordered(Ordering::Equal) => {
                let dec = self.add_decision_tag.then_some(Decision::Ambiguous);
                // Both go to ambiguous — encode as winner=a, loser=None, and
                // emit b separately by pushing a second ScoredFragment.
                // Simpler: store both in winner using a sentinel loser=Some(b)
                // with ambiguous decision. stage.rs handles this via is_ambiguous check.
                ScoredFragment {
                    winner: a,
                    loser: Some(b),
                    decision: dec,
                }
            }
            Res::Scored(delta) if delta > self.ambiguous_log_threshold => {
                let dec = self.phred_delta(delta);
                ScoredFragment {
                    winner: a,
                    loser: Some(b),
                    decision: dec,
                }
            }
            Res::Scored(delta) if delta < -self.ambiguous_log_threshold => {
                let dec = self.phred_delta(-delta);
                ScoredFragment {
                    winner: b,
                    loser: Some(a),
                    decision: dec,
                }
            }
            Res::Scored(_) => {
                let dec = self.add_decision_tag.then_some(Decision::Ambiguous);
                ScoredFragment {
                    winner: a,
                    loser: Some(b),
                    decision: dec,
                }
            }
        };
        Ok(sf)
    }

    fn phred_decision_first(&self) -> Option<Decision> {
        self.add_decision_tag.then_some(Decision::First)
    }

    fn phred_decision_last(&self) -> Option<Decision> {
        self.add_decision_tag.then_some(Decision::Last)
    }

    fn phred_delta(&self, abs_delta: f64) -> Option<Decision> {
        self.add_decision_tag.then(|| {
            let p = (10.0 * abs_delta / std::f64::consts::LN_10) as u32;
            Decision::ConfDelta(p.min(255) as u8)
        })
    }

    fn score_one(
        &mut self,
        state: &FragmentState<R>,
        mcfs: SmallVec<[MdCigFlags<'_>; READ_CT]>,
        aln_idx: usize,
    ) -> Result<f64> {
        let mut segment: SmallVec<[&R; READ_CT]> = SmallVec::new();
        let mut seg_mcfs: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();
        let mut mcfs_opt: SmallVec<[Option<MdCigFlags>; READ_CT]> =
            mcfs.into_iter().map(Some).collect();
        let mut dvnt: FragEvalVec<'_> = SmallVec::new();

        for idx in state.order_mates() {
            let rec = &state.get_records()[idx];
            let flags = state
                .flags(idx)
                .ok_or_else(|| anyhow!("No flags for record {idx}"))?;
            if flags.is_unmapped() {
                dvnt.push(SmallVec::new());
            } else {
                let tid = rec
                    .ref_seq_id()
                    .transpose()?
                    .ok_or_else(|| anyhow!("No reference sequence ID"))?;
                let start = rec
                    .alignment_start()
                    .transpose()?
                    .ok_or_else(|| anyhow!("No alignment start"))?
                    .get();
                let cig_len = mcfs_opt[idx]
                    .as_ref()
                    .ok_or_else(|| anyhow!("MdCigFlags missing for {idx}"))?
                    .get_cigar()
                    .len();
                let end = start + cig_len;
                let delta_vars = if let Some(store) = self.aln[aln_idx].variant_store() {
                    store.overlapping_multi(tid, start, end)
                } else {
                    SmallVec::new()
                };
                dvnt.push(delta_vars);
            }
            if !flags.is_secondary() {
                segment.push(rec);
                seg_mcfs.push(
                    mcfs_opt[idx]
                        .take()
                        .ok_or_else(|| anyhow!("MdCigFlags consumed for {idx}"))?,
                );
            } else if flags.is_last_segment() {
                break;
            }
        }

        Fragment::new(&self.penalties, segment, seg_mcfs)?
            .score(&mut self.scratch, &mut dvnt)
            .map_err(|e| {
                anyhow!(
                    "Error scoring fragment for alignment {aln_idx}: {}\n{}",
                    e,
                    state
                        .get_records()
                        .iter()
                        .map(stringify_record)
                        .collect::<Vec<_>>()
                        .join("\n")
                )
            })
    }

    pub(crate) fn print_counters(&self) {
        for i in 0..2 {
            eprintln!("hashlookup[filter:{}]: {}", i, self.branch_counters[i << 1]);
            eprintln!(
                "hashlookup[out:{}]: {}",
                i,
                self.branch_counters[1 + (i << 1)]
            );
            eprintln!("hashlookup[ambig:{}]: {}", i, self.branch_counters[16 + i]);
        }
    }
}

#[cfg(test)]
mod tests;
