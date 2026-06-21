//! Output staging for the hash-lookup backend.
//!
//! Because fragments complete out-of-order but the driving stream's record
//! order must be preserved in output, completed scored fragments are staged in
//! a `BTreeMap` keyed by their driving-stream sequence number. The staging
//! buffer is flushed whenever the minimum key equals the next expected
//! sequence number.

use crate::alignment::SimpleRec;
use crate::aln_stream::AlignmentStream;
use crate::filter_algorithm::hash_lookup::ScoredFragment;
use anyhow::Result;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::RecordBuf;
use crate::filter_algorithm::line_by_line::ordering::Decision;
use std::collections::BTreeMap;

pub(crate) struct StagedOutput<R: SimpleRec> {
    next_emit: u64,
    pending: BTreeMap<u64, ScoredFragment<R>>,
}

impl<R: SimpleRec> StagedOutput<R> {
    pub(crate) fn new() -> Self {
        Self { next_emit: 0, pending: BTreeMap::new() }
    }

    /// Stage a completed fragment. Call `flush` afterwards.
    pub(crate) fn push(&mut self, seq_nr: u64, frag: ScoredFragment<R>) {
        self.pending.insert(seq_nr, frag);
    }

    /// Flush as many consecutive fragments as possible in emission order.
    pub(crate) fn flush(
        &mut self,
        aln: &mut [Box<dyn AlignmentStream<R>>],
        branch_counters: &mut [u64; 32],
        add_decision_tag: bool,
    ) -> Result<()> {
        while let Some(&min_key) = self.pending.keys().next() {
            if min_key != self.next_emit {
                break;
            }
            let sf = self.pending.remove(&min_key).unwrap();
            emit_scored(sf, aln, branch_counters, add_decision_tag)?;
            self.next_emit += 1;
        }
        Ok(())
    }

    /// Flush everything remaining regardless of order (called at end-of-stream).
    pub(crate) fn flush_all(
        &mut self,
        aln: &mut [Box<dyn AlignmentStream<R>>],
        branch_counters: &mut [u64; 32],
        add_decision_tag: bool,
    ) -> Result<()> {
        let keys: Vec<u64> = self.pending.keys().copied().collect();
        for k in keys {
            let sf = self.pending.remove(&k).unwrap();
            emit_scored(sf, aln, branch_counters, add_decision_tag)?;
        }
        Ok(())
    }
}

fn emit_scored<R: SimpleRec>(
    sf: ScoredFragment<R>,
    aln: &mut [Box<dyn AlignmentStream<R>>],
    branch_counters: &mut [u64; 32],
    add_decision_tag: bool,
) -> Result<()> {
    emit_frag(sf.winner, true, sf.decision.as_ref(), aln, branch_counters, add_decision_tag)?;
    if let Some(loser) = sf.loser {
        emit_frag(loser, false, None, aln, branch_counters, add_decision_tag)?;
    }
    Ok(())
}

fn emit_frag<R: SimpleRec>(
    mut frag: crate::alignment::FragmentState<R>,
    is_winner: bool,
    decision: Option<&Decision>,
    aln: &mut [Box<dyn AlignmentStream<R>>],
    branch_counters: &mut [u64; 32],
    add_decision_tag: bool,
) -> Result<()> {
    let nr = frag.get_nr();
    let best_state = if is_winner { Some(true) } else { Some(false) };
    let is_ambiguous = decision.map(|d| matches!(d, Decision::Ambiguous)).unwrap_or(false);
    let effective_state = if is_ambiguous { None } else { best_state };

    for r in frag.drain_records() {
        let header = aln[nr].header();
        let mut rec: RecordBuf = r.as_record_buf(header)?;
        if is_winner && !is_ambiguous {
            match decision {
                Some(Decision::ConfDelta(v)) => {
                    rec.data_mut().insert(Tag::new(b'X', b'F'), Value::from(*v));
                }
                Some(Decision::VariantRescued(v)) => {
                    rec.data_mut().insert(Tag::new(b'X', b'R'), Value::from(*v));
                }
                _ => {}
            }
            branch_counters[1 + (nr << 1)] += 1;
        } else if is_ambiguous {
            branch_counters[16 + nr] += 1;
        } else {
            branch_counters[nr << 1] += 1;
        }
        aln[nr].write_record(rec, effective_state)?;
    }
    Ok(())
}
