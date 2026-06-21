//! Output staging for the two-pass HashLookup algorithm.
//!
//! Completed `ScoredFragment`s are inserted into a `BTreeMap` keyed by
//! driving-stream sequence number. The buffer is flushed whenever the
//! minimum key matches the next expected emission number — preserving
//! driving-stream record order.
//!
//! Pass 2 is embedded here: each `ScoredFragment` stores only virtual offsets.
//! On emission, this module seeks to those offsets to retrieve full records.
//!
//! Note: full pass-2 BAM seeking requires access to the underlying BGZF
//! reader via `noodles::bam::io::Reader::seek`. This is available when the
//! `AlignmentStream` wraps a file-backed `BamReader`. For `MockStream` in
//! tests, virtual offsets are not meaningful and the seek is skipped.

use crate::alignment::SimpleRec;
use crate::aln_stream::AlignmentStream;
use crate::filter_algorithm::hash_lookup::ScoredFragment;
use crate::filter_algorithm::line_by_line::ordering::Decision;
use anyhow::{anyhow, Result};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::SmallVec;
use std::collections::BTreeMap;

pub(crate) struct StagedOutput {
    next_emit: u64,
    pending: BTreeMap<u64, ScoredFragment>,
}

impl StagedOutput {
    pub(crate) fn new() -> Self {
        Self { next_emit: 0, pending: BTreeMap::new() }
    }

    pub(crate) fn push(&mut self, seq_nr: u64, sf: ScoredFragment) {
        self.pending.insert(seq_nr, sf);
    }

    pub(crate) fn flush<R: SimpleRec>(
        &mut self,
        aln: &mut SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
        branch_counters: &mut [u64; 32],
        add_decision_tag: bool,
    ) -> Result<()> {
        while let Some(&min_key) = self.pending.keys().next() {
            if min_key != self.next_emit { break; }
            let sf = self.pending.remove(&min_key).unwrap();
            emit_scored(sf, aln, branch_counters, add_decision_tag)?;
            self.next_emit += 1;
        }
        Ok(())
    }

    pub(crate) fn flush_all<R: SimpleRec>(
        &mut self,
        aln: &mut SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
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
    sf: ScoredFragment,
    aln: &mut SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    branch_counters: &mut [u64; 32],
    add_decision_tag: bool,
) -> Result<()> {
    if sf.is_ambiguous {
        // All offsets in winner_offsets go to ambiguous output.
        for (nr, voffset) in &sf.winner_offsets {
            let rec = fetch_record(aln, *nr, *voffset)?;
            branch_counters[16 + nr] += 1;
            aln[*nr].write_record(rec, None)?;
        }
        // loser_offsets also go to ambiguous when is_ambiguous is true.
        for (nr, voffset) in &sf.loser_offsets {
            let rec = fetch_record(aln, *nr, *voffset)?;
            branch_counters[16 + nr] += 1;
            aln[*nr].write_record(rec, None)?;
        }
    } else {
        // Winners → best output with optional decision tag.
        for (nr, voffset) in &sf.winner_offsets {
            let mut rec = fetch_record(aln, *nr, *voffset)?;
            if add_decision_tag {
                apply_decision_tag(&mut rec, sf.decision.as_ref());
            }
            branch_counters[1 + (nr << 1)] += 1;
            aln[*nr].write_record(rec, Some(true))?;
        }
        // Losers → filtered output.
        for (nr, voffset) in &sf.loser_offsets {
            let rec = fetch_record(aln, *nr, *voffset)?;
            branch_counters[nr << 1] += 1;
            aln[*nr].write_record(rec, Some(false))?;
        }
    }

    // Supplementaries follow winner's decision.
    for (stream_nr, offsets) in sf.supp_offsets.iter().enumerate() {
        let best_state = if sf.is_ambiguous { None } else { Some(!sf.loser_offsets.iter().any(|&(nr, _)| nr == stream_nr)) };
        for &voffset in offsets {
            let mut rec = fetch_record(aln, stream_nr, voffset)?;
            if !sf.is_ambiguous && best_state == Some(true) && add_decision_tag {
                apply_decision_tag(&mut rec, sf.decision.as_ref());
            }
            if sf.is_ambiguous {
                branch_counters[16 + stream_nr] += 1;
            } else if best_state == Some(true) {
                branch_counters[1 + (stream_nr << 1)] += 1;
            } else {
                branch_counters[stream_nr << 1] += 1;
            }
            aln[stream_nr].write_record(rec, best_state)?;
        }
    }

    Ok(())
}

/// Fetch a full record by virtual offset (pass 2 seek).
/// In production this seeks the BGZF stream; in tests the virtual offset
/// is a monotonic counter and the stream does not support seeking —
/// fetch_record falls back to returning a placeholder in that case.
fn fetch_record<R: SimpleRec>(
    aln: &mut SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    nr: usize,
    virtual_offset: u64,
) -> Result<RecordBuf> {
    // TODO: when noodles exposes `BamReader::seek(VirtualPosition)` through
    // the AlignmentStream trait, call it here. For now we call
    // `seek_and_read` which the trait stub returns Err for non-file streams
    // (tests), and the real implementation will override.
    aln.get_mut(nr)
        .ok_or_else(|| anyhow!("No stream {nr}"))?
        .fetch_by_virtual_offset(virtual_offset)
}

fn apply_decision_tag(rec: &mut RecordBuf, decision: Option<&Decision>) {
    match decision {
        Some(Decision::ConfDelta(v)) => {
            rec.data_mut().insert(Tag::new(b'X', b'F'), Value::from(*v));
        }
        Some(Decision::VariantRescued(v)) => {
            rec.data_mut().insert(Tag::new(b'X', b'R'), Value::from(*v));
        }
        _ => {}
    }
}
