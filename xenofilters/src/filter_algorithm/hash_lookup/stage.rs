//! Output staging: preserves driving-stream order via a `BTreeMap` keyed by
//! sequence number. Pass-2 seeks are performed here by calling
//! `fetch_by_virtual_offset` on the appropriate `AlignmentStream`.

use crate::alignment::SimpleRec;
use crate::aln_stream::AlignmentStream;
use crate::filter_algorithm::hash_lookup::ScoredFragment;
use crate::filter_algorithm::line_by_line::apply_decision_tag;
use crate::filter_algorithm::line_by_line::ordering::Decision;
use crate::Error;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::SmallVec;
use std::collections::BTreeMap;

// CONCURRENCY STUB — HashLookup pass-2 seek-IO thread
//
// In `emit_scored`, each `fetch_by_virtual_offset` call serialises on disk
// seeks.  To overlap seeks with scoring work, offload to a dedicated thread:
//
//   let (seek_tx, seek_rx) =
//       crossbeam_channel::unbounded::<(usize, u64, Sender<RecordBuf>)>();
//   thread::spawn(move || {
//       for (stream_nr, voffset, reply_tx) in seek_rx {
//           let rec = aln[stream_nr].fetch_by_virtual_offset(voffset).unwrap();
//           let _ = reply_tx.send(rec);
//       }
//   });
//
// Writers stay on the IO thread (no Mutex); the seek thread owns BAM handles.
// Bounded work-channel capacity provides back-pressure.
//
// N-STREAM NOTE: `FragmentTable` memory scales as O(in-flight × streams).
// Beyond 2 streams the table can exhaust RAM; profile before enabling.
pub(crate) struct StagedOutput {
    next_emit: u64,
    pending: BTreeMap<u64, ScoredFragment>,
}

impl StagedOutput {
    pub(crate) fn new() -> Self {
        Self {
            next_emit: 0,
            pending: BTreeMap::new(),
        }
    }

    pub(crate) fn push(&mut self, seq_nr: u64, sf: ScoredFragment) {
        self.pending.insert(seq_nr, sf);
    }

    pub(crate) fn flush<R: SimpleRec>(
        &mut self,
        aln: &mut SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
        routing_counters: &mut SmallVec<[u64; 8]>,
        add_decision_tag: bool,
    ) -> Result<(), Error> {
        while let Some(&min_key) = self.pending.keys().next() {
            if min_key != self.next_emit {
                break;
            }
            let sf = self.pending.remove(&min_key).unwrap();
            emit_scored(sf, aln, routing_counters, add_decision_tag)?;
            self.next_emit += 1;
        }
        Ok(())
    }

    pub(crate) fn flush_all<R: SimpleRec>(
        &mut self,
        aln: &mut SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
        routing_counters: &mut SmallVec<[u64; 8]>,
        add_decision_tag: bool,
    ) -> Result<(), Error> {
        let keys: Vec<u64> = self.pending.keys().copied().collect();
        for k in keys {
            let sf = self.pending.remove(&k).unwrap();
            emit_scored(sf, aln, routing_counters, add_decision_tag)?;
        }
        Ok(())
    }
}

fn emit_scored<R: SimpleRec>(
    sf: ScoredFragment,
    aln: &mut SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    routing_counters: &mut SmallVec<[u64; 8]>,
    add_decision_tag: bool,
) -> Result<(), Error> {
    if sf.is_ambiguous {
        for (nr, voffset) in sf.winner_offsets.iter().chain(sf.loser_offsets.iter()) {
            let rec = fetch(aln, *nr, *voffset)?;
            routing_counters[(nr * 4) + 2] += 1;
            aln[*nr].write_record(rec, None)?;
        }
    } else {
        for (nr, voffset) in &sf.winner_offsets {
            let mut rec = fetch(aln, *nr, *voffset)?;
            if add_decision_tag {
                apply_decision_tag(&mut rec, sf.decision.as_ref());
            }
            routing_counters[(nr * 4) + 1] += 1;
            aln[*nr].write_record(rec, Some(true))?;
        }
        for (nr, voffset) in &sf.loser_offsets {
            let rec = fetch(aln, *nr, *voffset)?;
            routing_counters[nr * 4] += 1;
            aln[*nr].write_record(rec, Some(false))?;
        }
    }

    // Supplementaries follow winner stream's decision.
    let winner_nr = sf.winner_nr;
    for (stream_nr, offsets) in sf.supp_offsets.iter().enumerate() {
        let base = stream_nr * 4;
        let is_winner_stream = stream_nr == winner_nr || sf.is_ambiguous;
        let best_state = if sf.is_ambiguous {
            None
        } else {
            Some(is_winner_stream)
        };
        for &voffset in offsets {
            let mut rec = fetch(aln, stream_nr, voffset)?;
            if is_winner_stream && !sf.is_ambiguous && add_decision_tag {
                apply_decision_tag(&mut rec, sf.decision.as_ref());
            }
            if sf.is_ambiguous {
                routing_counters[base + 2] += 1;
            } else if is_winner_stream {
                routing_counters[base + 1] += 1;
            } else {
                routing_counters[base] += 1;
            }
            aln[stream_nr].write_record(rec, best_state)?;
        }
    }
    Ok(())
}

fn fetch<R: SimpleRec>(
    aln: &mut SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    nr: usize,
    virtual_offset: u64,
) -> Result<RecordBuf, Error> {
    aln.get_mut(nr)
        .ok_or(Error::NoStream { nr })?
        .fetch_by_virtual_offset(virtual_offset)
}
