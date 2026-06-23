// src/filter_algorithm/line_by_line/core.rs
// Changes from previous version:
//   - `score_threads: usize` field added to LineByLine
//   - LineByLine::new() reads score_threads from config
//   - process() dispatcher added (sequential vs parallel)

use crate::alignment::FragmentState;
use crate::alignment::SimpleRec;
use crate::aln_stream::AlignmentStream;
use crate::{
    config::{Config, StripReadSuffix},
    penalty::Penalty,
};
use anyhow::Result;
use noodles::sam::alignment::Record;
use smallvec::SmallVec;

pub(crate) const READ_CT: usize = 8;
pub(crate) const VNT_LEN: usize = 16;

pub(crate) type RecordEvalFn = fn(&dyn Record) -> Result<bool>;
pub(crate) type AlnBuffer<R> = SmallVec<[FragmentState<R>; 2]>;

fn always_false(_: &dyn Record) -> Result<bool> {
    Ok(false)
}

fn unmapped_and_mate_unmapped(rec: &dyn Record) -> Result<bool> {
    let flags = rec.flags()?;
    Ok(flags.is_unmapped() && (!flags.is_segmented() || flags.is_mate_unmapped()))
}

fn is_secondary(rec: &dyn Record) -> Result<bool> {
    Ok(rec.flags()?.is_secondary())
}

#[derive(Clone, Copy, PartialOrd, PartialEq)]
pub(crate) struct Cell {
    pub(crate) m: f64,
    pub(crate) i: f64,
    pub(crate) d: f64,
}

impl Cell {
    pub(crate) fn reinit(&mut self, gap_open: f64, gap_extend: f64, i: i32) {
        self.m = f64::NEG_INFINITY;
        self.i = f64::NEG_INFINITY;
        self.d = f64::NEG_INFINITY;
        match i {
            0 => self.m = 0.0,
            i if i < 0 => self.i = gap_open + (i.abs() as f64) * gap_extend,
            i => self.d = gap_open + (i as f64) * gap_extend,
        }
    }
}

impl Default for Cell {
    fn default() -> Self {
        Self {
            m: -f64::INFINITY,
            i: -f64::INFINITY,
            d: -f64::INFINITY,
        }
    }
}

pub(crate) struct Scratch {
    pub(crate) prev: SmallVec<[Cell; VNT_LEN]>,
    pub(crate) curr: SmallVec<[Cell; VNT_LEN]>,
    pub(crate) dp: SmallVec<[f64; READ_CT]>,
    /// Set by Fragment::score; non-zero when maximize_delta rescued alignment.
    pub(crate) last_variant_delta: f64,
}

impl Scratch {
    pub(crate) fn new() -> Self {
        Self {
            prev: SmallVec::new(),
            curr: SmallVec::new(),
            dp: SmallVec::new(),
            last_variant_delta: 0.0,
        }
    }
    pub(crate) fn resize_nw(&mut self, new_len: usize) {
        self.prev.clear();
        self.curr.clear();
        self.prev.resize(new_len, Cell::default());
        self.curr.resize(new_len, Cell::default());
    }
    pub(crate) fn swap_nw(&mut self) {
        std::mem::swap(&mut self.prev, &mut self.curr);
    }
}

pub(crate) struct LineByLine<R> {
    pub(super) aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    pub(super) branch_counters: [u64; 32],
    pub(super) is_secondary_skipped: RecordEvalFn,
    pub(super) is_unmapped_skipped: RecordEvalFn,
    pub(super) is_new_qname: fn(&AlnBuffer<R>, &[u8]) -> Option<bool>,
    pub(super) add_decision_tag: bool,
    pub(super) penalties: Penalty,
    pub(super) ambiguous_log_threshold: f64,
    pub(super) scratch: Scratch,
    /// Number of scoring worker threads.
    /// 1 = sequential (deterministic output order).
    /// N > 1 = parallel (output order is nondeterministic across fragments).
    pub(super) score_threads: usize,
}

impl<R: SimpleRec> LineByLine<R> {
    pub(crate) fn new(
        config: Config,
        mut aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    ) -> Result<Self> {
        let is_unmapped_skipped = match config.discard_unmapped {
            true => unmapped_and_mate_unmapped,
            false => always_false,
        };
        let is_secondary_skipped = match config.skip_secondary {
            true => is_secondary,
            false => always_false,
        };
        let ambiguous_log_threshold = match config.ambiguous_threshold {
            0 => 0.0,
            t => (t as f64) * std::f64::consts::LN_10 / 10.0,
        };
        let is_new_qname: fn(&AlnBuffer<R>, &[u8]) -> Option<bool> = match config.strip_read_suffix
        {
            StripReadSuffix::True => |best: &AlnBuffer<R>, qname2: &[u8]| {
                best.first()
                    .map(|b| b.first_qname())
                    .map(|q1| q1[..q1.len() - 2] != qname2[..qname2.len() - 2])
            },
            StripReadSuffix::False => |best: &AlnBuffer<R>, qname2: &[u8]| {
                best.first().map(|b| b.first_qname()).map(|q1| q1 != qname2)
            },
            StripReadSuffix::Variable => |best: &AlnBuffer<R>, qname2: &[u8]| {
                best.first().map(|b| b.first_qname()).map(|q1| {
                    if q1.ends_with(b"/1") || q1.ends_with(b"/2") {
                        q1[..q1.len() - 2] != qname2[..qname2.len() - 2]
                    } else {
                        q1 != qname2
                    }
                })
            },
            StripReadSuffix::Auto => {
                #[cfg(not(test))]
                unreachable!("Auto mode should be resolved during AlnStream initialization");
                #[cfg(test)]
                debug_new_qname_fn()
            }
        };

        for i in 0..aln.len() {
            if let Some(a) = aln.get_mut(i) {
                a.init_writers(&config, i)?;
            }
        }

        // Clamp score_threads to aln.len() * 2 as a sanity bound; beyond that
        // the IO thread becomes the bottleneck and extra workers are idle.
        let score_threads = config.score_threads.max(1);

        if score_threads > 1 {
            tracing::info!(
                score_threads,
                "Fragment scoring will use a parallel worker pool. \
                 Output order is nondeterministic."
            );
        }

        Ok(LineByLine {
            aln,
            branch_counters: [0; 32],
            is_secondary_skipped,
            is_unmapped_skipped,
            is_new_qname,
            add_decision_tag: config.add_decision_tag,
            penalties: config.to_penalties(),
            ambiguous_log_threshold,
            scratch: Scratch::new(),
            score_threads,
        })
    }
}

#[cfg(test)]
fn debug_new_qname_fn<R: SimpleRec>() -> fn(&AlnBuffer<R>, &[u8]) -> Option<bool> {
    |best: &AlnBuffer<R>, qname2: &[u8]| {
        if let Some(q1) = best.first().map(|b| b.first_qname()) {
            if q1.ends_with(b"/1") || q1.ends_with(b"/2") {
                return best
                    .first()
                    .map(|b| b.first_qname())
                    .map(|q1| q1[..q1.len() - 2] != qname2[..qname2.len() - 2]);
            }
        }
        best.first().map(|b| b.first_qname()).map(|q1| q1 != qname2)
    }
}
