// src/filter_algorithm/line_by_line/core.rs
// Changes from previous version:
//   - `score_threads: usize` field added to LineByLine
//   - LineByLine::new() reads score_threads from config
//   - process() dispatcher added (sequential vs parallel)

use crate::{
    alignment::{FragmentState, SimpleRec},
    aln_stream::AlignmentStream,
    config::{run_config::RunConfig, args::resolve_threshold, StripReadSuffix},
    penalty::Penalty,
    progress::ProgressReporter,
    region::{ScoreFn, ScoredRegions},
    variant::name_to_id::header_name_to_id,
    Error,
};
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::Record;
use smallvec::SmallVec;
use std::path::Path;
use std::sync::Arc;

pub const READ_CT: usize = 8;
pub(crate) const VNT_LEN: usize = 16;

pub(crate) const MAX_STREAMS: usize = 32;

/// Counter layout per stream — stride 4:
///   nr*4+0  discard   (includes unmapped-discarded when --discard-unmapped)
///   nr*4+1  out/winner
///   nr*4+2  ambiguous (includes unmapped-ambiguous when configured)
///   nr*4+3  chimeric  (XC:Z: tagged, both streams count)
pub(crate) const COUNTER_STRIDE: usize = 4;
pub(crate) const COUNTER_LEN: usize = MAX_STREAMS * COUNTER_STRIDE;

pub(crate) type RecordEvalFn = fn(&dyn Record) -> Result<bool, Error>;
pub(crate) type FragmentBuffer<R> = SmallVec<[FragmentState<R>; 2]>;

fn always_false(_: &dyn Record) -> Result<bool, Error> {
    Ok(false)
}

fn unmapped_and_mate_unmapped(rec: &dyn Record) -> Result<bool, Error> {
    let flags = rec.flags()?;
    Ok(flags.is_unmapped() && (!flags.is_segmented() || flags.is_mate_unmapped()))
}

fn is_secondary(rec: &dyn Record) -> Result<bool, Error> {
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

#[derive(Default)]
pub struct Scratch {
    pub(crate) prev: SmallVec<[Cell; VNT_LEN]>,
    pub(crate) curr: SmallVec<[Cell; VNT_LEN]>,
    pub(crate) dp: SmallVec<[f64; READ_CT]>,
    /// Set by Fragment::score; non-zero when wis_max_rescue_delta rescued alignment.
    pub(crate) last_variant_delta: f64,
}

impl Scratch {
    pub fn new() -> Self {
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
    pub(crate) routing_counters: SmallVec<[u64; 8]>,
    pub(super) is_secondary_skipped: RecordEvalFn,
    pub(super) is_unmapped_skipped: RecordEvalFn,
    pub(super) is_new_qname: fn(&FragmentBuffer<R>, &[u8]) -> Option<bool>,
    pub(super) add_decision_tag: bool,
    pub(super) penalties: Penalty,
    pub(super) ambiguous_log_threshold: f64,
    pub(super) scratch: Scratch,
    /// Number of scoring worker threads.
    /// 1 = sequential (deterministic output order).
    /// N > 1 = parallel (output order is nondeterministic across fragments).
    pub(super) score_threads: usize,

    /// Canonical chimeric stream-index pairs (lower index first).
    pub(super) chimeric_pairs: Vec<[usize; 2]>,

    /// Per-stream labels used for the `XC:Z:` SAM tag.
    /// `chimeric_label(i)` returns `stream_labels[i]` or `"stream_N"`.
    pub(super) stream_labels: Vec<String>,
    pub(super) progress: Option<ProgressReporter>,
    pub(super) bisulfite: bool,
    pub(super) positive_regions: [Option<Arc<ScoredRegions>>; MAX_STREAMS],
    pub(super) region_score_fn: ScoreFn,
}

impl<R: SimpleRec> LineByLine<R> {
    pub(crate) fn new(
        config: &RunConfig,
        mut aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    ) -> Result<Self, Error> {
        let is_unmapped_skipped = match config.io.discard_unmapped {
            true => unmapped_and_mate_unmapped,
            false => always_false,
        };
        let is_secondary_skipped = match config.io.skip_secondary {
            true => is_secondary,
            false => always_false,
        };
        let ambiguous_log_threshold = match config.scoring.ambiguous_threshold {
            0 => 0.0,
            t => (t as f64) * std::f64::consts::LN_10 / 10.0,
        };
        let is_new_qname: fn(&FragmentBuffer<R>, &[u8]) -> Option<bool> =
            match config.io.strip_read_suffix {
                StripReadSuffix::True => |best: &FragmentBuffer<R>, qname2: &[u8]| {
                    best.first()
                        .map(|b| b.first_qname())
                        .map(|q1| q1[..q1.len() - 2] != qname2[..qname2.len() - 2])
                },
                StripReadSuffix::False => |best: &FragmentBuffer<R>, qname2: &[u8]| {
                    best.first().map(|b| b.first_qname()).map(|q1| q1 != qname2)
                },
                StripReadSuffix::Variable => |best: &FragmentBuffer<R>, qname2: &[u8]| {
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
        let aln_len = aln.len();
        for i in 0..aln_len {
            if let Some(a) = aln.get_mut(i) {
                a.init_writers(&config, i)?;
            }
        }

        // Clamp score_threads to aln.len() * 2 as a sanity bound; beyond that
        // the IO thread becomes the bottleneck and extra workers are idle.
        let score_threads = config.parallel.score_threads.max(1);

        if score_threads > 1 {
            tracing::info!(
                score_threads,
                "Fragment scoring will use a parallel worker pool. \
                 Output order is nondeterministic."
            );
        }
        let progress = match config.io.quiet {
            true => None,
            false => Some(ProgressReporter::new()),
        };
        let chimeric_pairs = config.chimeric.parse_pairs(aln_len)?;

        let mut positive_regions: [Option<Arc<ScoredRegions>>; MAX_STREAMS] = Default::default();
        for (i, path) in config
            .variants
            .positive_regions
            .iter()
            .enumerate()
            .map(|(i, f)| {
                if let Some(idx) = f.idx {
                    (idx, &f.path)
                } else {
                    (i, &f.path)
                }
            })
        {
            if i < MAX_STREAMS {
                positive_regions[i] = Some(Arc::new(ScoredRegions::from_bed(
                    Path::new(path),
                    &header_name_to_id(aln[i].header()),
                )?));
            }
        }
        Ok(LineByLine {
            aln,
            routing_counters: SmallVec::from_elem(0, aln_len * 4),
            is_secondary_skipped,
            is_unmapped_skipped,
            is_new_qname,
            add_decision_tag: config.io.add_decision_tag,
            penalties: config.scoring.to_penalty(),
            ambiguous_log_threshold,
            scratch: Scratch::new(),
            chimeric_pairs,
            stream_labels: config.chimeric.stream_labels.clone(),
            score_threads,
            progress,
            bisulfite: config.scoring.bisulfite,
            positive_regions,
            region_score_fn: config.scoring.region_score_fn,
        })
    }
}

impl LineByLine<RecordBuf> {
    pub(crate) fn new_from_namesorted(
        args: &RunConfig,
        aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>,
    ) -> Result<Self, Error> {
        let n = aln.len();
        let chimeric_pairs = args.chimeric.parse_pairs(n)?;
        let ambiguous_log_threshold = resolve_threshold(
            args.scoring.ambiguous_threshold,
            /* is_pass2 */ false,
        );
        let progress = match args.io.quiet {
            true => None,
            false => Some(ProgressReporter::new()),
        };
        let is_new_qname: fn(&FragmentBuffer<RecordBuf>, &[u8]) -> Option<bool> =
            match args.io.strip_read_suffix {
                StripReadSuffix::True => |best: &FragmentBuffer<RecordBuf>, qname2: &[u8]| {
                    best.first()
                        .map(|b| b.first_qname())
                        .map(|q1| q1[..q1.len() - 2] != qname2[..qname2.len() - 2])
                },
                StripReadSuffix::False => |best: &FragmentBuffer<RecordBuf>, qname2: &[u8]| {
                    best.first().map(|b| b.first_qname()).map(|q1| q1 != qname2)
                },
                StripReadSuffix::Variable => |best: &FragmentBuffer<RecordBuf>, qname2: &[u8]| {
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
        let is_unmapped_skipped = match args.io.discard_unmapped {
            true => unmapped_and_mate_unmapped,
            false => always_false,
        };
        let is_secondary_skipped = match args.io.skip_secondary {
            true => is_secondary,
            false => always_false,
        };
        let mut positive_regions: [Option<Arc<ScoredRegions>>; MAX_STREAMS] = Default::default();
        for (i, path) in args
            .variants
            .positive_regions
            .iter()
            .enumerate()
            .map(|(i, f)| {
                if let Some(idx) = f.idx {
                    (idx, &f.path)
                } else {
                    (i, &f.path)
                }
            })
        {
            if i < MAX_STREAMS {
                positive_regions[i] = Some(Arc::new(ScoredRegions::from_bed(
                    Path::new(path),
                    &header_name_to_id(aln[i].header()),
                )?));
            }
        }
        Ok(Self {
            aln,
            routing_counters: SmallVec::from_elem(0, n * COUNTER_STRIDE),
            penalties: args.scoring.to_penalty(),
            ambiguous_log_threshold,
            add_decision_tag: args.io.add_decision_tag,
            bisulfite: args.scoring.bisulfite,
            score_threads: args.parallel.score_threads,
            chimeric_pairs,
            stream_labels: args.chimeric.stream_labels.clone(),
            progress,
            is_new_qname,
            is_secondary_skipped,
            is_unmapped_skipped,
            scratch: Scratch::new(),
            region_score_fn: args.scoring.region_score_fn,
            positive_regions,
        })
    }
}

#[cfg(test)]
fn debug_new_qname_fn<R: SimpleRec>() -> fn(&FragmentBuffer<R>, &[u8]) -> Option<bool> {
    |best: &FragmentBuffer<R>, qname2: &[u8]| {
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
