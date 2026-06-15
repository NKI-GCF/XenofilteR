use crate::aln_stream::AlignmentStream;
use crate::alignment::FragmentState;
use crate::{config::{Config, StripReadSuffix}, penalty::Penalty};
use anyhow::Result;
use noodles::sam::alignment::Record;
use smallvec::SmallVec;
use crate::alignment::SimpleRec;

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

pub(crate) struct LineByLine<R> {
    pub(super) aln: SmallVec<[Box<dyn AlignmentStream<R>>; 2]>,
    pub(super) branch_counters: [u64; 32],
    pub(super) is_secondary_skipped: RecordEvalFn,
    pub(super) is_unmapped_skipped: RecordEvalFn,
    pub(super) is_new_qname: fn(&AlnBuffer<R>, &[u8]) -> Option<bool>,
    pub(super) add_decision_tag: bool,
    pub(super) penalties: Penalty,
    pub(super) ambiguous_log_threshold: f64,
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
            t => (t as f64) * std::f64::consts::LN_10 / 10.0, // convert from phred to probability
        };

        let is_new_qname = match config.strip_read_suffix {
            StripReadSuffix::True => |best: &AlnBuffer<R>, qname2: &[u8]| {
                best.first()
                    .map(|b| b.first_qname())
                    .map(|qname1| qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2])
            },
            StripReadSuffix::False => |best: &AlnBuffer<R>, qname2: &[u8]| {
                best.first()
                    .map(|b| b.first_qname())
                    .map(|qname1| qname1 != qname2)
            },
            StripReadSuffix::Variable => |best: &AlnBuffer<R>, qname2: &[u8]| {
                best.first().map(|b| b.first_qname()).map(|qname1| {
                    if qname1.ends_with(b"/1") || qname1.ends_with(b"/2") {
                        qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2]
                    } else {
                        qname1 != qname2
                    }
                })
            },
            StripReadSuffix::Auto => {
                #[cfg(not(test))]
                unreachable!("Auto mode should be handled during AlnStream initialization");
                #[cfg(test)]
                debug_new_qname_fn()
            }
        };
        for i in 0..aln.len() {
            if let Some(a) = aln.get_mut(i) {
                a.init_writers(&config, i)?;
            }
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
        })
    }
}

#[cfg(test)]
fn debug_new_qname_fn<R: SimpleRec>() -> fn(&AlnBuffer<R>, &[u8]) -> Option<bool> {
    |best: &AlnBuffer<R>, qname2: &[u8]| {
        if let Some(first_qname) = best.first().map(|b| b.first_qname()) {
            if first_qname.ends_with(b"/1") || first_qname.ends_with(b"/2") {
                return best.first().map(|b| b.first_qname()).map(|qname1| {
                    qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2]
                });
            }
        }
        best.first()
            .map(|b| b.first_qname())
            .map(|qname1| qname1 != qname2)
    }
}

