use crate::aln_stream::AlignmentStream;
use crate::alignment::FragmentState;
use crate::{config::{Config, StripReadSuffix}, penalty::Penalty};
use anyhow::{Result, ensure};
use noodles::sam::alignment::Record;
use smallvec::{SmallVec, smallvec};
use noodles::bam::record::Record as BamRecord;

pub(crate) type RecordEvalFn = fn(&dyn Record) -> Result<bool>;
pub(crate) type AlnBuffer = SmallVec<[FragmentState<BamRecord>; 2]>;

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

pub(crate) struct LineByLine {
    pub(super) aln: SmallVec<[Box<dyn AlignmentStream>; 2]>,
    pub(super) branch_counters: [u64; 32],
    pub(super) is_secondary_skipped: RecordEvalFn,
    pub(super) is_unmapped_skipped: RecordEvalFn,
    pub(super) is_new_qname: fn(&AlnBuffer, &[u8]) -> Option<bool>,
    pub(super) add_decision_tag: bool,
    pub(super) penalties: Penalty,
    pub(super) ambiguous_log_threshold: f64,
}

impl LineByLine {
    pub(crate) fn new(
        config: Config,
        mut aln: SmallVec<[Box<dyn AlignmentStream>; 2]>,
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
            StripReadSuffix::True => |best: &AlnBuffer, qname2: &[u8]| {
                best.first()
                    .map(|b| b.first_qname())
                    .map(|qname1| qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2])
            },
            StripReadSuffix::False => |best: &AlnBuffer, qname2: &[u8]| {
                best.first()
                    .map(|b| b.first_qname())
                    .map(|qname1| qname1 != qname2)
            },
            StripReadSuffix::Variable => |best: &AlnBuffer, qname2: &[u8]| {
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

    pub(crate) fn process(&mut self) -> Result<()> {
        let mut best: AlnBuffer = smallvec![];

        let mut i = 0;
        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.handle_record_is_fragment_finished(i, rec, &mut best)? {
                    break;
                }
            }
            let mut decision = None;
            if best.len() > 1 {
                let last_idx = best.len() - 1;
                let mut ord = best[0].partial_cmp(&best[last_idx]);
                #[cfg(test)]
                debug_print_best(&best, &best[last_idx], ord);

                if ord.is_none() {
                    let refs_first = best[0].md_cig_refs()?;
                    let refs_last  = best[last_idx].md_cig_refs()?;

                    // Compare first read of each fragment (index 0).
                    if let (Some(a), Some(b)) =
                        (refs_first.first(), refs_last.first())
                    {
                        ord = a.partial_cmp(b);
                        #[cfg(test)]
                        debug_print_best(&best, &best[last_idx], ord);
                    }
                }

                decision = self.handle_ordering(&mut best, ord)?;
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
}

#[cfg(test)]
fn debug_print_best(best: &AlnBuffer, last: &AlnBuffer, ord: Option<std::cmp::Ordering>) {
    assert_eq!(best[0].records[0].name(), last.records[0].name());
    eprintln!(
        "{}: {} vs {} => {:?}",
        std::str::from_utf8(best[0].records[0].name().as_ref().unwrap()).unwrap_or("<?>"),
        best[0].get_nr(),
        best.last().unwrap().get_nr(),
        ord
    );
}

#[cfg(test)]
fn debug_new_qname_fn() -> fn(&AlnBuffer, &[u8]) -> Option<bool> {
    |best: &AlnBuffer, qname2: &[u8]| {
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

