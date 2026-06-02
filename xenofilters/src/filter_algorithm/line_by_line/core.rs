use crate::aln_stream::AlignmentStream;
use crate::fragment_state::FragmentState;
use crate::{Config, Penalty, StripReadSuffix};
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
                #[cfg(test)]
                {
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
                #[cfg(not(test))]
                unreachable!("Auto mode should be handled during AlnStream initialization")
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
                let last = best.last().unwrap();
                let ord = best[0].partial_cmp(last);

                #[cfg(test)]
                assert_eq!(best[0].records[0].name(), last.records[0].name());
                #[cfg(test)]
                eprintln!(
                    "{}: {} vs {} => {:?}",
                    std::str::from_utf8(best[0].records[0].name().as_ref().unwrap()).unwrap_or("<?>"),
                    best[0].get_nr(),
                    last.get_nr(),
                    ord
                );
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
