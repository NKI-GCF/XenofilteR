use anyhow::{Result, ensure};
use rust_htslib::bam::record::Record;
use smallvec::{SmallVec, smallvec};

use crate::CONFIG;
use crate::aln_stream::AlnStream;
use crate::fragment::FragmentState;

// Usually two species compared so two alignments, two fragment states.

use std::sync::OnceLock;

type QnameCompareFn = fn(&AlnBuffer, &[u8]) -> Option<bool>;

static COMPARE_QNAME: OnceLock<QnameCompareFn> = OnceLock::new();
static IS_UNMAPPED_SKIPPED: OnceLock<fn(&Record) -> bool> = OnceLock::new();
static IS_SECONDARY_SKIPPED: OnceLock<fn(&Record) -> bool> = OnceLock::new();

type AlnState = (usize, FragmentState);
type AlnBuffer = SmallVec<[AlnState; 2]>;

pub struct LineByLine {
    aln: SmallVec<[AlnStream; 2]>,
    branch_counters: [u64; 32],
}

impl LineByLine {
    pub fn new(aln: SmallVec<[AlnStream; 2]>) -> Self {
        LineByLine {
            aln,
            branch_counters: [0; 32],
        }
    }

    fn write_record(&mut self, i: usize, rec: Record, best_state: Option<bool>) -> Result<()> {
        let is_unmapped_skipped = IS_UNMAPPED_SKIPPED.get_or_init(init_unmapped_skipper);
        if is_unmapped_skipped(&rec) {
            return Ok(());
        }
        match (i, best_state) {
            (i, Some(false)) => self.branch_counters[i << 1] += 1,
            (i, Some(true)) => self.branch_counters[1 + (i << 1)] += 1,
            (i, None) => self.branch_counters[16 + i] += 1,
        }
        self.aln[i].write_record(rec, best_state)
    }

    fn handle_record_is_fragment_finished(
        &mut self,
        i: usize,
        rec: Record,
        best: &mut AlnBuffer,
    ) -> Result<bool> {
        let is_secondary_skipped = IS_SECONDARY_SKIPPED.get_or_init(init_secondary_skipper);
        if is_secondary_skipped(&rec) {
            let compare_qname = COMPARE_QNAME.get_or_init(init_qname_comparer);
            if let Some(new_readname) = compare_qname(best, rec.qname()) {
                if new_readname {
                    // end of round
                    self.aln[i].un_next(rec);
                    if best.len() > 1 && best.last().unwrap().1 != best[0].1 {
                        let mut last = best.pop().unwrap();
                        if last.1 > best[0].1 {
                            // TODO: use best.splice(.., last) when stable
                            best.drain(..).try_for_each(|mut b| {
                                b.1.records
                                    .drain(..)
                                    .try_for_each(|r| self.write_record(b.0, r, None))
                            })?;
                            best.push(last);
                        } else {
                            last.1
                                .records
                                .drain(..)
                                .try_for_each(|r| self.write_record(i, r, None))?;
                        }
                    }
                    return Ok(true);
                }
                for (j, state) in best.iter_mut().rev() {
                    if *j == i {
                        state.records.push(rec);
                        return Ok(false);
                    }
                }
            }
            best.push((i, FragmentState::from_record(rec)));
        } // else skip secondary

        Ok(false)
    }
    fn handle_best(&mut self, best: &mut AlnBuffer) -> Result<()> {
        if best.len() > 1 {
            best.sort_unstable_by(|a, b| {
                b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal)
            });
            if best[0].1 > best[1].1 {
                best.drain(1..).try_for_each(|mut b| {
                    b.1.records
                        .drain(..)
                        .try_for_each(|r| self.write_record(b.0, r, None))
                })?;
            }
        }
        let best_state = (best.len() == 1).then_some(true);

        best.drain(..).try_for_each(|mut b| {
            b.1.records
                .drain(..)
                .try_for_each(|r| self.write_record(b.0, r, best_state))
        })
    }

    pub fn process(&mut self) -> Result<()> {
        let mut best: AlnBuffer = smallvec![];

        let mut i = 0;
        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.handle_record_is_fragment_finished(i, rec, &mut best)? {
                    i += 1;
                    if i == self.aln.len() {
                        self.handle_best(&mut best)?;
                        i = 0;
                    }
                }
            }
            i += 1;
        }
        self.handle_best(&mut best)?;
        while i > 0 {
            i -= 1;
            eprintln!(
                "Filtered from alignment {i}: {}",
                self.branch_counters[i << 1]
            );
            eprintln!(
                "Assigned to alignment {i}: {}",
                self.branch_counters[1 + (i << 1)]
            );
            eprintln!(
                "Ambiguous for alignment {i}: {}",
                self.branch_counters[16 + i]
            );
            ensure!(
                self.aln[i].next_rec()?.is_none(),
                "alignment {i} still has reads"
            );
        }
        Ok(())
    }
}

fn init_secondary_skipper() -> fn(&Record) -> bool {
    match CONFIG.get().unwrap().skip_secondary {
        true => |rec: &Record| rec.is_secondary(),
        false => |_| false,
    }
}

fn init_unmapped_skipper() -> fn(&Record) -> bool {
    match CONFIG.get().unwrap().discard_unmapped {
        true => |rec: &Record| rec.is_unmapped() && rec.is_mate_unmapped(),
        false => |_| false,
    }
}

fn init_qname_comparer() -> fn(&AlnBuffer, &[u8]) -> Option<bool> {
    match CONFIG.get().unwrap().strip_read_suffix {
        Some(true) => |best: &AlnBuffer, qname2: &[u8]| {
            best.first()
                .map(|b| b.1.first_qname())
                .map(|qname1| qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2])
        },
        Some(false) => |best: &AlnBuffer, qname2: &[u8]| {
            best.first()
                .map(|b| b.1.first_qname())
                .map(|qname1| qname1 != qname2)
        },
        // This variant is currently unreachable, as the config enforces Some(true/false)
        // but could be useful if we want to allow input from mixed sources in the future.
        None => |best: &AlnBuffer, qname2: &[u8]| {
            best.first().map(|b| b.1.first_qname()).map(|qname1| {
                if qname1.ends_with(b"/1") || qname1.ends_with(b"/2") {
                    qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2]
                } else {
                    qname1 != qname2
                }
            })
        },
    }
}
