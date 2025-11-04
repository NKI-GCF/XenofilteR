use anyhow::{Result, ensure};
use rust_htslib::bam::record::Record;
use smallvec::{SmallVec, smallvec};

use crate::aln_stream::AlnStream;
use crate::fragment::FragmentState;
use crate::CONFIG;

// Usually two species compared so two alignments, two fragment states.

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

    fn filter_records(&mut self, i: usize, mut fs: FragmentState) -> Result<()> {
        let config = CONFIG.get().unwrap();
        for r in fs.drain() {
            if config.discard_unmapped && r.is_unmapped() && r.is_mate_unmapped() {
                continue;
            }
            self.write_record(i, r, Some(false))?;
        }
        Ok(())
    }
    fn neqn(&self, best: &AlnBuffer, qname2: &[u8]) -> Option<bool> {
        let config = CONFIG.get().unwrap();
        best.first().map(|b| b.1.first_qname()).map(|qname1| {
            if config.strip_read_suffix.unwrap_or(false) {
                qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2]
            } else {
                qname1 != qname2
            }
        })
    }

    fn write_record(&mut self, i: usize, rec: Record, best_state: Option<bool>) -> Result<()> {
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
    ) -> bool {
        let config = CONFIG.get().unwrap();
        if !config.skip_secondary || !rec.is_secondary() {
            if let Some(new_readname) = self.neqn(best, rec.qname()) {
                if new_readname {
                    // end of round
                    self.aln[i].un_next(rec);
                    if best.len() > 1 && best.last().unwrap().1 != best[0].1 {
                        let last = best.pop().unwrap();
                        if last.1 > best[0].1 {
                            for b in best.drain(..) {
                                self.filter_records(b.0, b.1).unwrap();
                            }
                            best.push(last);
                        } else {
                            self.filter_records(i, last.1).unwrap();
                        }
                    }
                    return true;
                }
                for (j, state) in best.iter_mut().rev() {
                    if *j == i {
                        state.add_record(rec);
                        return false;
                    }
                }
            }
            best.push((
                i,
                FragmentState::from_record(rec),
            ));
        } // else skip secondary

        false
    }
    fn handle_best(&mut self, best: &mut AlnBuffer) -> Result<()> {
        if best.len() > 1 {
            best.sort_unstable_by(|a, b| {
                b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal)
            });
            if best[0].1 > best[1].1 {
                for (i, state) in best.drain(1..) {
                    self.filter_records(i, state)?;
                }
            }
        }
        let best_state = (best.len() == 1).then_some(true);

        for (i, mut state) in best.drain(..) {
            for rec in state.drain() {
                self.write_record(i, rec, best_state)?;
            }
        }
        Ok(())
    }

    pub fn process(&mut self) -> Result<()> {
        let mut best: AlnBuffer = smallvec![];

        let mut i = 0;
        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                if self.handle_record_is_fragment_finished(i, rec, &mut best) {
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
