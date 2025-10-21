use anyhow::{Result, ensure};
use smallvec::{SmallVec, smallvec};
use rust_htslib::bam::record::Record;

use crate::aln_stream::AlnStream;
use crate::Config;
use crate::fragment::FragmentState;

type AlnState = (usize, FragmentState);
type AlnBuffer = SmallVec<[AlnState; 2]>;

pub struct LineByLine {
    aln: SmallVec<[AlnStream; 2]>,
    log_likelihood_mismatch: [f64; 95],
    config: Config,
    branch_counters: [u64; 32]
}

impl LineByLine {
    pub fn new(config: Config, log_likelihood_mismatch: [f64; 95], aln: SmallVec<[AlnStream; 2]>) -> Result<Self> {
        Ok(LineByLine {
            aln,
            log_likelihood_mismatch,
            config,
            branch_counters: [0; 32],
        })
    }

    fn filter_records(&mut self, i: usize, mut fs: FragmentState) -> Result<()> {
        for r in fs.drain() {
            if self.config.discard_unmapped && r.is_unmapped() && r.is_mate_unmapped() {
                continue;
            }
            self.write_record(i, r, Some(false))?;
        }
        Ok(())
    }
    fn neqn(&self, best: &AlnBuffer, qname2: &[u8]) -> Option<bool> {
        best.first().map(|b| b.1.first_qname()).map(|qname1| {
            if self.config.strip_read_suffix.unwrap_or(false) {
                qname1[..qname1.len()-2] != qname2[..qname2.len()-2]
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
        };
        self.aln[i].write_record(rec, best_state)
    }

    fn handle_record_is_fragment_finished(&mut self, i: usize, rec: Record, best: &mut AlnBuffer) -> bool {

        if !self.config.skip_secondary || !rec.is_secondary() {
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
                    return true
                }
                for (j, state) in best.iter_mut().rev() {
                    if *j == i {
                        state.add_record(rec);
                        return false;
                    }
                }
            }
            best.push((i, FragmentState::from_record(rec, self.log_likelihood_mismatch)));
        } // else skip secondary

        false
    }
    fn handle_best(&mut self, best: &mut AlnBuffer) -> Result<()> {
        if best.len() > 1 {
            best.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            if best[0].1 > best[1].1 {
                best.truncate(1); // Keep only the single best alignment.
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
        loop {
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
            if i == self.aln.len() {
                self.handle_best(&mut best)?;
                break;
            }
        }
        for i in 0..self.aln.len() {
            eprintln!("Filtered from alignment {i}: {}", self.branch_counters[i << 1]);
            eprintln!("Assigned to alignment {i}: {}", self.branch_counters[1 + (i << 1)]);
            eprintln!("Ambiguous for alignment {i}: {}", self.branch_counters[16 + i]);
        }

        for i in 0..self.aln.len() {
            ensure!(self.aln[i].next_rec()?.is_none(), "alignment {i} still has reads");
        }
        Ok(())
    }
}


