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
}

impl LineByLine {
    pub fn new(config: Config, log_likelihood_mismatch: [f64; 95], aln: SmallVec<[AlnStream; 2]>) -> Result<Self> {
        Ok(LineByLine {
            aln,
            log_likelihood_mismatch,
            config,
        })
    }

    fn filter_records(&mut self, i: usize, mut fs: FragmentState) -> Result<()> {
        for r in fs.drain() {
            if self.config.discard_unmapped && r.is_unmapped() && r.is_mate_unmapped() {
                continue;
            }
            self.aln[i].write_record(r, Some(false))?;
        }
        Ok(())
    }
    fn handle_record(&mut self, i: usize, fst_rec: Record, best: &mut AlnBuffer) -> bool {

        if !(self.config.skip_secondary && fst_rec.is_secondary()) {
            if let Some(new_readname) = best.first().map(|b| b.1.first_qname() != fst_rec.qname()) {
                if new_readname {
                    // end of round
                    self.aln[i].un_next(fst_rec);
                    return true
                }
                // same readname, add to existing best
                for (j, state) in best.iter_mut().rev() {
                    if *j == i {
                        state.add_record(fst_rec);
                        break;
                    }
                }
            } else {
                // first record
                best.push((i, FragmentState::from_record(fst_rec, self.log_likelihood_mismatch)));
            }
        }
        false
    }

    pub fn process(&mut self) -> Result<()> {

        let mut best: AlnBuffer = smallvec![];

        let mut i = 0;
        while let Some(fst_rec) = self.aln[i].next_rec()? {
            if self.handle_record(i, fst_rec, &mut best) {
                if i == 1 && best[0].1 > best[1].1 {
                    best.swap(0, 1);
                }
                i += 1;
                if i == self.aln.len() {
                    while best.len() > 1 && best.last().unwrap().1 < best[0].1 {
                        let (n, state) = best.pop().unwrap();
                        self.filter_records(n, state)?;
                    }
                    let best_state = (best.len() == 1).then_some(true);

                    for (i, mut state) in best.drain(..) {
                        for rec in state.drain() {
                            self.aln[i].write_record(rec, best_state)?;
                        }
                    }
                    i = 0;
                }
            }
        }

        for i in 0..self.aln.len() {
            ensure!(self.aln[i].next_rec()?.is_none(), "alignment {i} still has reads");
        }
        Ok(())
    }
}


