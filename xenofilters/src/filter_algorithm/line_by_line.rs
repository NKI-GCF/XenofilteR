use anyhow::{Result, bail, ensure};
use smallvec::{SmallVec, smallvec};

use crate::aln_stream::AlnStream;
use crate::Config;
use crate::fragment_state::FragmentState;

pub struct LineByLine {
    aln: SmallVec<[AlnStream; 2]>,
    config: Config,
}

impl LineByLine {
    pub fn new(config: Config, aln: SmallVec<[AlnStream; 2]>) -> Result<Self> {
        Ok(LineByLine {
            aln,
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

    pub fn process(&mut self) -> Result<()> {
        let mut best: SmallVec<[(usize, FragmentState); 1]> = smallvec![];
        
        while let Some(fst_rec) = self.aln[0].next_rec()? {
            if self.config.skip_secondary && fst_rec.is_secondary() {
                continue;
            }
            ensure!(best.is_empty(), "Internal error: best not empty at start of loop");
            best[0] = (0, FragmentState::from_record(fst_rec));

            while let Some(nxt_rec) = self.aln[0].next_rec()? {
                if self.config.skip_secondary && nxt_rec.is_secondary() {
                    continue;
                }
                if nxt_rec.qname() == best[0].1.first_qname() {
                    best[0].1.add_record(nxt_rec);
                } else {
                    self.aln[0].un_next(nxt_rec);
                    break;
                }
            }
            for i in 1..self.aln.len() {
                let mut fragment_state = if let Some(alt_rec) = self.aln[i].next_rec()? {
                    FragmentState::from_record(alt_rec)
                } else {
                    bail!("alignment {i} has no more reads");
                };
                while let Some(alt_rec) = self.aln[i].next_rec()? {
                    if self.config.skip_secondary && alt_rec.is_secondary() {
                        continue;
                    }
                    if alt_rec.qname() == best[0].1.first_qname() {
                        fragment_state.add_record(alt_rec);
                    } else {
                        self.aln[i].un_next(alt_rec);
                        break;
                    }
                }
                if fragment_state > best[0].1 {
                    for (i, state) in best.drain(..) {
                        self.filter_records(i, state)?;
                    }
                    best[0] = (i, fragment_state);
                } else if fragment_state == best[0].1 {
                    best.push((i, fragment_state));
                } else {
                    self.filter_records(i, fragment_state)?;
                }
            }

            // None state for ambiguous
            let best_state = (best.len() == 1).then_some(true);

            for (i, mut state) in best.drain(..) {
                for rec in state.drain() {
                    self.aln[i].write_record(rec, best_state)?;
                }
            }
        }
        Ok(())
    }
}


