use crate::alignment::stitched_fragment;
use crate::aln_stream::AlnStream;
use crate::fragment::FragmentState;
use crate::{Config, Penalties};
use anyhow::{Result, ensure};
use rust_htslib::bam::record::Record;
use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;

// Usually two species compared so two alignments, two fragment states.

type RecordEvalFn = fn(&Record) -> bool;
type AlnBuffer = SmallVec<[FragmentState; 2]>;

fn always_false(_: &Record) -> bool {
    false
}

fn unmapped_and_mate_unmapped(rec: &Record) -> bool {
    rec.is_unmapped() && rec.is_mate_unmapped()
}

fn is_secondary(rec: &Record) -> bool {
    rec.is_secondary()
}

pub struct LineByLine {
    aln: SmallVec<[AlnStream; 2]>,
    branch_counters: [u64; 32],
    is_secondary_skipped: RecordEvalFn,
    is_unmapped_skipped: RecordEvalFn,
    is_new_qname: fn(&AlnBuffer, &[u8]) -> Option<bool>,
    penalties: Penalties,
}

impl LineByLine {
    pub fn new(config: Config, aln: SmallVec<[AlnStream; 2]>) -> Self {
        let is_unmapped_skipped = match config.discard_unmapped {
            true => unmapped_and_mate_unmapped,
            false => always_false,
        };

        let is_secondary_skipped = match config.skip_secondary {
            true => is_secondary,
            false => always_false,
        };
        let is_new_qname = match config.strip_read_suffix {
            Some(true) => |best: &AlnBuffer, qname2: &[u8]| {
                eprintln!(
                    "Comparing {} to {} with suffix stripped",
                    String::from_utf8_lossy(best.first().unwrap().first_qname()),
                    String::from_utf8_lossy(qname2)
                );
                eprintln!("(stripped to {}) vs (stripped to {})",
                    String::from_utf8_lossy(&best.first().unwrap().first_qname()[..best.first().unwrap().first_qname().len()-2]),
                    String::from_utf8_lossy(&qname2[..qname2.len()-2])
                );
                best.first()
                    .map(|b| b.first_qname())
                    .map(|qname1| qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2])
            },
            Some(false) => |best: &AlnBuffer, qname2: &[u8]| {
                best.first()
                    .map(|b| b.first_qname())
                    .map(|qname1| qname1 != qname2)
            },
            None => |best: &AlnBuffer, qname2: &[u8]| {
                best.first().map(|b| b.first_qname()).map(|qname1| {
                    if qname1.ends_with(b"/1") || qname1.ends_with(b"/2") {
                        qname1[..qname1.len() - 2] != qname2[..qname2.len() - 2]
                    } else {
                        qname1 != qname2
                    }
                })
            },
        };
        LineByLine {
            aln,
            branch_counters: [0; 32],
            is_secondary_skipped,
            is_unmapped_skipped,
            is_new_qname,
            penalties: config.to_penalties(),
        }
    }

    fn write_record(&mut self, i: usize, rec: Record, best_state: Option<bool>) -> Result<()> {
        if (self.is_unmapped_skipped)(&rec) {
            return Ok(());
        }
        match (i, best_state) {
            (i, Some(false)) => self.branch_counters[i << 1] += 1,
            (i, Some(true)) => self.branch_counters[1 + (i << 1)] += 1,
            (i, None) => self.branch_counters[16 + i] += 1,
        }
        if let Some(aln) = self.aln.get_mut(i) {
            aln.write_record(rec, best_state)
        } else {
            Ok(())
        }
    }

    fn handle_ordering(&mut self, best: &mut AlnBuffer, ord: Option<Ordering>) -> Result<()> {
        match ord {
            Some(Ordering::Greater) => {
                let all_before_last = best.len() - 1;
                best.drain(0..all_before_last).try_for_each(|mut b| {
                    let nr = b.get_nr();
                    b.records
                        .drain(..)
                        .try_for_each(|r| self.write_record(nr, r, None))
                })
            }
            Some(Ordering::Less) => {
                let mut last = best.pop().unwrap();
                let nr = last.get_nr();
                last.records
                    .drain(..)
                    .try_for_each(|r| self.write_record(nr, r, None))
            }
            Some(Ordering::Equal) => {
                // Neither is written. (FIXME: make configurable?)
                Ok(())
            }
            None => {
                // None of the alignments were fully unmapped or perfect matches,
                // so we need to score them to find the best.
                let first = &best.first().unwrap();
                let first_ord = first.order_mates();

                let last = &best.last().unwrap();
                let last_ord = last.order_mates();

                let first_score =
                    stitched_fragment(&self.penalties, &first.records, first_ord)?.score()?;
                let last_score =
                    stitched_fragment(&self.penalties, &last.records, last_ord)?.score()?;

                let mut ord = first_score.partial_cmp(&last_score);
                if ord.is_none() {
                    ord = Some(Ordering::Equal);
                }
                self.handle_ordering(best, ord)
            }
        }
    }

    fn handle_record_is_fragment_finished(
        &mut self,
        i: usize,
        rec: Record,
        best: &mut AlnBuffer,
    ) -> Result<bool> {
        if !(self.is_secondary_skipped)(&rec) {
            if let Some(new_readname) = (self.is_new_qname)(best, rec.qname()) {
                if new_readname {
                    // end of round for this alignment
                    self.aln[i].un_next(rec);
                    // FIXME: comparing more than 2 alignments?
                    if best.len() > 1 {
                        let ord = best.last().unwrap().partial_cmp(&best[0]);
                        self.handle_ordering(best, ord)?;
                    }
                    return Ok(true);
                }
                for state in best.iter_mut().rev() {
                    if state.get_nr() == i {
                        state.records.push(rec);
                        return Ok(false);
                    }
                }
            }
            best.push(FragmentState::from_record(rec, i));
        } // else skip secondary

        Ok(false)
    }

    fn handle_best(&mut self, best: &mut AlnBuffer) -> Result<()> {
        if best.len() > 1 {
            best.sort_unstable_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
            if best[0] > best[1] {
                best.drain(1..).try_for_each(|mut b| {
                    let nr = b.get_nr();
                    b.records
                        .drain(..)
                        .try_for_each(|r| self.write_record(nr, r, None))
                })?;
            }
        }
        let best_state = (best.len() == 1).then_some(true);

        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.records
                .drain(..)
                .try_for_each(|r| self.write_record(nr, r, best_state))
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Config;
    use rust_htslib::bam::record::Record;

    fn mock_rec(qname: &[u8]) -> Record {
        let mut r = Record::new();
        r.set(qname, None, &[], &[]);
        r
    }

    #[test]
    fn test_qname_suffix_logic() {
        let mut config = Config::default();
        
        // Mode: None (Auto-detect /1 or /2)
        config.strip_read_suffix = None;
        let lbl = LineByLine::new(config.clone(), smallvec![]);
        let best: AlnBuffer = smallvec![FragmentState::from_record(mock_rec(b"read/1"), 0)];
        assert_eq!((lbl.is_new_qname)(&best, b"read/1"), Some(false));
        assert_eq!((lbl.is_new_qname)(&best, b"other/1"), Some(true));
        assert_eq!((lbl.is_new_qname)(&best, b"read/2"), Some(false));

        // Mode: Some(true) (Always strip last 2)
        config.strip_read_suffix = Some(true);
        let lbl = LineByLine::new(config.clone(), smallvec![]);
        assert_eq!((lbl.is_new_qname)(&best, b"read_suffix"), Some(true)); // "read" != "read_suff"

        // Mode: Some(false) (Exact match)
        config.strip_read_suffix = Some(false);
        let lbl = LineByLine::new(config, smallvec![]);
        assert_eq!((lbl.is_new_qname)(&best, b"read/1"), Some(false));
        assert_eq!((lbl.is_new_qname)(&best, b"read/2"), Some(true));
    }

    #[test]
    fn test_branch_counters_and_skipping() {
        let mut config = Config::default();
        config.discard_unmapped = true;
        config.skip_secondary = true;
        
        let mut lbl = LineByLine::new(config, smallvec![]);
        
        let mut unmapped = mock_rec(b"u");
        unmapped.set_unmapped();
        unmapped.set_mate_unmapped();
        
        let mut secondary = mock_rec(b"s");
        secondary.set_secondary();

        // Should return early (skipped)
        assert!(lbl.write_record(0, unmapped, Some(true)).is_ok());
        assert_eq!(lbl.branch_counters[1], 0);

        // handle_record_is_fragment_finished should skip secondary
        let mut best: AlnBuffer = smallvec![];
        let finished = lbl.handle_record_is_fragment_finished(0, secondary, &mut best).unwrap();
        assert!(!finished);
        assert!(best.is_empty());
    }

    #[test]
    fn test_handle_ordering_logic() {
        let lbl_setup = LineByLine::new(Config::default(), smallvec![]);
        // Test logic in handle_ordering requires mocked AlnStream for write_record calls.
        // Direct testing of branch_counters incrementation via write_record:
        let mut lbl = lbl_setup;
        
        lbl.write_record(0, mock_rec(b"r1"), Some(true)).unwrap();  // Assigned alignment 0
        lbl.write_record(0, mock_rec(b"r2"), Some(false)).unwrap(); // Filtered from alignment 0
        lbl.write_record(0, mock_rec(b"r3"), None).unwrap();        // Ambiguous alignment 0
        
        assert_eq!(lbl.branch_counters[1], 1);  // 1 + (0 << 1)
        assert_eq!(lbl.branch_counters[0], 1);  // (0 << 1)
        assert_eq!(lbl.branch_counters[16], 1); // 16 + 0
    }

    #[test]
    fn test_fragment_finished_transitions() {
        let lbl = LineByLine::new(Config::default(), smallvec![]);
        let mut lbl = lbl;
        let mut best: AlnBuffer = smallvec![FragmentState::from_record(mock_rec(b"R1"), 0)];
        
        // Same QName: continues fragment
        let fin = lbl.handle_record_is_fragment_finished(0, mock_rec(b"R1"), &mut best).unwrap();
        assert!(!fin);
        assert_eq!(best[0].records.len(), 2);

        // Different QName: finishes fragment
        // Note: this will attempt to call aln[i].un_next(), requiring a mock AlnStream.
    }
}

