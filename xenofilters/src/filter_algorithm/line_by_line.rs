use crate::alignment::stitched_fragment;
use crate::aln_stream::AlignmentStream;
use crate::fragment::FragmentState;
use crate::{Config, Penalties, StripReadSuffix};
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
    rec.is_unmapped() && (rec.is_paired() == false || rec.is_mate_unmapped())
}

fn is_secondary(rec: &Record) -> bool {
    rec.is_secondary()
}

pub struct LineByLine {
    aln: SmallVec<[Box<dyn AlignmentStream>; 2]>,
    branch_counters: [u64; 32],
    is_secondary_skipped: RecordEvalFn,
    is_unmapped_skipped: RecordEvalFn,
    is_new_qname: fn(&AlnBuffer, &[u8]) -> Option<bool>,
    penalties: Penalties,
}

impl LineByLine {
    pub fn new(config: Config, aln: SmallVec<[Box<dyn AlignmentStream>; 2]>) -> Self {
        #[cfg(test)]
        eprintln!("Unmapped discard: {}, Secondary skip: {}, Strip suffix: {:?}", config.discard_unmapped, config.skip_secondary, config.strip_read_suffix);
        let is_unmapped_skipped = match config.discard_unmapped {
            true => unmapped_and_mate_unmapped,
            false => always_false,
        };

        let is_secondary_skipped = match config.skip_secondary {
            true => is_secondary,
            false => always_false,
        };
        let is_new_qname = match config.strip_read_suffix {
            StripReadSuffix::True => |best: &AlnBuffer, qname2: &[u8]| {
                /*eprintln!(
                    "Comparing {} to {} with suffix stripped",
                    String::from_utf8_lossy(best.first().unwrap().first_qname()),
                    String::from_utf8_lossy(qname2)
                );
                eprintln!(
                    "(stripped to {}) vs (stripped to {})",
                    String::from_utf8_lossy(
                        &best.first().unwrap().first_qname()
                            [..best.first().unwrap().first_qname().len() - 2]
                    ),
                    String::from_utf8_lossy(&qname2[..qname2.len() - 2])
                );*/
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
                            eprintln!(
                                "{} vs {}",
                                std::str::from_utf8(first_qname).unwrap_or("<?>"),
                                std::str::from_utf8(qname2).unwrap_or("<?>")
                            );
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
        match (i, best_state) {
            (i, Some(false)) => self.branch_counters[i << 1] += 1,
            (i, Some(true)) => self.branch_counters[1 + (i << 1)] += 1,
            (i, None) => {
                if (self.is_unmapped_skipped)(&rec) {
                    self.branch_counters[24 + i] += 1;
                    return Ok(());
                }
                self.branch_counters[16 + i] += 1;
            },
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
                        .try_for_each(|r| self.write_record(nr, r, Some(false)))
                })
            }
            Some(Ordering::Less) => {
                let mut last = best.pop().unwrap();
                let nr = last.get_nr();
                last.records
                    .drain(..)
                    .try_for_each(|r| self.write_record(nr, r, Some(false)))
            }
            Some(Ordering::Equal) => Ok(()),
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
                #[cfg(test)]
                eprintln!("Scoring to break tie among {} fragments: {} vs {}", best.len(), first_score, last_score);

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
                    self.aln[i].un_next(rec)?;
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
        /*#[cfg(test)]
        {
            eprintln!("Handling best for {} fragments", best.len());
            for b in best.iter() {
                eprintln!(
                    "Fragment from alignment {}: {} records",
                    b.get_nr(),
                    b.records.len()
                );
            }
        }*/
        /*if best.len() > 1 {
            best.sort_unstable_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
            if best[0] > best[1] {
                #[cfg(test)]
                panic!("THIS BRANCH SHOULD NOT HAPPEN");
                best.drain(1..).try_for_each(|mut b| {
                    let nr = b.get_nr();
                    b.records
                        .drain(..)
                        .try_for_each(|r| self.write_record(nr, r, Some(false)))
                })?;
            }
        }*/
        let best_state = (best.len() == 1).then_some(true);

        best.drain(..).try_for_each(|mut b| {
            let nr = b.get_nr();
            b.records
                .drain(..)
                .try_for_each(|r| self.write_record(nr, r, best_state))
        })
    }

    pub fn print_counters(&self, i: usize) {
        eprintln!(
            "[{}]: Filtered from alignment {i}: {}",
            i << 1,
            self.branch_counters[i << 1]
        );
        eprintln!(
            "[{}]: Assigned to alignment {i}: {}",
            1 + (i << 1),
            self.branch_counters[1 + (i << 1)]
        );
        eprintln!(
            "[{}]: Ambiguous for alignment {i}: {}",
            16 + i,
            self.branch_counters[16 + i]
        );
        eprintln!(
            "[{}]: Unmapped filtered for alignment {i}: {}",
            24 + i,
            self.branch_counters[24 + i]
        );
    }

    pub fn process(&mut self) -> Result<()> {
        let mut best: AlnBuffer = smallvec![];

        let mut i = 0;
        while i != self.aln.len() {
            while let Some(rec) = self.aln[i].next_rec()? {
                /*#[cfg(test)]
                eprintln!(
                    "Processing record from alignment {i}: {}, best_len: {}",
                    String::from_utf8_lossy(rec.qname()),
                    best.len()
                );*/
                if self.handle_record_is_fragment_finished(i, rec, &mut best)? {
                    break;
                }
            }
            if best.len() > 1 {
                let ord = best.last().unwrap().partial_cmp(&best[0]);
                #[cfg(test)]
                eprintln!("{} vs {} => {:?}", best.last().unwrap().get_nr(), best[0].get_nr(), ord);
                self.handle_ordering(&mut best, ord)?;
                assert!(!best.is_empty());
            }
            i += 1;
            if i == self.aln.len() {
                if best.is_empty() {
                    break;
                }
                eprintln!("Processing best buffer of size {}", best.len());
                self.handle_best(&mut best)?;
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
mod tests {
    use super::*;
    use crate::tests::{MockStream, create_record};
    use crate::{Config, StripReadSuffix};
    use anyhow::Result;

    fn setup_mock_streams() -> SmallVec<[Box<dyn AlignmentStream>; 2]> {
        let mut stream1 = MockStream::new(
            0,
            vec![
                create_record(b"R1", "10M", &[], &[], "10", false).unwrap(), // perfect => out
                create_record(b"R2", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => filtered
                create_record(b"R3", "10M", &[], &[], "10", false).unwrap(),  // perfect => out
                create_record(b"R4", "*", &[], &[], "10", false).unwrap(), // unmapped => filtered
                create_record(b"R5", "10M", &[], &[], "10", false).unwrap(), // perfect => ambiguous
                create_record(b"R6", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => ambiguous
                create_record(b"R7", "*", &[], &[], "10", false).unwrap(), // unmapped => ambiguous
                create_record(b"R8", "6M4S", &[], &[], "6", false).unwrap(), // mismatch => out
            ],
        );
        let mut stream2 = MockStream::new(
            1,
            vec![
                create_record(b"R1", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => filtered
                create_record(b"R2", "10M", &[], &[], "10", false).unwrap(),  // perfect => out
                create_record(b"R3", "*", &[], &[], "10", false).unwrap(), // unmapped => filtered
                create_record(b"R4", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => out
                create_record(b"R5", "10M", &[], &[], "10", false).unwrap(), // perfect => ambiguous
                create_record(b"R6", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => ambiguous
                create_record(b"R7", "*", &[], &[], "9", false).unwrap(), // unmapped => ambiguous
                create_record(b"R8", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => filtered
            ],
        );
        smallvec![
            Box::new(stream1) as Box<dyn AlignmentStream>,
            Box::new(stream2) as Box<dyn AlignmentStream>
        ]
    }

    fn setup_mock_streams_r4() -> SmallVec<[Box<dyn AlignmentStream>; 2]> {
        let mut stream1 = MockStream::new(
            0,
            vec![
                create_record(b"R4", "*", &[], &[], "10", false).unwrap(), // unmapped => filtered
            ],
        );
        let mut stream2 = MockStream::new(
            1,
            vec![
                create_record(b"R4", "5M5S", &[], &[], "5", false).unwrap(), // mismatch => out
            ],
        );
        smallvec![
            Box::new(stream1) as Box<dyn AlignmentStream>,
            Box::new(stream2) as Box<dyn AlignmentStream>
        ]
    }

    // %s/\vmock_rec\((b".*?")/create_record(\1, "10M", &[], &[], "10", false)?/g
    #[test]
    fn test_qname_suffix_logic() -> Result<()> {
        let mut config = Config::default();

        config.strip_read_suffix = StripReadSuffix::Auto;
        let lbl = LineByLine::new(config.clone(), smallvec![]);
        let best: AlnBuffer = smallvec![FragmentState::from_record(
            create_record(b"read/1", "10M", &[], &[], "10", false)?,
            0
        )];
        assert_eq!((lbl.is_new_qname)(&best, b"read/1"), Some(false));
        assert_eq!((lbl.is_new_qname)(&best, b"other/1"), Some(true));
        assert_eq!((lbl.is_new_qname)(&best, b"read/2"), Some(false));

        // Mode: Some(true) (Always strip last 2)
        config.strip_read_suffix = StripReadSuffix::True;
        let lbl = LineByLine::new(config.clone(), smallvec![]);
        assert_eq!((lbl.is_new_qname)(&best, b"read_suffix"), Some(true)); // "read" != "read_suff"

        // Mode: Some(false) (Exact match)
        config.strip_read_suffix = StripReadSuffix::False;
        let lbl = LineByLine::new(config, smallvec![]);
        assert_eq!((lbl.is_new_qname)(&best, b"read/1"), Some(false));
        assert_eq!((lbl.is_new_qname)(&best, b"read/2"), Some(true));
        Ok(())
    }

    #[test]
    fn test_branch_counters_and_skipping() -> Result<()> {
        let mut config = Config::default();
        config.discard_unmapped = true;
        config.skip_secondary = true;

        let mut lbl = LineByLine::new(config.clone(), smallvec![]);

        let mut unmapped_fwd = create_record(b"u", "*", &[], &[], "10", false)?;
        unmapped_fwd.set_unmapped();
        unmapped_fwd.set_paired();
        unmapped_fwd.set_mate_unmapped();

        let mut unmapped_rev = unmapped_fwd.clone();
        unmapped_rev.set_reverse();

        let mut secondary = create_record(b"s", "*", &[], &[], "10", false)?;
        secondary.set_secondary();

        let mut unmapped_single = create_record(b"u2", "*", &[], &[], "10", false)?;
        unmapped_single.set_unmapped();

        // Should return early (skipped)
        assert!(lbl.write_record(0, unmapped_fwd.clone(), None).is_ok());
        assert!(lbl.write_record(0, unmapped_rev.clone(), None).is_ok());
        assert!(lbl.write_record(0, unmapped_single, Some(false)).is_ok());
        lbl.print_counters(0);
        assert_eq!(lbl.branch_counters[24], 2); // unmapped:0: 2
        assert_eq!(lbl.branch_counters[0], 1); // filter:0:
        // handle_record_is_fragment_finished should skip secondary
        let mut best: AlnBuffer = smallvec![];
        let finished = lbl
            .handle_record_is_fragment_finished(0, secondary, &mut best)
            .unwrap();
        assert!(!finished);
        assert!(best.is_empty());

        config.discard_unmapped = false;
        let mut lbl = LineByLine::new(config, smallvec![]);
        assert!(lbl.write_record(0, unmapped_fwd, None).is_ok());
        assert!(lbl.write_record(0, unmapped_rev, None).is_ok());
        lbl.print_counters(0);
        assert_eq!(lbl.branch_counters[16], 2); // ambiguous:0: 2


        Ok(())
    }

    #[test]
    fn test_handle_ordering_logic() -> Result<()> {
        let lbl_setup = LineByLine::new(Config::default(), smallvec![]);
        // Test logic in handle_ordering requires mocked AlnStream for write_record calls.
        // Direct testing of branch_counters incrementation via write_record:
        let mut lbl = lbl_setup;

        lbl.write_record(
            0,
            create_record(b"r1", "10M", &[], &[], "10", false)?,
            Some(true),
        )
        .unwrap(); // Assigned alignment 0
        lbl.write_record(
            0,
            create_record(b"r2", "10M", &[], &[], "10", false)?,
            Some(false),
        )
        .unwrap(); // Filtered from alignment 0
        lbl.write_record(0, create_record(b"r3", "10M", &[], &[], "10", false)?, None)
            .unwrap(); // Ambiguous alignment 0

        assert_eq!(lbl.branch_counters[1], 1); // 1 + (0 << 1)
        assert_eq!(lbl.branch_counters[0], 1); // (0 << 1)
        assert_eq!(lbl.branch_counters[16], 1); // 16 + 0
        Ok(())
    }

    #[test]
    fn test_fragment_finished_transitions() -> Result<()> {
        let lbl = LineByLine::new(Config::default(), smallvec![]);
        let mut lbl = lbl;
        let mut best: AlnBuffer = smallvec![FragmentState::from_record(
            create_record(b"R1", "10M", &[], &[], "10", false)?,
            0
        )];

        // Same QName: continues fragment
        let fin = lbl
            .handle_record_is_fragment_finished(
                0,
                create_record(b"R1", "10M", &[], &[], "10", false)?,
                &mut best,
            )
            .unwrap();
        assert!(!fin);
        assert_eq!(best[0].records.len(), 2);

        // Different QName: finishes fragment
        // Note: this will attempt to call aln[i].un_next(), requiring a mock AlnStream.
        Ok(())
    }
    #[test]
    fn test_process_multi_stream_sync_r4() -> Result<()> {
        let mut config = Config::default();
        config.discard_unmapped = true;
        let mut lbl = LineByLine::new(config, setup_mock_streams_r4());

        // R4 -> stream 1 (mismatch vs unmapped)

        lbl.process()?;

        assert_eq!(lbl.branch_counters[0], 1); // filter:0: R4
        assert_eq!(lbl.branch_counters[1], 0); // out:0:
        assert_eq!(lbl.branch_counters[16], 0); // ambiguous:0:
        assert_eq!(lbl.branch_counters[24], 0); // unmapped:0:
        assert_eq!(lbl.branch_counters[2], 0); // filter:1:
        assert_eq!(lbl.branch_counters[3], 1); // out:1: R4
        assert_eq!(lbl.branch_counters[17], 0); // ambiguous:1:
        assert_eq!(lbl.branch_counters[25], 0); // unmapped:1:
        Ok(())
    }


    #[test]
    fn test_process_multi_stream_sync() -> Result<()> {
        let mut config = Config::default();
        config.discard_unmapped = true;
        let mut lbl = LineByLine::new(config, setup_mock_streams());

        // Streams now contain R1..R7
        // Expected winners:
        // R1 -> stream 0 (perfect vs mismatch)
        // R2 -> stream 1 (perfect vs mismatch)
        // R3 -> stream 0 (perfect vs unmapped)
        // R4 -> stream 1 (mismatch vs unmapped)
        // R5 -> tie perfect/perfect -> ambiguity
        // R6 -> tie mismatch/mismatch -> ambiguity
        // R7 -> unmapped/unmapped -> ambiguity
        // R8 -> stream 0 (mismatch vs more mismatches, but stream 1 filtered)

        lbl.process()?;

        assert_eq!(lbl.branch_counters[0], 2); // filter:0: R2, R4
        assert_eq!(lbl.branch_counters[1], 3); // out:0: R1, R3, R8
        assert_eq!(lbl.branch_counters[16], 2); // ambiguous:0: R5, R6
        assert_eq!(lbl.branch_counters[24], 0); // unmapped:0: R7
        assert_eq!(lbl.branch_counters[2], 3); // filter:1: R1, R3, R8
        assert_eq!(lbl.branch_counters[3], 2); // out:1: R2, R4
        assert_eq!(lbl.branch_counters[17], 2); // ambiguous:1: R5, R6
        assert_eq!(lbl.branch_counters[25], 1); // unmapped:1: R7

        Ok(())
    }

    #[test]
    fn test_handle_ordering_drain_logic() -> Result<()> {
        let mut lbl = LineByLine::new(Config::default(), setup_mock_streams());
        let mut best: AlnBuffer = smallvec![
            FragmentState::from_record(create_record(b"R1", "10M", &[], &[], "10", false)?, 0),
            FragmentState::from_record(create_record(b"R1", "5M5S", &[], &[], "10", false)?, 1),
        ];

        // stream 0 better than stream 1
        lbl.handle_ordering(&mut best, Some(Ordering::Greater))?;

        assert_eq!(best.len(), 1);
        assert_eq!(best[0].get_nr(), 0);
        assert_eq!(lbl.branch_counters[16 + 1], 1); // stream 1 filtered

        Ok(())
    }

    #[test]
    fn test_complex_fragment_grouping() -> Result<()> {
        let mut lbl = LineByLine::new(Config::default(), setup_mock_streams());
        let mut best: AlnBuffer = smallvec![];

        // paired-end style: same QNAME twice
        let r1_1 = create_record(b"READX", "10M", &[], &[], "10", false)?;
        let r1_2 = create_record(b"READX", "5M5S", &[], &[], "10", false)?;

        lbl.handle_record_is_fragment_finished(0, r1_1, &mut best)?;
        assert_eq!(best.len(), 1);
        assert_eq!(best[0].records.len(), 1);

        lbl.handle_record_is_fragment_finished(0, r1_2, &mut best)?;
        assert_eq!(best.len(), 1);
        assert_eq!(best[0].records.len(), 2);

        Ok(())
    }

    #[test]
    fn test_line_by_line_full_flow() -> Result<()> {
        let mut rec1 = Record::new();
        rec1.set(b"R1", None, &[], &[]);
        let mut rec2 = Record::new();
        rec2.set(b"R2", None, &[], &[]);

        // Mocking AlnStream behavior
        // Note: You may need to wrap MockStream in AlnStream enum/trait if required by your types
        // This targets handle_record_is_fragment_finished coverage
        let config = Config::default();
        let mut lbl = LineByLine::new(config, smallvec![]);

        let mut best: AlnBuffer = smallvec![];

        // 1. First read
        let fin = lbl.handle_record_is_fragment_finished(0, rec1.clone(), &mut best)?;
        assert!(!fin);
        assert_eq!(best.len(), 1);

        // 2. QName change triggers finish
        let fin = lbl.handle_record_is_fragment_finished(0, rec2, &mut best)?;
        // This will attempt to call aln[0].un_next(), ensure your mock handles this.
        assert!(fin);
        Ok(())
    }

    #[test]
    fn test_scoring_path_coverage() -> Result<()> {
        let config = Config::default();
        // Mock stream needs to exist to avoid indexing panics
        let mut lbl = LineByLine::new(config, smallvec![]);

        // Populate records with valid CIGAR/MD data to avoid null-pointer panics
        let mut best: AlnBuffer = smallvec![
            FragmentState::from_record(create_record(b"R1", "10M", &[], &[], "10", false)?, 0),
            FragmentState::from_record(create_record(b"R1", "10M", &[], &[], "10", false)?, 1),
        ];

        // This triggers the full scoring pipeline:
        // handle_ordering -> stitched_fragment -> UnifiedOpIterator -> score()
        let result = lbl.handle_ordering(&mut best, None);

        // it will have successfully traversed the scoring logic.
        assert!(result.is_ok());
        Ok(())
    }
}
