use crate::alignment::AlignmentError;
use crate::alignment::{UnifiedOp, UnifiedOpIterator};
use crate::{MAX_Q, Penalties};
use anyhow::{Result, anyhow};
use rust_htslib::bam::record::{Cigar, Record};

#[allow(dead_code)]
// A new struct to hold the combined alignment data
pub struct StitchedFragment<'a> {
    pub seq: Box<dyn Iterator<Item = u8> + 'a>,
    pub qual: Box<dyn Iterator<Item = u8> + 'a>,
    pub tid: i32,
    pub pos: i64,
    pub ops: Box<dyn Iterator<Item = Result<UnifiedOp, AlignmentError>> + 'a>,
    penalties: &'a Penalties,
}

impl<'a> StitchedFragment<'a> {
    pub fn score(&mut self) -> Result<f64, AlignmentError> {
        let mut total_score = 0.0;
        for op in self.ops.by_ref() {
            //#[cfg(test)]
            //eprintln!("Operation: {:?}", op);
            match op? {
                UnifiedOp::Match(len) => {
                    let len_usize = len as usize;
                    for _ in 0..len_usize {
                        if let Some(q) = self.qual.next() {
                            let q_idx = q as usize;
                            if q_idx < MAX_Q {
                                total_score += self.penalties.log_likelihood_match[q_idx];
                            }
                        }
                    }
                }
                UnifiedOp::Mis(len) => {
                    let len_usize = len as usize;
                    for _ in 0..len_usize {
                        if let Some(q) = self.qual.next() {
                            let q_idx = q as usize;
                            if q_idx < MAX_Q {
                                total_score += self.penalties.log_likelihood_mismatch[q_idx];
                            }
                        }
                    }
                }
                UnifiedOp::Ins(len) => {
                    let len_usize = len as usize;
                    // Insertions do not consume reference bases, only read bases
                    for _ in 0..len_usize {
                        self.qual.next();
                    }
                    total_score +=
                        self.penalties.gap_open + (len as f64 - 1.0) * self.penalties.gap_extend;
                }
                UnifiedOp::Del(len) => {
                    // Deletions consume reference bases, not read bases
                    total_score +=
                        self.penalties.gap_open + (len as f64 - 1.0) * self.penalties.gap_extend;
                }
                UnifiedOp::RefSkip(_len) => {
                    // RefSkips do not affect score
                }
                UnifiedOp::Trans(penalty) => {
                    total_score -= penalty;
                }
            }
        }
        Ok(total_score)
    }
}

const fn revcmp_encoded(b: u8) -> u8 {
    match b {
        1 => 8,  // A -> T
        8 => 1,  // T -> A
        2 => 4,  // C -> G
        4 => 2,  // G -> C
        3 => 3,  // M -> M (A/C)
        5 => 5,  // R -> R (A/G)
        6 => 6,  // S -> S (C/G)
        7 => 7,  // V -> V (A/C/G)
        9 => 9,  // W -> W (A/T)
        10 => 10,// Y -> Y (C/T)
        11 => 11,// H -> H (A/C/T)
        12 => 12,// K -> K (G/T)
        13 => 13,// D -> D (A/G/T)
        14 => 14,// B -> B (C/G/T)
        15 => 15,// N -> N
        _ => b,  // = or garbage
    }
}

fn get_read_iterators<'a>(
    rec: &'a Record,
) -> (
    Box<dyn Iterator<Item = u8> + 'a>,
    Box<dyn Iterator<Item = u8> + 'a>,
) {
    let seq_encoded = rec.seq().encoded;
    let qual = rec.qual();
    //#[cfg(test)]
    //eprintln!("Seq encoded: {:?}, {:?}", seq_encoded, qual);

    if rec.is_reverse() {
        (
            Box::new(seq_encoded.iter().rev().map(|&b| revcmp_encoded(b))),
            Box::new(qual.iter().rev().copied()),
        )
    } else {
        (
            Box::new(seq_encoded.iter().copied()),
            Box::new(qual.iter().copied()),
        )
    }
}

/// Calculates the penalty for a translocation.
/// This penalty is a trade-off between the cost of the unaligned/soft-clipped bases
/// and the quality of the aligned segment.
/// (Penalty for unaligned bases) - (Match Log-Likelihood Score)
fn calculate_translocation_penalty(penalties: &Penalties, record: &Record) -> Result<f64> {
    let cigar_view = record.cigar();

    let mut qual_iter = record.qual().iter().copied();

    let mut clipped_quality_penalty = 0.0;
    let mut total_match_log_likelihood = 0.0;

    for op in cigar_view.iter() {
        let len = op.len() as usize;
        match op {
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                for q in qual_iter.by_ref().take(len) {
                    let q_idx = q as usize;
                    if q_idx < MAX_Q {
                        total_match_log_likelihood += penalties.log_likelihood_match[q_idx];
                    }
                }
            }
            Cigar::Ins(_) => {
                qual_iter.by_ref().nth(len.saturating_sub(1));
            }
            Cigar::SoftClip(_) => {
                for q in qual_iter.by_ref().take(len) {
                    let q_idx = q as usize;
                    if q_idx < MAX_Q {
                        clipped_quality_penalty += penalties.log_likelihood_mismatch[q_idx].abs();
                    }
                }
            }
            Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    let penalty = clipped_quality_penalty - total_match_log_likelihood;
    Ok(penalty.max(0.0))
}

pub fn stitched_fragment<'a>(
    penalties: &'a Penalties,
    records: &'a [Record],
    order: Vec<usize>,
) -> Result<StitchedFragment<'a>> {
    // Hard clipped may occur in non-primary, so because we need all seq and qual:
    let mut primary_it = order
        .iter()
        .map(|&i| &records[i])
        .filter(|r| !r.is_supplementary() && !r.is_secondary());

    let mut stitched = if let Some(anchor) = primary_it.next() {
        let (mut seq, mut qual) = get_read_iterators(anchor);
        if let Some(mate) = primary_it.next() {
            let (s, q) = get_read_iterators(mate);
            seq = Box::new(seq.chain(s));
            qual = Box::new(qual.chain(q));
        }
        StitchedFragment {
            seq,
            qual,
            tid: anchor.tid(),
            pos: anchor.pos(),
            ops: Box::new(std::iter::empty()),
            penalties,
        }
    } else {
        return Err(anyhow!("No primary alignment found"));
    };

    let mut first_record = true;

    for record in order.iter().map(|&i| &records[i]) {
        if record.is_secondary() {
            if !record.is_first_in_template() {
                break;
            }
            continue;
        }

        if first_record {
            first_record = false;
        } else if record.is_supplementary() {
            let penalty = calculate_translocation_penalty(stitched.penalties, record)?;
            stitched.ops = Box::new(
                stitched
                    .ops
                    .chain(std::iter::once(Ok(UnifiedOp::Trans(penalty)))),
            )
        }

        let ops = UnifiedOpIterator::new(record)?;

        stitched.ops = Box::new(stitched.ops.chain(ops));
    }
    Ok(stitched)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Config;
    use crate::tests::create_record;
    use rust_htslib::bam::record::CigarString;

    fn setup_penalties() -> Penalties {
        let c = Config::default();
        let mut p = c.to_penalties();
        p.log_likelihood_match = [0.0; 93]; // 0 log-lik for match
        p.log_likelihood_mismatch = [-1.0; 93];
        p.gap_open = -2.0;
        p.gap_extend = -0.5;
        p
    }

    #[test]
    fn test_stitched_fragment_scoring() {
        let penalties = setup_penalties();
        let mut fragment = StitchedFragment {
            seq: Box::new(std::iter::empty()),
            qual: Box::new(vec![30, 30, 30, 30].into_iter()),
            tid: 0,
            pos: 100,
            ops: Box::new(
                vec![
                    Ok(UnifiedOp::Match(2)),
                    Ok(UnifiedOp::Mis(1)),
                    Ok(UnifiedOp::Ins(1)), // Gap open
                    Ok(UnifiedOp::Del(2)), // Gap open + 1 extend
                ]
                .into_iter(),
            ),
            penalties: &penalties,
        };

        let score = fragment.score().unwrap();
        // Match(2)*0 + Mis(1)*-1 + Ins(1)*-2 + Del(2)*(-2 + -0.5) = -5.5
        assert_eq!(score, -5.5);
    }

    #[test]
    fn test_translocation_penalty_logic() {
        let penalties = setup_penalties();
        // Create a record with 10M (perfect) and 5S (soft clip)
        // High quality on soft clip = high penalty
        // High quality on matches = offsets penalty (lowers it)
        let mut record = Record::new();
        let vec_cig = vec![Cigar::Match(10), Cigar::SoftClip(5)];
        let cigar = CigarString(vec_cig);
        record.set(b"read1", Some(&cigar), &[b'A'; 15], &[30; 15]);

        let penalty = calculate_translocation_penalty(&penalties, &record).unwrap();

        // 5 * |-1.0| (mismatch for softclips) - 10 * 0.0 (match for matches) = 5.0
        assert!(penalty > 0.0);
        assert_eq!(penalty, 5.0);
    }

    #[test]
    fn test_stitched_fragment_creation() -> Result<()> {
        let penalties = setup_penalties();

        let record1 = create_record(b"read1", "5M3S", &[b'A'; 8], &[30; 8], "5", false)?;
        let record2 = create_record(b"read1", "4M4S", &[b'A'; 8], &[30; 8], "4", false)?;

        let records = vec![record1, record2];
        let order = vec![0, 1];

        let stitched = stitched_fragment(&penalties, &records, order).unwrap();

        // Check that the stitched fragment has the correct TID and POS
        assert_eq!(stitched.tid, records[0].tid());
        assert_eq!(stitched.pos, records[0].pos());
        Ok(())
    }
    #[test]
    fn test_get_read_iterators() -> Result<()> {
        let record = create_record(
            b"read1",
            "4M",
            &[b'A', b'C', b'G', b'T'],
            &[30, 31, 32, 33],
            "4",
            false,
        )?;

        let (seq_iter, qual_iter) = get_read_iterators(&record);

        let seq: Vec<u8> = seq_iter.collect();
        let qual: Vec<u8> = qual_iter.collect();

        assert_eq!(seq, record.seq().encoded.to_vec());
        assert_eq!(qual, vec![30, 31, 32, 33]);
        Ok(())
    }

    #[test]
    fn test_get_read_iterators_reverse() -> Result<()> {
        let record = create_record(
            b"read1",
            "4M",
            &[b'A', b'C', b'G', b'T'],
            &[30, 31, 32, 33],
            "4",
            true,
        )?;

        let (seq_iter, qual_iter) = get_read_iterators(&record);

        let seq: Vec<u8> = seq_iter.collect();
        let qual: Vec<u8> = qual_iter.collect();

        let expected: Vec<u8> = record
            .seq()
            .encoded
            .iter()
            .rev()
            .map(|&b| revcmp_encoded(b))
            .collect();

        assert_eq!(
            seq,
            expected,
        );
        assert_eq!(qual, vec![33, 32, 31, 30]);
        Ok(())
    }

    #[test]
    fn test_revcmp_encoded() {
        // A <-> T
        assert_eq!(revcmp_encoded(1), 8);
        assert_eq!(revcmp_encoded(8), 1);

        // C <-> G
        assert_eq!(revcmp_encoded(2), 4);
        assert_eq!(revcmp_encoded(4), 2);

        // N stays N
        assert_eq!(revcmp_encoded(15), 15);

        // '=' or garbage stays unchanged
        assert_eq!(revcmp_encoded(0), 0);
        assert_eq!(revcmp_encoded(42), 42);
    }

    #[test]
    fn test_calculate_translocation_penalty() -> Result<()> {
        let penalties = setup_penalties();
        let record = create_record(b"read1", "6M4S", &[b'A'; 10], &[30; 10], "6", false)?;

        let penalty = calculate_translocation_penalty(&penalties, &record).unwrap();
        // 4 * |-1.0| - 6 * 0.0 = 4.0
        assert_eq!(penalty, 4.0);
        Ok(())
    }

    #[test]
    fn test_calculate_translocation_penalty_no_softclip() -> Result<()> {
        let penalties = setup_penalties();
        let record = create_record(b"read1", "10M", &[b'A'; 10], &[30; 10], "10", false)?;

        let penalty = calculate_translocation_penalty(&penalties, &record).unwrap();
        // No soft clips, so penalty should be 0.0
        assert_eq!(penalty, 0.0);
        Ok(())
    }

    #[test]
    fn test_calculate_translocation_penalty_high_quality_softclip() -> Result<()> {
        let penalties = setup_penalties();
        let record = create_record(
            b"read1",
            "5M5S",
            &[b'A'; 10],
            &[40; 10], // High quality scores
            "5",
            false,
        )?;

        let penalty = calculate_translocation_penalty(&penalties, &record).unwrap();
        // 5 * |-1.0| - 5 * 0.0 = 5.0
        assert_eq!(penalty, 5.0);
        Ok(())
    }
}
