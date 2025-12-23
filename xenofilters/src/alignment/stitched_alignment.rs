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

const fn revcmp(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'N' => b'N',
        _ => base, // Non-standard bases are returned as-is
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

    if rec.is_reverse() {
        (
            Box::new(seq_encoded.iter().rev().map(|&b| revcmp(b))),
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
            ops: Box::new(vec![
                Ok(UnifiedOp::Match(2)),
                Ok(UnifiedOp::Mis(1)),
                Ok(UnifiedOp::Ins(1)), // Gap open
                Ok(UnifiedOp::Del(2)), // Gap open + 1 extend
            ].into_iter()),
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
}
