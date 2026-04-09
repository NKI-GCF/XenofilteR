#[cfg(test)]
use crate::alignment::stringify_record;
use crate::alignment::{AlignmentError, UnifiedOp, UnifiedOpIterator};
use crate::{MAX_Q, Penalty};
use anyhow::{Result, anyhow};
use rust_htslib::bam::record::{Cigar, Record};
use smallvec::SmallVec;
use crate::variant::Variant;

#[derive(Debug)]
struct Segment<'a> {
    rec: &'a Record,
    ref_start: i64,         // 0-based, inclusive
    ref_end: i64,           // exclusive
}

#[derive(Default)]
struct ReadEnd<'a> {
    segment: SmallVec<[Segment<'a>; 4]>,
    variant: SmallVec<[(f64, &'a dyn Variant); 8]>,
}

pub(crate) struct Fragment<'a> {
    read_end: SmallVec<[ReadEnd<'a>; 2]>,
    penalty: &'a Penalty,
}

impl<'a> Fragment<'a> {
    pub fn score(&mut self) -> Result<f64, AlignmentError> {
        let mut total = 0.0;
        for read_end in self.read_end.iter_mut() {
            total += read_end.score(self.penalty)?;
        }
        Ok(total)
    }
}

struct At {
    pub score: f64,
    pub pos: i64,
    pub offset: usize,
    pub seg_idx: usize,
}

impl At {
    fn increment_and_score(&mut self, score: f64) {
        self.score += score;
        self.pos += 1;
        self.offset += 1;
    }
}
impl<'a> ReadEnd<'a> {
    fn score(&mut self, penalty: &Penalty) -> Result<f64, AlignmentError> {
        let mut at = At { score: 0.0, pos: -1, offset: 0, seg_idx: 0 };

        // 1. Base Reference Scoring
        for seg in &self.segment {
            let op_iter = UnifiedOpIterator::new(seg.rec)
                .map_err(|e| anyhow!("Failed to create UnifiedOpIterator: {}", e))?;
            
            at.pos = seg.ref_start;
            at.offset = 0; 

            for op_res in op_iter {
                let op = op_res?;
                match op {
                    UnifiedOp::Mis(len) => {
                        for i in 0..len as usize {
                            let q = seg.rec.qual()[at.offset + i];
                            at.increment_and_score(penalty.log_likelihood_mismatch[(q as usize).min(MAX_Q - 1)]);
                        }
                    }
                    UnifiedOp::Match(len) => {
                        for i in 0..len as usize {
                            let q = seg.rec.qual()[at.offset + i];
                            at.increment_and_score(penalty.log_likelihood_match[(q as usize).min(MAX_Q - 1)]);
                        }
                    }
                    UnifiedOp::Del(len) => {
                        at.score += penalty.gap_open + (len as f64 - 1.0) * penalty.gap_extend;
                        at.pos += len as i64;
                    }
                    UnifiedOp::Ins(len) => {
                        at.score += penalty.gap_open + (len as f64 - 1.0) * penalty.gap_extend;
                        at.offset += len as usize;
                    }
                    UnifiedOp::RefSkip(len) => at.pos += len as i64,
                    UnifiedOp::Relocate { pos, penalty: penalty_score } => {
                        at.score -= penalty_score;
                        at.pos = pos;
                    }
                }
            }
            at.seg_idx += 1;
        }

        let base_reference_score = at.score;

        // 2. Evaluate Variants (Avoiding the Double Borrow)
        let mut deltas = Vec::with_capacity(self.variant.len());
        
        for &(_, variant) in &self.variant {
            // Pass ONLY the segment slice, not `self`
            let (read_bases, read_quals, base_penalty_incurred) = Self::extract_variant_context(&self.segment, variant, penalty)?;

            let mut delta = 0.0;
            if variant.matches_alt(&read_bases) {
                let variant_score = variant.score_alt_match(penalty, &read_quals);
                delta = variant_score - base_penalty_incurred;
            } else if variant.matches_ref(&read_bases) {
                let variant_score = variant.score_ref_match(penalty, &read_quals);
                delta = variant_score - base_penalty_incurred;
            }
            deltas.push(delta);
        }

        // Apply deltas back to self.variant
        for (v, d) in self.variant.iter_mut().zip(deltas) {
            v.0 = d;
        }

        // 3. Combinatorial Resolution
        let best_improvement = self.maximize_delta();
        
        Ok(base_reference_score + best_improvement)
    }

    /// Extracted into an associated function to decouple it from `&self`
    fn extract_variant_context(
        segments: &[Segment], 
        variant: &'a dyn Variant, 
        penalty: &Penalty
    ) -> Result<(Vec<u8>, Vec<u8>, f64), AlignmentError> {
        let mut bases = Vec::new();
        let mut quals = Vec::new();
        let mut incurred = 0.0;
        
        let v_start = variant.pos();
        let v_end = variant.end(); // exclusive
        
        for seg in segments {
            // Optimization: Skip segment entirely if it doesn't overlap variant bounds
            if seg.ref_end <= v_start || seg.ref_start >= v_end {
                continue;
            }
            
            let op_iter = UnifiedOpIterator::new(seg.rec)
                .map_err(|e| anyhow!("Failed to create UnifiedOpIterator: {}", e))?;
                
            let mut ref_pos = seg.ref_start;
            let mut offset = 0;
            
            for op_res in op_iter {
                let op = op_res?;
                
                match op {
                    UnifiedOp::Match(len) | UnifiedOp::Mis(len) => {
                        for i in 0..len as usize {
                            let curr_pos = ref_pos + i as i64;
                            if curr_pos >= v_start && curr_pos < v_end {
                                let q = seg.rec.qual()[offset + i];
                                let b = seg.rec.seq().encoded[offset + i]; 
                                bases.push(b);
                                quals.push(q);
                                
                                // Reverse exactly what was added in Phase 1
                                if let UnifiedOp::Mis(_) = op {
                                    incurred += penalty.log_likelihood_mismatch[(q as usize).min(MAX_Q - 1)];
                                } else {
                                    incurred += penalty.log_likelihood_match[(q as usize).min(MAX_Q - 1)];
                                }
                            }
                        }
                        ref_pos += len as i64;
                        offset += len as usize;
                    }
                    UnifiedOp::Ins(len) => {
                        // Insertions consume 0 ref bases, so they sit exactly at ref_pos
                        if ref_pos >= v_start && ref_pos < v_end {
                            for i in 0..len as usize {
                                let q = seg.rec.qual()[offset + i];
                                let b = seg.rec.seq().encoded[offset + i];
                                bases.push(b);
                                quals.push(q);
                            }
                            incurred += penalty.gap_open + (len as f64 - 1.0) * penalty.gap_extend;
                        }
                        offset += len as usize;
                    }
                    UnifiedOp::Del(len) => {
                        let end_pos = ref_pos + len as i64;
                        // If deletion overlaps the variant window
                        if ref_pos < v_end && end_pos > v_start {
                            incurred += penalty.gap_open + (len as f64 - 1.0) * penalty.gap_extend;
                        }
                        ref_pos += len as i64;
                    }
                    UnifiedOp::RefSkip(len) => ref_pos += len as i64,
                    UnifiedOp::Relocate { pos, .. } => ref_pos = pos,
                }
            }
        }
        
        Ok((bases, quals, incurred))
    }

    fn maximize_delta(&mut self) -> f64 {
        // Discard variants that worsen the score
        self.variant.retain(|v| v.0 > 0.0);
        if self.variant.is_empty() {
            return 0.0;
        }

        // Must sort by end coordinate for Weighted Interval Scheduling
        self.variant.sort_unstable_by_key(|v| v.1.end());

        let mut dp = vec![0.0; self.variant.len()];
        dp[0] = self.variant[0].0;

        for i in 1..self.variant.len() {
            let current_weight = self.variant[i].0;
            let current_start = self.variant[i].1.pos();
            
            let mut incl_weight = current_weight;
            
            // Rust's built-in binary search to find the latest non-overlapping interval
            // partition_point returns the index of the first element that evaluates to `false`.
            let latest_non_overlapping = self.variant[..i].partition_point(|v| v.1.end() <= current_start);
            
            if latest_non_overlapping > 0 {
                incl_weight += dp[latest_non_overlapping - 1];
            }

            dp[i] = dp[i - 1].max(incl_weight);
        }

        *dp.last().unwrap_or(&0.0)
    }
}

const fn revcmp_encoded(b: u8) -> u8 {
    match b {
        1 => 8,   // A -> T
        8 => 1,   // T -> A
        2 => 4,   // C -> G
        4 => 2,   // G -> C
        3 => 3,   // M -> M (A/C)
        5 => 5,   // R -> R (A/G)
        6 => 6,   // S -> S (C/G)
        7 => 7,   // V -> V (A/C/G)
        9 => 9,   // W -> W (A/T)
        10 => 10, // Y -> Y (C/T)
        11 => 11, // H -> H (A/C/T)
        12 => 12, // K -> K (G/T)
        13 => 13, // D -> D (A/G/T)
        14 => 14, // B -> B (C/G/T)
        15 => 15, // N -> N
        _ => b,   // = or garbage
    }
}

/// Calculates the penalty for a translocation.
/// This penalty is a trade-off between the cost of the unaligned/soft-clipped bases
/// and the quality of the aligned segment.
/// (Penalty for unaligned bases) - (Match Log-Likelihood Score)
fn calculate_translocation_penalty(penalty: &Penalty, record: &Record) -> Result<f64> {
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
                        total_match_log_likelihood += penalty.log_likelihood_match[q_idx];
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
                        clipped_quality_penalty += penalty.log_likelihood_mismatch[q_idx].abs();
                    }
                }
            }
            Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    let penalty = clipped_quality_penalty - total_match_log_likelihood;
    Ok(penalty.max(0.0))
}

pub fn build_fragment<'a>(
    penalty: &'a Penalty,
    records: &'a [Record],
    order: SmallVec<[usize; 2]>,
    variants_per_record: SmallVec<[SmallVec<[&'a dyn Variant; 8]>; 2]>,
) -> Result<Fragment<'a>> {
    let mut read_end = SmallVec::<[ReadEnd<'a>; 2]>::new();
    let mut current_read_end = ReadEnd { segment: SmallVec::new(), variant: SmallVec::new() };

    for &idx in &order {
        let record = &records[idx];
        if record.is_secondary() {
            if !record.is_first_in_template() { break; }
            continue;
        }

        let seg = Segment {
            rec: record,
            ref_start: record.pos() as i64,
            ref_end: record.cigar().end_pos() as i64,
        };

        current_read_end.segment.push(seg);
        if let Some(vars) = variants_per_record.get(idx) {
            current_read_end.variant.extend(vars.iter().map(|&v| (0.0, v)));
        }

        // When we hit a mate switch, finish the current read-end
        if !record.is_supplementary() && current_read_end.segment.len() > 1 {
            read_end.push(std::mem::take(&mut current_read_end));
        }
    }

    if !current_read_end.segment.is_empty() {
        read_end.push(current_read_end);
    }

    Ok(Fragment { read_end, penalty })
}

#[cfg(test)]
pub(crate) mod tests;
