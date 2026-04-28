use crate::Penalty;
use anyhow::Result;
use rust_htslib::bam::record::Record;
use smallvec::SmallVec;
use crate::variant::{VntPerRec, VariantEval};
use crate::at::At;
use crate::alignment::{AlignmentError, SegmentVec};

pub(crate) struct Fragment<'r> {
    segment: SegmentVec<'r>,
}

impl<'r> Fragment<'r> {
    pub(crate) fn score<'v, 's>(&'s self, pen: &'r Penalty, mut vnt_per_rec: VntPerRec<'v>) -> Result<f64, AlignmentError> 
        where 'r: 's, 'v: 's
    {
        let mut at = At::new(pen, &self.segment);
        let mut score = at.score(&mut vnt_per_rec[0])?;

        for i in 1..self.segment.len() {
            at.update_seg_i(i);
            score += at.score(&mut vnt_per_rec[i])?;
        }
        Ok(score + self.maximize_delta(vnt_per_rec))
    }

    fn maximize_delta<'v>(&self, delta: VntPerRec<'v>) -> f64 {
        let mut variants: SmallVec<[&VariantEval<'v>; 4]> = delta
            .iter()
            .flatten()
            .filter(|v| v.delta() > 0.0)
            .collect();

        if variants.is_empty() {
            return 0.0;
        }

        // Sort by end coordinate for Weighted Interval Scheduling
        variants.sort_unstable_by_key(|v| v.end());

        let mut dp = vec![0.0; variants.len()];
        dp[0] = variants[0].delta();

        for i in 1..variants.len() {
            let current_weight = variants[i].delta();
            let current_start = variants[i].start();

            let mut incl_weight = current_weight;

            let latest = variants[..i].partition_point(|v| v.end() <= current_start);
            if latest > 0 {
                incl_weight += dp[latest - 1];
            }

            dp[i] = dp[i - 1].max(incl_weight);
        }

        *dp.last().unwrap_or(&0.0)
    }
}

pub(crate) struct Segment<'a> {
    pub(crate) rec: &'a Record,
}

impl <'a> Segment<'a> {
    pub(crate) fn ref_start(&self) -> i64 {
        self.rec.pos()
    }
    pub(crate) fn ref_end(&self) -> i64 {
        self.rec.cigar().end_pos()
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

pub fn build_fragment<'r>(
    records: &'r [Record],
    order: SmallVec<[usize; 2]>,
) -> Result<Fragment<'r>> {
    let mut segment = SegmentVec::<'r>::new();

    for &idx in &order {
        let rec = &records[idx];
        if !rec.is_secondary() {
            segment.push(Segment { rec });
        } else if !rec.is_first_in_template() {
            break;
        }
    }
    Ok(Fragment { segment })
}

#[cfg(test)]
pub(crate) mod tests;
