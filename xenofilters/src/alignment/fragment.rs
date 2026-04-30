use crate::Penalty;
use anyhow::Result;
use rust_htslib::bam::record::Record;
use smallvec::SmallVec;
use crate::variant::{VntPerRec, VariantEval};
use crate::at::At;
use crate::alignment::AlignmentError;

pub(crate) struct Fragment<'r> {
    segment: SmallVec<[&'r Record; 2]>,
}

impl<'r> Fragment<'r> {
    pub(crate) fn new(segment: SmallVec<[&'r Record; 2]>) -> Self {
        Self { segment }
    }
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
        let mut variants: SmallVec<[&VariantEval; 4]> = delta
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

#[cfg(test)]
pub(crate) mod tests;
