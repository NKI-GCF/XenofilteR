use crate::{MAX_Q, Penalty};
use crate::variant::DeltaPerRec;
use anyhow::{Result, anyhow};
use crate::alignment::{UnifiedOpIterator, UnifiedOp};
use crate::variant::VariantEval;
use rust_htslib::bam::record::Record;
use smallvec::SmallVec;

#[derive(Clone, Copy)]
struct Cell {
    m: f64, // Match/Mismatch
    i: f64, // Insertion (gap in Alt)
    d: f64, // Deletion (gap in Read)
}

impl Default for Cell {
    fn default() -> Self {
        Self { m: -f64::INFINITY, i: -f64::INFINITY, d: -f64::INFINITY }
    }
}

pub(crate) struct At<'r, 's> {
    pen: &'r Penalty,
    seg: &'s SmallVec<[&'r Record; 2]>,
    seg_i: usize,
    refpos: i64,
    nt_i: usize,
    prev: Vec<Cell>,
    curr: Vec<Cell>,
}

impl<'r, 's> At<'r, 's> {
    pub(crate) fn new(pen: &'r Penalty, seg: &'s SmallVec<[&'r Record; 2]>) -> Self {
        Self {
            pen,
            seg,
            seg_i: 0,
            refpos: seg[0].pos(),
            nt_i: 0,
            prev: Vec::new(),
            curr: Vec::new(),
        }
    }
    pub(crate) fn update_seg_i(&mut self, seg_i: usize) {
        self.seg_i = seg_i;
        self.refpos = self.seg[seg_i].pos();
        self.nt_i = 0;
    }
    pub(crate) fn refpos(&self) -> i64 {
        self.refpos
    }
    pub(crate) fn nt_i(&self) -> usize {
        self.nt_i
    }
    fn segment_orientation_mismatch(&self, seg_i: usize) -> bool {
         seg_i != self.seg_i && self.seg[seg_i].is_reverse() != self.seg[self.seg_i].is_reverse()
    }
    fn q(&self, seg_i: usize, mut nt_i: usize) -> usize {
        let qual = self.seg[seg_i].qual();
        if self.segment_orientation_mismatch(seg_i) {
            nt_i = qual.len() - 1 - nt_i;
        }
        (qual[nt_i] as usize).min(MAX_Q - 1)
    }
    fn base(&self, seg_i: usize, nt_i: usize) -> u8 {
        if self.segment_orientation_mismatch(seg_i) {
            let encoded = self.seg[seg_i].seq().encoded;
            match encoded[encoded.len() - 1 - nt_i] {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            }
        } else {
            self.seg[seg_i].seq().encoded[nt_i]
        }
    }
    /// Score a variant's alt allele against read bases from a specific segment.
    /// `ref_start` and `ref_end` define the reference window to score within this segment.
    fn score_variant_in_seg(
        &mut self,
        vnt: &VariantEval,
        seg_i: usize,
        ref_start: i64,
        ref_end: i64,
    ) -> Option<f64> {
        let vnt_start = vnt.start();

        // Clamp the ref window to the variant's ref span
        let eff_ref_start = ref_start.max(vnt_start);
        let eff_ref_end = ref_end.min(vnt.ref_end());
        if eff_ref_start >= eff_ref_end { return None; }

        // Derive the alt slice from the ref offset (valid for MNPs; see caveat below)
        let vnt = vnt.vnt();
        let alt = vnt.alt_allele();
        let p_variant = vnt.p_variant();
        let ref_consumed = (eff_ref_start - vnt_start) as usize;
        let ref_slice_len = (eff_ref_end - eff_ref_start) as usize;
        let alt_slice = &alt[ref_consumed.min(alt.len())..(ref_consumed + ref_slice_len).min(alt.len())];

        if alt_slice.is_empty() && ref_slice_len == 0 { return None; }

        Some(self.score_variant_two_row_seg(alt_slice, ref_slice_len, seg_i, eff_ref_start, p_variant))
    }

    fn score_variant_two_row_seg(
        &mut self,
        alt: &[u8],
        ref_len: usize,
        seg_i: usize,
        ref_start: i64,
        p_variant: f64,  // passed in from vnt.vnt().p_variant()
    ) -> f64 {
        let n = alt.len();
        let m = n.max(ref_len);

        self.prev.resize(m + 1, Cell::default());
        self.curr.resize(m + 1, Cell::default());

        // Initialise first row: gaps in the read (deletions relative to alt)
        self.prev[0].m = 0.0;
        self.prev[0].i = f64::NEG_INFINITY;
        self.prev[0].d = f64::NEG_INFINITY;
        for j in 1..=m {
            self.prev[j].i = self.pen.gap_open + (j as f64 * self.pen.gap_extend);
            self.prev[j].m = f64::NEG_INFINITY;
            self.prev[j].d = f64::NEG_INFINITY;
        }

        let seg_ref_start = self.seg[seg_i].pos();
        let nt_i_base = (ref_start - seg_ref_start) as usize;

        for i in 1..=n {
            let alt_base = alt[i - 1];

            self.curr[0].d = self.pen.gap_open + (i as f64 * self.pen.gap_extend);
            self.curr[0].m = f64::NEG_INFINITY;
            self.curr[0].i = f64::NEG_INFINITY;

            for j in 1..=m {
                let nt_i = nt_i_base + j - 1;
                let read_base = self.base(seg_i, nt_i);
                let q = self.q(seg_i, nt_i);

                let lm = self.pen.log_likelihood_match[q];
                let lmm = self.pen.log_likelihood_mismatch[q];

                // Weight match/mismatch score by variant probability,
                // mirroring score_alt_match / score_ref_match
                let per_base_score = if alt_base == read_base {
                    p_variant * lm + (1.0 - p_variant) * lmm
                } else {
                    (1.0 - p_variant) * lm + p_variant * lmm
                };

                self.curr[j].m = per_base_score
                    + self.prev[j - 1].m
                        .max(self.prev[j - 1].i)
                        .max(self.prev[j - 1].d);

                self.curr[j].i = (self.curr[j - 1].m + self.pen.gap_open + self.pen.gap_extend)
                    .max(self.curr[j - 1].i + self.pen.gap_extend);

                self.curr[j].d = (self.prev[j].m + self.pen.gap_open + self.pen.gap_extend)
                    .max(self.prev[j].d + self.pen.gap_extend);
            }

            std::mem::swap(&mut self.prev, &mut self.curr);
        }

        self.prev[m].m.max(self.prev[m].i).max(self.prev[m].d)
    }

    fn score_variants_in_window(
        &mut self,
        vnt: &mut DeltaPerRec<'_>,
        seg_i: usize,
        start: i64,
        end: i64,
        ref_score: f64,
    ) {
        let mut i = 0;
        while i < vnt.len() && vnt[i].start() < end {
            if let Some(alt_score) = self.score_variant_in_seg(&vnt[i], seg_i, start, end) {
                vnt[i].update(ref_score, alt_score);
                let fully_processed = vnt[i].ref_end() <= end && vnt[i].alt_end() <= end;
                if !fully_processed {
                    i += 1;
                    continue;
                }
            }
            vnt.remove(i);
        }
    }

    pub(crate) fn score(&mut self, vnt: &mut DeltaPerRec<'_>) -> Result<f64> {
        let mut score = 0.0;

        // Score variant bases that fall within segments before the current one.
        for prior_seg_i in 0..self.seg_i {
            if self.seg[prior_seg_i].is_last_in_template() != self.seg[self.seg_i].is_last_in_template() {
                continue; // skip segments from different reads in the pair
            }
            let seg_ref_start = self.seg[prior_seg_i].pos();
            let seg_ref_end = seg_ref_start + self.seg[prior_seg_i].cigar().len() as i64;

            // prior (or following) segment bases contribute to alt scoring only, no incurrence.
            self.score_variants_in_window(vnt, prior_seg_i, seg_ref_start, seg_ref_end, 0.0);
        }

        let op_iter = UnifiedOpIterator::new(self.seg[self.seg_i])
            .map_err(|e| anyhow!("Failed to create UnifiedOpIterator: {}", e))?;

        for op_res in op_iter {
            let op = op_res?;
            let ref_start = self.refpos;
            let mut ref_score = 0.0;
            match op {
                UnifiedOp::Match(len) => {
                    for _ in 0..len {
                        ref_score += self.pen.log_likelihood_match[self.q(self.seg_i, self.nt_i)];
                        self.nt_i  += 1;
                        self.refpos += 1;
                    }
                },
                UnifiedOp::Mis(len) => {
                    for _ in 0..len {
                        ref_score += self.pen.log_likelihood_mismatch[self.q(self.seg_i, self.nt_i)];
                        self.nt_i  += 1;
                        self.refpos += 1;
                    }
                }
                UnifiedOp::Del(len) => {
                    self.refpos += len as i64;
                    ref_score += self.pen.gap_open + (len as f64) * self.pen.gap_extend;
                }
                UnifiedOp::Ins(len) => {
                    self.nt_i += len as usize;
                    ref_score += self.pen.gap_open + (len as f64) * self.pen.gap_extend;
                },
                UnifiedOp::Relocate { penalty_score, pos } => {
                    self.refpos = pos;
                    ref_score += penalty_score;
                }
                UnifiedOp::RefSkip(len) => {
                    self.refpos += len as i64;
                }
            }

            self.score_variants_in_window(vnt, self.seg_i, ref_start, self.refpos, ref_score);
            score += ref_score;
        }

        for next_seg_i in (self.seg_i + 1)..self.seg.len() {
            if self.seg[next_seg_i].is_last_in_template() != self.seg[self.seg_i].is_last_in_template() {
                break; // skip segments from different reads in the pair
            }
            let seg_ref_start = self.seg[next_seg_i].pos();
            let seg_ref_end = seg_ref_start + self.seg[next_seg_i].seq().len() as i64;

            self.score_variants_in_window(vnt, next_seg_i, seg_ref_start, seg_ref_end, 0.0);
        }

        Ok(score)
    }
}
