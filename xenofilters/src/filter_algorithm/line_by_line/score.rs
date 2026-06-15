use super::core::LineByLine;
use crate::alignment::{Fragment, stringify_record, SimpleRec};
use crate::alignment::FragmentState;
use anyhow::{Result, anyhow};
use smallvec::SmallVec;
use crate::alignment::MdCigFlags;

#[derive(Clone, Copy, PartialOrd, PartialEq)]
pub(crate) struct Cell {
    pub(crate) m: f64, // Match/Mismatch
    pub(crate) i: f64, // Insertion (gap in Alt)
    pub(crate) d: f64, // Deletion (gap in Read)
}

impl Cell {
    pub(crate) fn reinit(&mut self, gap_open: f64, gap_extend: f64, i: i32) {
        self.m = f64::NEG_INFINITY;
        self.i = f64::NEG_INFINITY;
        self.d = f64::NEG_INFINITY;
        match i {
            0 => self.m = 0.0,
            i if i < 0 => self.i = gap_open + (i.abs() as f64) * gap_extend,
            i => self.d = gap_open + (i as f64) * gap_extend,
        }
    }
}

impl Default for Cell {
    fn default() -> Self {
        Self { m: -f64::INFINITY, i: -f64::INFINITY, d: -f64::INFINITY }
    }
}

#[derive(Clone)]
pub(crate) struct NeedlemanWunsch {
    pub(crate) prev: Vec<Cell>,
    pub(crate) curr: Vec<Cell>,
}

impl NeedlemanWunsch {
    pub(super) fn new() -> Self {
        Self { prev: Vec::new(), curr: Vec::new() }
    }

    pub(crate) fn resize(&mut self, new_len: usize) {
        self.prev.clear();
        self.curr.clear();
        self.prev.resize(new_len, Cell::default());
        self.curr.resize(new_len, Cell::default());
    }

    pub(crate) fn swap(&mut self) {
        std::mem::swap(&mut self.prev, &mut self.curr);
    }
}

impl<R: SimpleRec> LineByLine<R> {
    pub(super) fn score_candidate(
        &mut self,
        state: &FragmentState<R>,
        aln_idx: usize,
    ) -> Result<f64> {
        self.nw_scratch.resize(state.get_records().len());
        let mut segment: SmallVec<[&R; 8]> = SmallVec::new();
        let mut md_cig_flags = SmallVec::with_capacity(state.get_records().len());
        let aln = self.aln.get(aln_idx).ok_or_else(|| anyhow!("No alignment for index {aln_idx}"))?;
        let mut dvnt_per_rec = SmallVec::with_capacity(state.get_records().len());

        for idx in state.order_mates() {
            let rec = &state.get_records()[idx];
let flags = state.flags(idx).ok_or_else(|| anyhow!("No flags for record index {idx} in alignment {aln_idx}"))?;
            if flags.is_unmapped() {
                dvnt_per_rec.push(SmallVec::new());
            } else {
                let tid = rec.ref_seq_id().transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no reference sequence ID"))?;
                let start = rec.alignment_start().transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no alignment start"))?
                    .get();
                let end = start + rec.cigar().len();
                let delta_vars = aln.variant_store()
                    .map(|s| s.overlapping_multi(tid, start, end))
                    .unwrap_or_default();
                dvnt_per_rec.push(delta_vars);
            }
            if !flags.is_secondary() {
                segment.push(rec);
                md_cig_flags.push(MdCigFlags::try_from_record(rec, flags)?);
            } else if flags.is_last_segment() {
                break;
            }
        }
        Fragment::new(&self.penalties, segment, md_cig_flags)?.score(&mut self.nw_scratch, dvnt_per_rec).map_err(|e| {
            anyhow!(
                "Error scoring fragment for alignment {aln_idx}: {}\n{}",
                e,
                state
                    .get_records()
                    .iter()
                    .map(stringify_record)
                    .collect::<Vec<String>>()
                    .join("\n")
            )
        })
    }
}
