use super::core::LineByLine;
use crate::alignment::{Fragment, stringify_record, QualityAt};
use crate::alignment::FragmentState;
use anyhow::{Result, anyhow};
use noodles::sam::alignment::Record;
use smallvec::SmallVec;
use crate::alignment::MdCigFlags;
use crate::variant::Eval;

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
pub(crate) struct NeedlemanWunsch<'v> {
    pub(crate) prev: Vec<Cell>,
    pub(crate) curr: Vec<Cell>,
    pub(crate) dvnt_per_rec: SmallVec<[SmallVec<[Eval<'v>; 0]>; 8]>,
}

impl<'v> NeedlemanWunsch<'v> {
    pub(super) fn new(capacity: usize) -> Self {
        Self { prev: Vec::new(), curr: Vec::new(), dvnt_per_rec: SmallVec::with_capacity(capacity) }
    }

    pub(crate) fn resize(&mut self, new_len: usize) {
        self.prev.resize(new_len, Cell::default());
        self.curr.resize(new_len, Cell::default());
    }

    pub(crate) fn swap(&mut self) {
        std::mem::swap(&mut self.prev, &mut self.curr);
    }
}

impl<R: Record + PartialEq + QualityAt> LineByLine<R> {
    pub(super) fn score_candidate(
        &self,
        state: &FragmentState<R>,
        aln_idx: usize,
    ) -> Result<f64> {
        let mut nw = NeedlemanWunsch::new(state.get_records().len());
        let mut segment: SmallVec<[&R; 8]> = SmallVec::new();
        let mut md_cig_flags = SmallVec::with_capacity(state.get_records().len());
        let aln = self.aln.get(aln_idx).ok_or_else(|| anyhow!("No alignment for index {aln_idx}"))?;

        for idx in state.order_mates(&self.aln) {
            let rec = &state.get_records()[idx];
let flags = state.flags(idx).ok_or_else(|| anyhow!("No flags for record index {idx} in alignment {aln_idx}"))?;
            if flags.is_unmapped() {
                nw.dvnt_per_rec.push(SmallVec::new());
            } else {
                let header = aln.header();
                let tid = rec.reference_sequence_id(header).transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no reference sequence ID"))?;
                let start = rec.alignment_start().transpose()?
                    .ok_or_else(|| anyhow!("Mapped record has no alignment start"))?
                    .get();
                let end = start + rec.cigar().len();
                let delta_vars = aln.variant_store()
                    .map(|s| s.overlapping_multi(tid, start, end))
                    .unwrap_or_default();
                nw.dvnt_per_rec.push(delta_vars);
            }
            if !flags.is_secondary() {
                segment.push(rec);
                md_cig_flags.push(MdCigFlags::try_from_record(flags, rec)?);
            } else if flags.is_last_segment() {
                break;
            }
        }
        Fragment::new(&self.penalties, segment, md_cig_flags)?.score(nw).map_err(|e| {
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
