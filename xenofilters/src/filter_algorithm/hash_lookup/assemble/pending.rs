use crate::region::tabix_query::{TabixBed, TabixVcf};
use smallvec::SmallVec;
use super::FragmentTable;
use super::stream::{StreamAccumulator, StreamKind, RecordKind};
use crate::Error;

// ---------------------------------------------------------------------------
// PendingFragment -- unchanged structurally; type updated to RecordKind
// ---------------------------------------------------------------------------

pub struct PendingFragment {
    driving_buf: StreamAccumulator,
    lookup_buf: StreamAccumulator,
    pub(crate) driving: StreamKind,
    pub(crate) lookup: StreamKind,
    pub(crate) supplementary_offsets: [SmallVec<[u64; 1]>; 2],
    pub(crate) seq_nr: u64,
    pub(crate) is_paired: Option<bool>,
}

impl PendingFragment {
    pub(crate) fn new(seq_nr: u64) -> Self {
        Self {
            driving_buf: StreamAccumulator::default(),
            lookup_buf: StreamAccumulator::default(),
            driving: StreamKind::Empty,
            lookup: StreamKind::Empty,
            supplementary_offsets: [SmallVec::new(), SmallVec::new()],
            seq_nr,
            is_paired: None,
        }
    }

    pub(crate) fn push(
        &mut self,
        rec: RecordKind,
        nr: usize,
        bed: Option<&TabixBed>,
        vcf: Option<&TabixVcf>,
    ) -> Result<bool, Error> {
        if self.is_paired.is_none() {
            self.is_paired = Some(rec.flags().is_segmented());
        }
        if rec.is_supplementary() {
            self.supplementary_offsets[nr].push(rec.virtual_offset());
            return self.check_complete(bed, vcf);
        }
        if nr == 0 {
            self.driving_buf.push(rec);
        } else {
            self.lookup_buf.push(rec);
        }
        self.check_complete(bed, vcf)
    }

    fn expected_primaries(&self) -> usize {
        self.is_paired.map_or(1, |p| p as usize + 1)
    }

    fn check_complete(
        &mut self,
        bed: Option<&TabixBed>,
        vcf: Option<&TabixVcf>,
    ) -> Result<bool, Error> {
        let exp = self.expected_primaries();
        if self.driving.is_empty() && self.driving_buf.primary_count >= exp {
            let buf = std::mem::take(&mut self.driving_buf);
            self.driving = buf.classify(bed, vcf)?;
        }
        if self.lookup.is_empty() && self.lookup_buf.primary_count >= exp {
            let buf = std::mem::take(&mut self.lookup_buf);
            self.lookup = buf.classify(bed, vcf)?;
        }
        Ok(!self.driving.is_empty() && !self.lookup.is_empty())
    }
}

pub(crate) fn insert(
    table: &mut FragmentTable,
    rec: RecordKind,
    name: String,
    nr: usize,
    seq_counter: &mut u64,
    bed: Option<&TabixBed>,
    vcf: Option<&TabixVcf>,
) -> Result<(String, bool), Error> {
    let entry = table.entry(name.clone()).or_insert_with(|| {
        let sn = *seq_counter;
        *seq_counter += 1;
        PendingFragment::new(sn)
    });
    let complete = entry.push(rec, nr, bed, vcf)?;
    Ok((name, complete))
}
