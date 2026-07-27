use super::FragmentTable;
use super::stream::{RecordKind, StreamAccumulator, StreamKind};
use crate::Error;
use crate::region::{ambiguous::AmbiguousRegions, diagnostic::SegregateVariants};
use smallvec::SmallVec;

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
        bed: Option<&AmbiguousRegions>,
        vcf: Option<&SegregateVariants>,
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
        bed: Option<&AmbiguousRegions>,
        vcf: Option<&SegregateVariants>,
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
    key: Box<[u8]>,
    nr: usize,
    seq_counter: &mut u64,
    bed: Option<&AmbiguousRegions>,
    vcf: Option<&SegregateVariants>,
) -> Result<(Box<[u8]>, bool), Error> {
    let entry = table.entry(key.clone()).or_insert_with(|| {
        let sn = *seq_counter;
        *seq_counter += 1;
        PendingFragment::new(sn)
    });
    let complete = entry.push(rec, nr, bed, vcf)?;
    Ok((key, complete))
}
