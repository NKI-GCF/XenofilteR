//! `SimpleRec` — the minimal interface `Fragment` requires from a BAM/SAM record.
//!
//! Implemented for both `noodles::bam::Record` (zero-copy, lazy-decoded) and
//! `noodles::sam::alignment::RecordBuf` (eager, owned).
//! Static dispatch everywhere — no `dyn SimpleRec` in the hot path.

use noodles::sam::{Header, alignment::record_buf::RecordBuf};

pub trait SimpleRec:
    noodles::sam::alignment::Record + PartialEq
{
    fn quality_at(&self, i: usize) -> Option<u8>;
    fn ref_seq_id(&self) -> Option<Result<usize, std::io::Error>>;
    fn as_record_buf(&self, header: &Header) -> Result<RecordBuf, std::io::Error>;
}

impl SimpleRec for noodles::bam::Record {
    fn quality_at(&self, i: usize) -> Option<u8> {
        self.quality_scores().as_ref().get(i).copied()
    }
    fn ref_seq_id(&self) -> Option<Result<usize, std::io::Error>> {
        self.reference_sequence_id()
    }
    fn as_record_buf(&self, header: &Header) -> Result<RecordBuf, std::io::Error> {
        RecordBuf::try_from_alignment_record(header, self)
    }
}

impl SimpleRec for RecordBuf {
    fn quality_at(&self, i: usize) -> Option<u8> {
        self.quality_scores().as_ref().get(i).copied()
    }
    fn ref_seq_id(&self) -> Option<Result<usize, std::io::Error>> {
        self.reference_sequence_id().map(Ok)
    }
    fn as_record_buf(&self, _header: &Header) -> Result<RecordBuf, std::io::Error> {
        Ok(self.clone())
    }
}
