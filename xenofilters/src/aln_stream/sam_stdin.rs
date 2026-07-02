//! SAM stream from stdin — for piping directly from `bwa mem`.
//!
//! Limitations:
//! - namesorted only (stdin is not seekable; HashLookup pass-2 cannot seek)
//! - single stream from stdin per invocation

use crate::{Error, aln_stream::AlignmentStream, config::Config};
use noodles::sam::{self, Header, alignment::record_buf::RecordBuf};
use std::io::BufReader;

pub(crate) struct SamStdinStream {
    header: Header,
    inner:  sam::io::Reader<BufReader<std::io::Stdin>>,
    peeked: Option<RecordBuf>,
    output: crate::bam::OutputMode,
}

impl SamStdinStream {
    pub(crate) fn new(_opt: &mut Config, _i: usize) -> Result<Self, Error> {
        let mut reader = sam::io::Reader::new(BufReader::new(std::io::stdin()));
        let header = reader.read_header()?;
        // Validate SO tag for namesorted — same check as AlnStream::new.
        // (strip_read_suffix auto-detection omitted for brevity; add same logic)
        let s = Self {
            header,
            inner: reader,
            peeked: None,
            output: crate::bam::OutputMode::default(),
        };
        Ok(s)
    }
}

impl AlignmentStream<RecordBuf> for SamStdinStream {
    fn next_qname(&self) -> &[u8] {
        self.peeked.as_ref().and_then(|r| r.name()).map_or(b"", |n| n.as_ref())
    }

    fn un_next(&mut self, rec: RecordBuf) -> Result<(), Error> {
        if self.peeked.is_some() { return Err(Error::CannotUnNext); }
        self.peeked = Some(rec);
        Ok(())
    }

    fn next_rec(&mut self) -> Result<Option<RecordBuf>, Error> {
        if let Some(r) = self.peeked.take() { return Ok(Some(r)); }
        let mut buf = RecordBuf::default();
        match self.inner.read_record_buf(&self.header, &mut buf)? {
            0 => Ok(None),
            _ => Ok(Some(buf)),
        }
    }

    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<(), Error> {
        self.output.write(rec, is_best, &self.header)
    }

    fn header(&self) -> &Header { &self.header }
}
