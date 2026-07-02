//! CRAM input stream. namesorted only; no virtual-offset seeking.
//!
//! Requires: `--reference <fasta>` pointing to an fai-indexed FASTA.
//! NEEDS VERIFICATION: noodles 0.111.0 `cram::io::Reader` exact API
//! and `fasta::io::indexed_reader::Builder` availability.

use crate::{Error, aln_stream::AlignmentStream, config::Config, variant::StoreTrait};
use noodles::{cram, fasta, sam::{Header, alignment::record_buf::RecordBuf}};
use std::{path::Path, sync::Arc};

pub(crate) struct CramStream {
    header:  Header,
    records: Box<dyn Iterator<Item = Result<RecordBuf, Error>> + Send>,
    peeked:  Option<RecordBuf>,
    output:  crate::bam::OutputMode,
}

impl CramStream {
    pub(crate) fn new(path: &Path, reference: &Path) -> Result<Self, Error> {
        // NEEDS VERIFICATION: exact builder API for noodles 0.111.0 cram/fasta.
        let repo_reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(reference)?;
        let indexed_reader = fasta::repository::adapters::IndexedReader::new(repo_reader);
        let repo = fasta::Repository::new(indexed_reader);

        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repo)
            .build_from_path(path)?;
        let header = reader.read_file_header()?;

        // Collect records into a boxed iterator. Allocates per-record but
        // avoids lifetime entanglement between reader and struct.
        // For large files, replace with a channel-backed producer thread.
        let header2 = header.clone();
        let records: Vec<Result<RecordBuf, Error>> = reader
            .records(&header2)
            .map(|r| {
                r.map_err(Error::from).and_then(|cram_rec| {
                    // NEEDS VERIFICATION: cram::Record → RecordBuf conversion
                    RecordBuf::try_from_alignment_record(&header2, &cram_rec)
                        .map_err(Error::from)
                })
            })
            .collect();

        Ok(Self {
            header,
            records: Box::new(records.into_iter()),
            peeked: None,
            output: crate::bam::OutputMode::default(),
        })
    }
}

impl AlignmentStream<RecordBuf> for CramStream {
    fn next_qname(&self) -> &[u8] {
        self.peeked.as_ref().and_then(|r| r.name()).map_or(b"", |n| n.as_ref())
    }
    fn un_next(&mut self, rec: RecordBuf) -> Result<(), Error> {
        if self.peeked.is_some() { return Err(Error::CannotUnNext); }
        self.peeked = Some(rec); Ok(())
    }
    fn next_rec(&mut self) -> Result<Option<RecordBuf>, Error> {
        if let Some(r) = self.peeked.take() { return Ok(Some(r)); }
        self.records.next().transpose()
    }
    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<(), Error> {
        self.output.write(rec, is_best, &self.header)
    }
    fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<(), Error> { Ok(()) }
    fn variant_store(&self) -> Option<Arc<dyn StoreTrait>> { None }
    fn header(&self) -> &Header { &self.header }
}
