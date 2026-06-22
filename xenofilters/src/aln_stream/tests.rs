use crate::bam::{BamFormat, BamOutput, OutputMode};
use crate::config::{Config, StripReadSuffix};
use crate::tests::create_record;
use crate::variant::StoreTrait;
use crate::{AlignmentStream, AlnStream};
use anyhow::{anyhow, Result};
use noodles::sam::{alignment::record_buf::RecordBuf, Header};
use std::num::NonZeroUsize;
use std::sync::Arc;

pub(crate) struct MockStream {
    pub(crate) reads: Vec<RecordBuf>,
    pub(crate) original_reads: Vec<RecordBuf>,
    written: Vec<(RecordBuf, Option<bool>)>,
    aln_stream: AlnStream<RecordBuf>,
    i: usize,
}

impl AlignmentStream<RecordBuf> for MockStream {
    fn next_qname(&self) -> &[u8] {
        self.aln_stream.next_qname()
    }
    fn un_next(&mut self, rec: RecordBuf) -> Result<()> {
        self.un_next(rec)
    }
    fn next_rec(&mut self) -> Result<Option<RecordBuf>> {
        self.next_rec()
    }
    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<()> {
        self.write_record(rec, is_best)
    }
    fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<()> {
        Ok(())
    }
    fn variant_store(&self) -> Option<Arc<dyn StoreTrait>> {
        None
    }
    fn header(&self) -> &Header {
        &self.aln_stream.header
    }
}

impl MockStream {
    pub(crate) fn new(i: usize, reads: Vec<RecordBuf>) -> Self {
        let original_reads = reads.clone();
        let aln_stream = AlnStream {
            bam: None,
            next: None,
            sample_variants: None,
            population_variants: None,
            header: Header::default(),
            output_mode: OutputMode::MultiFile {
                output: None,
                filt: None,
                ambiguous: None,
            },
            threads: NonZeroUsize::MIN,
        };
        Self {
            reads,
            original_reads,
            written: Vec::new(),
            aln_stream,
            i,
        }
    }

    fn next_rec(&mut self) -> Result<Option<RecordBuf>> {
        if let Some(rec) = self.aln_stream.next_rec()? {
            return Ok(Some(rec));
        }
        if self.reads.is_empty() {
            return Ok(None);
        }
        let rec = self.reads.remove(0);
        self.aln_stream.un_next(rec)?;
        self.aln_stream.next_rec()
    }

    fn un_next(&mut self, rec: RecordBuf) -> Result<()> {
        let name = rec.name().expect("Invalid Name");
        eprintln!(
            "Un-next({}) read: {}",
            self.i,
            std::str::from_utf8(name.as_ref()).unwrap_or("Invalid UTF-8")
        );
        self.aln_stream.un_next(rec)
    }

    fn write_record(&mut self, rec: RecordBuf, state: Option<bool>) -> Result<()> {
        self.written.push((rec, state));
        Ok(())
    }
    fn fetch_by_virtual_offset(&mut self, virtual_offset: u64) -> Result<RecordBuf> {
        self.original_reads
            // Tests usually use small 0-based offsets matching the read index
            .get(virtual_offset as usize)
            .cloned()
            // Fallback to the first read just in case the mock offset is weird
            .or_else(|| self.original_reads.first().cloned())
            .ok_or_else(|| anyhow!("MockStream has no original reads to fetch"))
    }
}

fn empty_aln_stream() -> AlnStream<RecordBuf> {
    AlnStream {
        bam: None,
        next: None,
        sample_variants: None,
        population_variants: None,
        header: Header::default(),
        output_mode: OutputMode::MultiFile {
            output: None,
            filt: None,
            ambiguous: None,
        },
        threads: NonZeroUsize::MIN,
    }
}

#[test]
fn test_aln_stream_new_mismatch_strip_suffix_true_instead_of_false() {
    let mut config = Config {
        alignment: vec!["tests/data/test_input_1_a.bam".to_string()],
        stdout_format: BamFormat::Sam,
        strip_read_suffix: StripReadSuffix::True,
        ..Default::default()
    };
    assert!(AlnStream::<RecordBuf>::new(&mut config, 0).is_err());
}

#[test]
fn test_aln_stream_next_rec() -> Result<()> {
    let records = vec![
        create_record(b"read1/1", "10M", &[], &[30; 10], "10", false)?,
        create_record(b"read2/1", "10M", &[], &[30; 10], "10", false)?,
    ];
    let mut mock_stream = MockStream::new(0, records.clone());
    for expected in records {
        let rec = mock_stream.next_rec()?.unwrap();
        assert_eq!(rec.name(), expected.name());
    }
    assert!(mock_stream.next_rec()?.is_none());
    Ok(())
}

#[test]
fn test_aln_stream_un_next() -> Result<()> {
    let records = vec![
        create_record(b"read1/1", "10M", &[], &[30; 10], "10", false)?,
        create_record(b"read2/1", "10M", &[], &[30; 10], "10", false)?,
    ];
    let mut mock_stream = MockStream::new(0, records);

    let rec1 = mock_stream.next_rec()?.unwrap();
    assert_eq!(rec1.name().unwrap().as_ref() as &[u8], b"read1/1");
    mock_stream.un_next(rec1)?;

    let rec2 = mock_stream.next_rec()?.unwrap();
    assert_eq!(rec2.name().unwrap().as_ref() as &[u8], b"read1/1");
    let rec3 = mock_stream.next_rec()?.unwrap();
    assert_eq!(rec3.name().unwrap().as_ref() as &[u8], b"read2/1");
    Ok(())
}

#[test]
fn test_next_qname_empty_when_no_next_record() {
    assert_eq!(empty_aln_stream().next_qname(), b"");
}

#[test]
fn test_next_qname_returns_pending_records_name() -> Result<()> {
    let mut stream = empty_aln_stream();
    let rec = create_record(b"r1", "5M", &[], &[30; 5], "5", false)?;
    stream.un_next(rec)?;
    assert_eq!(stream.next_qname(), b"r1");
    Ok(())
}

#[test]
fn test_un_next_errors_when_already_occupied() -> Result<()> {
    let mut stream = empty_aln_stream();
    stream.un_next(create_record(b"r1", "5M", &[], &[30; 5], "5", false)?)?;
    assert!(stream
        .un_next(create_record(b"r2", "5M", &[], &[30; 5], "5", false)?)
        .is_err());
    Ok(())
}

#[test]
fn test_next_rec_none_when_empty_and_no_bam_reader() -> Result<()> {
    assert!(empty_aln_stream().next_rec()?.is_none());
    Ok(())
}

#[test]
fn test_write_record_is_noop_without_attached_writers() -> Result<()> {
    let mut stream = empty_aln_stream();
    let rec = create_record(b"r", "5M", &[], &[30; 5], "5", false)?;
    assert!(stream.write_record(rec.clone(), Some(true)).is_ok());
    assert!(stream.write_record(rec.clone(), Some(false)).is_ok());
    assert!(stream.write_record(rec, None).is_ok());
    Ok(())
}

#[test]
fn test_variant_store_none_when_unset() {
    assert!(empty_aln_stream().variant_store().is_none());
}
