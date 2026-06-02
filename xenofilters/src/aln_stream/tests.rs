use crate::tests::{BamFormat, create_record};
use crate::{AlignmentStream, AlnStream, Config, StripReadSuffix};
use anyhow::Result;
use noodles::sam::alignment::record_buf::RecordBuf as BamRecord;
use crate::variant::VariantStoreTrait;
use noodles::sam::header::Header;
use noodles::sam::alignment::Record as RecordTrait;

pub(crate) struct MockStream {
    pub(crate) reads: Vec<BamRecord>,
    written: Vec<(BamRecord, Option<bool>)>,
    aln_stream: AlnStream,
    i: usize,
}

impl AlignmentStream for MockStream {
    fn next_qname(&self) -> &[u8] {
        self.aln_stream.next_qname()
    }

    fn un_next(&mut self, rec: BamRecord) -> Result<()> {
        self.un_next(rec)
    }

    fn next_rec(&mut self) -> Result<Option<BamRecord>> {
        self.next_rec()
    }

    fn write_record(&mut self, rec: &dyn RecordTrait, is_best: Option<bool>) -> Result<()> {
        self.write_record(rec, is_best)
    }
    fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<()> {
        Ok(())
    }
    fn variant_store(&self) -> Option<&dyn VariantStoreTrait> {
        None
    }
    fn header(&self) -> &Header {
        &self.aln_stream.header
    }
}

impl MockStream {
    pub(crate) fn new(i: usize, reads: Vec<BamRecord>) -> Self {
        let aln_stream = AlnStream {
            ambiguous: None,
            bam: None,
            filt: None,
            next: None,
            output: None,
            sample_variants: None,
            population_variants: None,
            header: Header::default(),
        };
        Self {
            reads,
            written: Vec::new(),
            aln_stream,
            i
        }
    }
    fn next_rec(&mut self) -> Result<Option<BamRecord>> {
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
    fn un_next(&mut self, rec: BamRecord) -> Result<()> {
        let name = rec.name().expect("Invalid Name");
        eprintln!(
            "Un-next({}) read: {}",
            self.i,
            std::str::from_utf8(name.as_ref()).unwrap_or("Invalid UTF-8")
        );
        self.aln_stream.un_next(rec)
    }
    fn write_record(&mut self, rec: &dyn RecordTrait, state: Option<bool>) -> Result<()> {
        self.written.push((rec, state));
        Ok(())
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

    let result = AlnStream::new(&mut config, 0);
    assert!(result.is_err());
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
    let name1 = rec1.name().unwrap();
    assert_eq!(name1.as_ref(), b"read1/1");

    mock_stream.un_next(rec1)?;

    let rec2 = mock_stream.next_rec()?.unwrap();
    let name2 = rec2.name().unwrap();
    assert_eq!(name2.as_ref(), b"read1/1");

    let rec3 = mock_stream.next_rec()?.unwrap();
    let name3 = rec3.name().unwrap();
    assert_eq!(name3.as_ref(), b"read2/1");

    Ok(())
}
