// src/aln_stream/tests.rs
use crate::aln_stream::OutputMode;
use crate::bam::BamFormat;
use crate::config::{Config, StripReadSuffix};
use crate::tests::create_record;
use crate::variant::StoreTrait;
use crate::{AlignmentStream, AlnStream};
use anyhow::Result;
use noodles::sam::header::record::value::{map::ReadGroup, Map};
use noodles::sam::{alignment::record_buf::RecordBuf, header::Header};
use std::sync::Arc;

pub(crate) struct MockStream {
    pub(crate) reads: Vec<RecordBuf>,
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
        Self {
            reads,
            written: Vec::new(),
            aln_stream: empty_aln_stream(),
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
        self.aln_stream.un_next(rec.into())?;
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
}

fn empty_aln_stream() -> AlnStream<RecordBuf> {
    AlnStream {
        bam: None,
        output_mode: OutputMode::MultiFile {
            ambiguous: None,
            filt: None,
            output: None,
        },
        next: None,
        variant_store: None,
        header: Header::default(),
        threads: 1,
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
    stream.un_next(create_record(b"r1", "5M", &[], &[30; 5], "5", false)?)?;
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

#[test]
fn test_init_writers_configures_merged_mode() -> Result<()> {
    let temp_dir = tempfile::tempdir()?;
    let merged_path = temp_dir.path().join("merged.bam");

    let mut config = Config::default();
    config.merged_output = Some(merged_path.clone());
    config.no_program_line = true; // simplify header checks

    let mut stream = empty_aln_stream();

    // Inject a dummy read group to test that header expansion is triggered
    // during init_writers.
    stream
        .header
        .read_groups_mut()
        .insert("rg_test".parse()?, Map::<ReadGroup>::default());

    // Initialize writers for stream 0
    stream.init_writers(&config, 0)?;

    // Verify the state transitioned to OutputMode::Merged
    match &stream.output_mode {
        OutputMode::Merged(merged_out) => {
            let keys: Vec<String> = merged_out
                .header()
                .read_groups()
                .keys()
                .map(|k| k.to_string())
                .collect();

            // 1 original + 2 derived suffixes = 3 total Read Groups expected
            assert_eq!(keys.len(), 3, "Header was not properly expanded");
            assert!(keys.contains(&"rg_test".to_string()));
            assert!(keys.contains(&format!("rg_test{}", crate::bam::SUFFIX_FILTERED)));
            assert!(keys.contains(&format!("rg_test{}", crate::bam::SUFFIX_AMBIGUOUS)));
        }
        OutputMode::MultiFile { .. } => {
            panic!("init_writers failed to set OutputMode::Merged when configured");
        }
    }

    // Verify the file was physically created
    assert!(merged_path.exists(), "Merged file was not created on disk");

    Ok(())
}

#[test]
fn test_init_writers_defaults_to_multi_file() -> Result<()> {
    let temp_dir = tempfile::tempdir()?;
    let out_path = temp_dir.path().join("out.bam");

    let mut config = Config::default();
    config.output = vec![out_path.clone()];
    config.no_program_line = true;

    let mut stream = empty_aln_stream();
    stream.init_writers(&config, 0)?;

    // Verify the state is MultiFile with the output writer populated
    match &stream.output_mode {
        OutputMode::MultiFile {
            output,
            filt,
            ambiguous,
        } => {
            assert!(
                output.is_some(),
                "Standard output writer should be initialized"
            );
            assert!(filt.is_none(), "Filtered output should be None");
            assert!(ambiguous.is_none(), "Ambiguous output should be None");
        }
        OutputMode::Merged(_) => {
            panic!("init_writers incorrectly fell back to Merged mode");
        }
    }

    assert!(out_path.exists());
    Ok(())
}
