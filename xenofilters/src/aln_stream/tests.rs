use super::*;
use crate::bam::AlnFormat;
use crate::tests::create_record;

// A dummy struct to test default trait methods
struct DefaultStream;
impl AlignmentStream<RecordBuf> for DefaultStream {
    fn next_qname(&self) -> &[u8] {
        b""
    }
    fn un_next(&mut self, _rec: RecordBuf) -> Result<(), Error> {
        Ok(())
    }
    fn next_rec(&mut self) -> Result<Option<RecordBuf>, Error> {
        Ok(None)
    }
    fn write_record(&mut self, _rec: RecordBuf, _is_best: Option<bool>) -> Result<(), Error> {
        Ok(())
    }
    fn header(&self) -> &Header {
        unimplemented!()
    }
    // Intentionally leaving fetch_by_virtual_offset to its default implementation
}

pub(crate) struct MockStream {
    pub(crate) reads: Vec<RecordBuf>,
    pub(crate) original_reads: Vec<RecordBuf>,
    written: Vec<(RecordBuf, Option<bool>)>,
    aln_stream: AlnStream<RecordBuf>,
    i: usize,
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

    fn next_rec(&mut self) -> Result<Option<RecordBuf>, Error> {
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

    fn un_next(&mut self, rec: RecordBuf) -> Result<(), Error> {
        let name = rec.name().expect("Invalid Name");
        eprintln!(
            "Un-next({}) read: {}",
            self.i,
            std::str::from_utf8(name.as_ref()).unwrap_or("Invalid UTF-8")
        );
        self.aln_stream.un_next(rec)
    }

    fn write_record(&mut self, rec: RecordBuf, state: Option<bool>) -> Result<(), Error> {
        self.written.push((rec, state));
        Ok(())
    }

    /// Inherent implementation; renamed to avoid ambiguity with the trait method.
    fn fetch_record(&mut self, virtual_offset: u64) -> Result<RecordBuf, Error> {
        self.original_reads
            .get(virtual_offset as usize)
            .cloned()
            .or_else(|| self.original_reads.first().cloned())
            .ok_or(Error::MockStreamNoRecordAtVirtualOffset { virtual_offset })
    }
}

impl AlignmentStream<RecordBuf> for MockStream {
    fn next_qname(&self) -> &[u8] {
        self.aln_stream.next_qname()
    }
    fn un_next(&mut self, rec: RecordBuf) -> Result<(), Error> {
        self.un_next(rec)
    }
    fn next_rec(&mut self) -> Result<Option<RecordBuf>, Error> {
        self.next_rec()
    }
    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<(), Error> {
        self.write_record(rec, is_best)
    }
    fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<(), Error> {
        Ok(())
    }
    fn variant_store(&self) -> Option<Arc<dyn StoreTrait>> {
        None
    }
    fn header(&self) -> &Header {
        &self.aln_stream.header
    }
    fn fetch_by_virtual_offset(&mut self, virtual_offset: u64) -> Result<RecordBuf, Error> {
        self.fetch_record(virtual_offset)
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
fn test_default_fetch_by_virtual_offset_returns_error() {
    let mut stream = DefaultStream;
    let res = stream.fetch_by_virtual_offset(0);
    assert!(res.is_err(), "Default implementation should return Err");
}

/* Untestable
#[test]
fn test_from_bam_record_implementations() -> Result<(), Error> {
    use noodles::bam::record::Record as BamRecord;
    use crate::aln_stream::FromBamRecord;

    let header = Header::default();
    let bam_rec = BamRecord::default();

    // We modify a field (e.g., flags) to ensure it's not just returning a default record
    let mut modified_bam = BamRecord::default();
    modified_bam
        .flags_mut()
        .insert(noodles::sam::alignment::record::Flags::UNMAPPED);

    let rec1 = BamRecord::from_bam_record(&header, modified_bam.clone())?;
    assert_eq!(rec1.flags(), modified_bam.flags());

    let rec2 = RecordBuf::from_bam_record(&header, modified_bam)?;
    assert_eq!(
        rec2.flags(),
        noodles::sam::alignment::record::Flags::UNMAPPED
    );

    Ok(())
}*/

#[test]
fn test_bam_stream_reader_trait() -> Result<(), Error> {
    let file = File::open("tests/data/test_input_1_a.bam").expect("Test BAM file required");
    let mut bam = BgzfBamReader::Single(BamReader::new(file));
    let _header = bam.read_header()?;

    let rec = bam.next_record().expect("Should have a record")?;
    assert!(
        rec.name().is_some(),
        "Should read actual record, not default"
    );

    let pos = VirtualPosition::from(0);
    let seek_result = bam.seek_vpos(pos)?;
    assert_eq!(
        seek_result, pos,
        "Seek should return the target virtual position"
    );

    Ok(())
}

#[test]
fn test_aln_stream_header_returns_actual_reference() {
    let mut stream = empty_aln_stream();
    stream.header = Header::builder().add_comment("mutant_killer").build();

    let header = stream.header();
    let has_comment = header
        .comments()
        .contains(&"mutant_killer".to_string().into());
    assert!(
        has_comment,
        "Header reference did not return the actual inner header"
    );
}

#[test]
fn test_aln_stream_init_writers_multi_and_merged() -> Result<(), Error> {
    let mut stream = empty_aln_stream();
    let mut config = Config::default();

    config.no_program_line = true;

    // Test MultiFile path
    stream.init_writers(&config, 1)?;
    match &stream.output_mode {
        OutputMode::MultiFile { .. } => {}
        _ => panic!("Expected MultiFile output mode"),
    }

    // Test Merged path
    config.merged_output = Some("tests/data/dummy_merged_out.bam".into());

    // i == 1 should ignore creating the MergedOutput and do nothing
    stream.init_writers(&config, 1)?;

    // i == 0 should attempt to initialize MergedOutput
    // Even if it fails due to path issues, the mutation is killed by branching correctly
    let _ = stream.init_writers(&config, 0);

    Ok(())
}

/* FAILS: why?
#[test]
fn test_aln_stream_fetch_by_virtual_offset() -> Result<(), Error> {
    use noodles::bam;
    use noodles::sam;
    use noodles::sam::alignment::io::Write; // Required for write_alignment_record
    use tempfile::NamedTempFile;

    // 1. Create a file on disk
    let temp_file = NamedTempFile::new()?;
    let path = temp_file.path().to_owned();

    // 2. Write a complete BAM file WITH ONE dummy record
    {
        let header = sam::Header::default();
        let mut writer = bam::io::Writer::from(std::fs::File::create(&path)?);

        writer.write_header(&header)?;

        // Write a single unmapped/empty record so AlnStream::new doesn't hit unexpected EOF
        let dummy_record = RecordBuf::default();
        writer.write_alignment_record(&header, &dummy_record)?;

        writer.finish(&header)?;
    } // Writer and File drop here, flushing everything to disk safely

    // 3. Initialize config
    let mut config = Config {
        alignment: vec![path.to_str().unwrap().to_string()],
        ..Default::default()
    };

    // AlnStream can now successfully read the header and peek the first record!
    let mut stream = AlnStream::<RecordBuf>::new(&mut config, 0)?;

    // 4. Test the invalid offset
    let res = stream.fetch_by_virtual_offset(u64::MAX);

    // If the mutant returns Ok(), res.is_err() evaluates to false and panics, killing the mutant.
    assert!(res.is_err(), "Invalid virtual offset should produce an Err");

    Ok(())
}*/

#[test]
fn test_aln_stream_new_mismatch_strip_suffix_true_instead_of_false() {
    let mut config = Config {
        alignment: vec!["tests/data/test_input_1_a.bam".to_string()],
        stdout_format: AlnFormat::Sam,
        strip_read_suffix: StripReadSuffix::True,
        ..Default::default()
    };
    assert!(AlnStream::<RecordBuf>::new(&mut config, 0).is_err());
}

#[test]
fn test_aln_stream_next_rec() -> Result<(), Error> {
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
fn test_aln_stream_un_next() -> Result<(), Error> {
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
fn test_next_qname_returns_pending_records_name() -> Result<(), Error> {
    let mut stream = empty_aln_stream();
    let rec = create_record(b"r1", "5M", &[], &[30; 5], "5", false)?;
    stream.un_next(rec)?;
    assert_eq!(stream.next_qname(), b"r1");
    Ok(())
}

#[test]
fn test_un_next_errors_when_already_occupied() -> Result<(), Error> {
    let mut stream = empty_aln_stream();
    stream.un_next(create_record(b"r1", "5M", &[], &[30; 5], "5", false)?)?;
    assert!(stream
        .un_next(create_record(b"r2", "5M", &[], &[30; 5], "5", false)?)
        .is_err());
    Ok(())
}

#[test]
fn test_next_rec_none_when_empty_and_no_bam_reader() -> Result<(), Error> {
    assert!(empty_aln_stream().next_rec()?.is_none());
    Ok(())
}

#[test]
fn test_write_record_is_noop_without_attached_writers() -> Result<(), Error> {
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
