use crate::tests::{BamFormat, create_record};
use crate::{AlignmentStream, AlnStream, Config, StripReadSuffix};
use anyhow::Result;
use rust_htslib::bam::record::Record;
use crate::variant::VariantStoreTrait;

pub(crate) struct MockStream {
    pub(crate) reads: Vec<Record>,
    written: Vec<(Record, Option<bool>)>,
    aln_stream: AlnStream,
    i: usize,
}

impl AlignmentStream for MockStream {
    fn next_qname(&self) -> &[u8] {
        self.aln_stream.next_qname()
    }

    fn un_next(&mut self, rec: Record) -> Result<()> {
        self.un_next(rec)
    }

    fn next_rec(&mut self) -> Result<Option<Record>> {
        self.next_rec()
    }

    fn write_record(&mut self, rec: Record, is_best: Option<bool>) -> Result<()> {
        self.write_record(rec, is_best)
    }
    fn init_writers(&mut self, _opt: &Config, _i: usize) -> Result<()> {
        Ok(())
    }
    fn variant_store(&self) -> Option<&dyn VariantStoreTrait> {
        None
    }
}

impl MockStream {
    pub(crate) fn new(i: usize, reads: Vec<Record>) -> Self {
        let aln_stream = AlnStream {
            ambiguous: None,
            bam: None,
            filt: None,
            next: None,
            output: None,
            sample_variants: None,
            population_variants: None,
        };
        Self {
            reads,
            written: Vec::new(),
            aln_stream,
            i
        }
    }
    fn next_rec(&mut self) -> Result<Option<Record>> {
        if let Some(rec) = self.aln_stream.next_rec()? {
            /*
            #[cfg(test)]
            eprintln!(
                "re-nexted({}): {}",
                self.i,
                std::str::from_utf8(rec.qname()).unwrap_or("Invalid UTF-8")
            );*/
            return Ok(Some(rec));
        }
        if self.reads.is_empty() {
            return Ok(None);
        }
        let rec = self.reads.remove(0);
        self.aln_stream.un_next(rec)?;
        self.aln_stream.next_rec()
    }
    fn un_next(&mut self, rec: Record) -> Result<()> {
        eprintln!(
            "Un-next({}) read: {}",
            self.i,
            std::str::from_utf8(rec.qname()).unwrap_or("Invalid UTF-8")
        );
        self.aln_stream.un_next(rec)
    }
    fn write_record(&mut self, rec: Record, state: Option<bool>) -> Result<()> {
        self.written.push((rec, state));
        Ok(())
    }
}

/*#[test]
fn test_aln_stream_new_ok() -> Result<()> {
    let mut config = Config {
        alignment: vec!["tests/data/test_input_1_a.bam".to_string()],
        stdout_format: BamFormat::Sam,
        ..Default::default()
    };

    let aln_stream = AlnStream::new(&mut config, 0)?;
    let qname = std::str::from_utf8(aln_stream.next_qname())?;
    assert_eq!(qname, "r000");
    Ok(())
}*/

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
        assert_eq!(rec.qname(), expected.qname());
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
    assert_eq!(rec1.qname(), b"read1/1");

    mock_stream.un_next(rec1)?;

    let rec2 = mock_stream.next_rec()?.unwrap();
    assert_eq!(rec2.qname(), b"read1/1");

    let rec3 = mock_stream.next_rec()?.unwrap();
    assert_eq!(rec3.qname(), b"read2/1");

    Ok(())
}
/*#[test]
fn test_aln_stream_new_force_strip_suffix_ok() -> Result<()> {
    // None of the test BAMs have read suffixes..
    let mut config = Config {
        alignment: vec!["tests/data/test_input_1_b.bam".to_string()],
        stdout_format: BamFormat::Sam,
        strip_read_suffix: StripReadSuffix::True,
        ..Default::default()
    };

    let aln_stream = AlnStream::new(&mut config, 0)?;
    //assert!(aln_stream.bam.is_paired());
    assert_eq!(config.strip_read_suffix, StripReadSuffix::True);
    Ok(())
}
#[test]
fn test_aln_stream_new_auto_detect_true_ok() -> Result<()> {
    // None of the test BAMs have read suffixes..
    let mut config = Config {
        alignment: vec!["tests/data/test_input_1_b.bam".to_string()],
        stdout_format: BamFormat::Sam,
        strip_read_suffix: StripReadSuffix::Auto,
        ..Default::default()
    };

    let aln_stream = AlnStream::new(&mut config, 0)?;
    //assert!(aln_stream.bam.is_paired());
    assert_eq!(config.strip_read_suffix, StripReadSuffix::True);
    Ok(())
}
fn test_aln_stream_new_mismatch_strip_suffix_false_instead_of_true() {
    // None of the test BAMs have read suffixes..
    let mut config = Config {
        alignment: vec!["tests/data/test_input_1_b.bam".to_string()],
        stdout_format: BamFormat::Sam,
        strip_read_suffix: StripReadSuffix::False,
        ..Default::default()
    };

    let result = AlnStream::new(&mut config, 0);
    assert!(result.is_err());
}*/
