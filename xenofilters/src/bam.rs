mod format;
mod io;
mod merged_output;

use crate::Error;
pub(crate) use format::BamFormat;
pub(crate) use io::{BamOutput, MergedOutput, out_from_file, path_unicode_ok};
pub(crate) use merged_output::{SUFFIX_AMBIGUOUS, SUFFIX_FILTERED, expand_header, rewrite_rg};
use noodles::sam::Header;
use noodles::sam::alignment::record_buf::RecordBuf;

pub(crate) enum OutputMode {
    /// Separate files for winners / discarded / ambiguous (original behaviour).
    MultiFile {
        output: Option<BamOutput>,
        filt: Option<BamOutput>,
        ambiguous: Option<BamOutput>,
    },
    /// Single merged file; non-winners get a `RG:Z` suffix before writing.
    Merged(MergedOutput),
}

impl Default for OutputMode {
    fn default() -> Self {
        OutputMode::Merged(MergedOutput::default())
    }
}

impl OutputMode {
    pub(crate) fn write(
        &mut self,
        mut rec: RecordBuf,
        is_best: Option<bool>,
        header: &Header,
    ) -> Result<(), Error> {
        match self {
            OutputMode::MultiFile {
                output,
                filt,
                ambiguous,
            } => {
                let dest = match is_best {
                    Some(true) => output.as_mut(),
                    Some(false) => filt.as_mut(),
                    None => ambiguous.as_mut(),
                };
                if let Some(w) = dest {
                    w.write_alignment_record(header, &rec)?;
                }
            }
            OutputMode::Merged(merged) => {
                // Winners: write as-is.
                // Non-winners: rewrite RG:Z tag, then write.
                let suffix = match is_best {
                    Some(true) => None,
                    Some(false) => Some(SUFFIX_FILTERED),
                    None => Some(SUFFIX_AMBIGUOUS),
                };
                if let Some(sfx) = suffix {
                    rewrite_rg(&mut rec, sfx)?;
                }
                merged.write_alignment_record(&rec)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_output_mode_write_commits_to_file() -> Result<(), Error> {
        let temp_path =
            std::env::temp_dir().join(format!("test_bam_write_{}.bam", std::process::id()));
        let header = Header::default();

        // out_from_file signature: (path, header, add_pg, threads)
        let out = out_from_file(temp_path.as_path(), &header, false, 1)?;

        let mut mode = OutputMode::MultiFile {
            output: Some(out),
            filt: None,
            ambiguous: None,
        };

        // Write a test record
        let rec = RecordBuf::default();
        mode.write(rec, Some(true), &header)?;

        // Drop the mode to ensure the underlying writers flush their buffers to disk
        drop(mode);

        // Read the BAM file back and count the records
        let mut reader = noodles::bam::io::Reader::new(std::fs::File::open(&temp_path)?);
        reader.read_header()?;
        let count = reader.records().count();

        // Clean up the temporary file
        let _ = std::fs::remove_file(&temp_path);

        // If the mutant replaced write with Ok(()), count will be 0.
        assert_eq!(
            count, 1,
            "Mutant killed: write() did not actually write the record"
        );

        Ok(())
    }
}
