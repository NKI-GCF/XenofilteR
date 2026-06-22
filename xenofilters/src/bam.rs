mod format;
mod io;
mod merged_output;

use anyhow::Result;
pub(crate) use format::BamFormat;
pub(crate) use io::{out_from_file, path_unicode_ok, BamOutput, MergedOutput};
pub(crate) use merged_output::{rewrite_rg, SUFFIX_AMBIGUOUS, SUFFIX_FILTERED};
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::Header;

pub(crate) enum OutputMode {
    /// Separate files for winners / filtered / ambiguous (original behaviour).
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
    fn write(&mut self, mut rec: RecordBuf, is_best: Option<bool>, header: &Header) -> Result<()> {
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
