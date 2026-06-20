//! `src/aln_stream.rs`
//!
//! Alignment stream abstraction.
//!
//! Two output modes are supported:
//!
//! **Multi-file mode** (default): winners, filtered, and ambiguous reads each
//! go to separate `BamOutput` files as before.
//!
//! **Merged-output mode** (`--merged-output`): all three destinations share
//! one `MergedOutput` file.  Non-winning records have their `RG:Z` tag
//! rewritten with a suffix before writing.

use crate::alignment::SimpleRec;
use crate::bam::{
    expand_header, out_from_file, path_unicode_ok, rewrite_rg, BamOutput, MergedOutput,
    SUFFIX_AMBIGUOUS, SUFFIX_FILTERED,
};
use crate::config::{Config, StripReadSuffix};
use crate::variant::{parse_population_record, parse_sample_record, Store, StoreTrait};
use anyhow::{anyhow, ensure, Result};
use noodles::bam::{io::Reader as BamReader, record::Record};
use noodles::bgzf::io::Reader as BgzfReader;
use noodles::sam::{alignment::record_buf::RecordBuf, header::Header};
use std::fs::File;
use std::io::Read as ioRead;
use std::path::PathBuf;
use std::sync::Arc;

pub(crate) trait AlignmentStream<R: SimpleRec> {
    fn next_qname(&self) -> &[u8];
    fn un_next(&mut self, rec: R) -> Result<()>;
    fn next_rec(&mut self) -> Result<Option<R>>;
    /// Write a record.
    ///
    /// `is_best`:
    /// - `Some(true)`  → winning alignment
    /// - `Some(false)` → filtered (lost) alignment
    /// - `None`        → ambiguous
    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<()>;
    fn init_writers(&mut self, opt: &Config, i: usize) -> Result<()>;
    fn variant_store(&self) -> Option<Arc<dyn StoreTrait>>;
    fn header(&self) -> &Header;
}

// -- OutputMode ----------------------------------------------------------------

/// Selects the output routing strategy for a stream.
enum OutputMode {
    /// Separate files for winners / filtered / ambiguous (original behaviour).
    MultiFile {
        output: Option<BamOutput>,
        filt: Option<BamOutput>,
        ambiguous: Option<BamOutput>,
    },
    /// Single merged file; non-winners get a `RG:Z` suffix before writing.
    Merged(MergedOutput),
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

// -- AlnStream -----------------------------------------------------------------

pub(crate) struct AlnStream<R> {
    pub(crate) bam: Option<BamReader<BgzfReader<File>>>,
    /// One-record look-ahead buffer (supports `un_next`).
    next: Option<R>,
    output_mode: OutputMode,
    variant_store: Option<Arc<dyn StoreTrait>>,
    header: Header,
    threads: usize,
}

impl<R> AlnStream<R>
where
    R: SimpleRec + FromBamRecord,
{
    pub(crate) fn new(opt: &mut Config, i: usize) -> Result<Self> {
        let bam_str = opt.alignment[i].as_str();
        tracing::debug!(stream = i, path = bam_str, "Opening BAM reader");

        let file = File::open(bam_str)
            .map_err(|e| anyhow!("Cannot open alignment file '{bam_str}': {e}"))?;
        let mut bam = BamReader::new(file);
        let header = bam.read_header()?;

        // -- Reject coordinate-sorted input --------------------------------
        let mut header_reader = bam.header_reader();
        header_reader.read_magic_number()?;
        let mut raw_sam_header_reader = header_reader.raw_sam_header_reader()?;
        let mut buf = Vec::new();
        raw_sam_header_reader.read_to_end(&mut buf)?;

        for parts in buf
            .split(|&b| b == b'\n')
            .map(|s| s.split(|&b| b == b'\t').collect::<Vec<_>>())
        {
            ensure!(
                parts.len() < 3
                    || parts[0] != b"@HD"
                    || (parts[2] != b"SO:coordinate" && parts[2] != b"GO:reference"),
                "Input file '{bam_str}' appears to be coordinate-sorted. \
                 Re-sort with `samtools sort -n`."
            );
        }

        // -- Peek at the first record --------------------------------------
        let test_record = match bam.records().next() {
            Some(Ok(rec)) => rec,
            Some(Err(e)) => return Err(anyhow!("Error reading from '{bam_str}': {e}")),
            None => return Err(anyhow!("'{bam_str}' contains no records")),
        };

        let name = test_record
            .name()
            .ok_or_else(|| anyhow!("First record in '{bam_str}' has no read name"))?;

        // -- Auto-detect strip-suffix mode ---------------------------------
        opt.strip_read_suffix = match opt.strip_read_suffix {
            StripReadSuffix::True => {
                ensure!(
                    name.ends_with(b"/1") || name.ends_with(b"/2"),
                    "Stream {i} ('{bam_str}'): --strip-read-suffix=true but \
                     first read name has no /1 or /2 suffix."
                );
                StripReadSuffix::True
            }
            StripReadSuffix::False => {
                ensure!(
                    !name.ends_with(b"/1") && !name.ends_with(b"/2"),
                    "Stream {i} ('{bam_str}'): --strip-read-suffix=false but \
                     first read name ends with a /1 or /2 suffix."
                );
                StripReadSuffix::False
            }
            StripReadSuffix::Auto => {
                if name.ends_with(b"/1") || name.ends_with(b"/2") {
                    tracing::debug!(stream = i, "Auto-detected /1 /2 read-name suffixes");
                    StripReadSuffix::True
                } else {
                    StripReadSuffix::False
                }
            }
            StripReadSuffix::Variable => StripReadSuffix::Variable,
        };

        // -- Paired-end consistency check ----------------------------------
        opt.is_paired = if i == 0 && opt.is_paired.is_none() {
            let paired = test_record.flags().is_segmented();
            tracing::debug!(paired, "Auto-detected paired-end status from stream 0");
            Some(paired)
        } else {
            ensure!(
                opt.is_paired == Some(test_record.flags().is_segmented()),
                "Stream {i} ('{bam_str}') has different paired-end status than stream 0."
            );
            opt.is_paired
        };

        // -- Variant store -------------------------------------------------
        let variant_store: Option<Arc<dyn StoreTrait>> = {
            let sample = opt
                .sample_variants
                .get(i)
                .filter(|s| !s.is_empty())
                .map(|s| -> Result<Arc<dyn StoreTrait>> {
                    tracing::debug!(stream = i, path = s, "Loading sample variants");
                    Ok(Arc::new(Store::new(
                        &PathBuf::from(s.as_str()),
                        parse_sample_record,
                    )?))
                })
                .transpose()?;

            if sample.is_some() {
                sample
            } else {
                opt.population_variants
                    .get(i)
                    .filter(|s| !s.is_empty())
                    .map(|s| -> Result<Arc<dyn StoreTrait>> {
                        tracing::debug!(stream = i, path = s, "Loading population variants");
                        Ok(Arc::new(Store::new(
                            &PathBuf::from(s.as_str()),
                            parse_population_record,
                        )?))
                    })
                    .transpose()?
            }
        };

        // -- Validate output paths early -----------------------------------
        let next_rec = R::from_bam_record(&header, test_record)?;
        if opt.merged_output.is_none() {
            opt.output.get(i).map(path_unicode_ok).transpose()?;
            opt.filtered_output
                .get(i)
                .map(path_unicode_ok)
                .transpose()?;
            opt.ambiguous_output
                .get(i)
                .map(path_unicode_ok)
                .transpose()?;
        } else {
            path_unicode_ok(opt.merged_output.as_ref().unwrap())?;
        }

        Ok(AlnStream {
            bam: Some(bam),
            next: Some(next_rec),
            // OutputMode is set later in init_writers once all streams are open.
            output_mode: OutputMode::MultiFile {
                output: None,
                filt: None,
                ambiguous: None,
            },
            variant_store,
            header,
            threads: opt.threads,
        })
    }
}

pub(crate) trait FromBamRecord: Sized {
    fn from_bam_record(header: &Header, rec: Record) -> std::io::Result<Self>;
}

impl FromBamRecord for Record {
    fn from_bam_record(_header: &Header, rec: Record) -> std::io::Result<Self> {
        Ok(rec)
    }
}

impl FromBamRecord for RecordBuf {
    fn from_bam_record(header: &Header, rec: Record) -> std::io::Result<Self> {
        RecordBuf::try_from_alignment_record(header, &rec)
    }
}

impl<R> AlignmentStream<R> for AlnStream<R>
where
    R: SimpleRec + FromBamRecord,
{
    fn next_qname(&self) -> &[u8] {
        self.next
            .as_ref()
            .and_then(|r| r.name())
            .map_or(b"", |n| n.as_ref())
    }

    fn un_next(&mut self, rec: R) -> Result<()> {
        if self.next.is_some() {
            return Err(anyhow!(
                "un_next called while look-ahead buffer is occupied"
            ));
        }
        self.next = Some(rec);
        Ok(())
    }

    fn next_rec(&mut self) -> Result<Option<R>> {
        let header = &self.header;
        Ok(self
            .next
            .take()
            .map(Ok)
            .or_else(|| {
                self.bam.as_mut().and_then(|b| {
                    b.records()
                        .next()
                        .map(|r| r.and_then(|r| R::from_bam_record(header, r)))
                })
            })
            .transpose()?)
    }

    fn write_record(&mut self, rec: RecordBuf, is_best: Option<bool>) -> Result<()> {
        let header = &self.header;
        self.output_mode.write(rec, is_best, header)
    }

    fn init_writers(&mut self, opt: &Config, i: usize) -> Result<()> {
        let add_pg = !opt.no_program_line;
        let threads = self.threads;

        if let Some(merged_path) = &opt.merged_output {
            // -- Merged-output mode ----------------------------------------
            // The header is expanded once here; all streams share the same
            // physical file but each opens an independent MergedOutput handle.
            // In practice only one stream's init_writers creates the file;
            // subsequent streams for the same path would truncate it.
            //
            // Recommendation: use --merged-output with a single logical file;
            // the header expansion covers all streams' RG lines.
            let expanded = expand_header(self.header.clone());
            self.output_mode =
                OutputMode::Merged(MergedOutput::new(merged_path, expanded, add_pg, threads)?);
            tracing::debug!(
                stream = i,
                path   = %merged_path.display(),
                "Merged output mode enabled"
            );
        } else {
            // -- Multi-file mode -------------------------------------------
            let output = opt
                .output
                .get(i)
                .map(|f| out_from_file(f, &self.header, add_pg, threads))
                .transpose()?;
            let filt = opt
                .filtered_output
                .get(i)
                .map(|f| out_from_file(f, &self.header, add_pg, threads))
                .transpose()?;
            let ambiguous = opt
                .ambiguous_output
                .get(i)
                .map(|f| out_from_file(f, &self.header, add_pg, threads))
                .transpose()?;

            tracing::debug!(
                stream = i,
                output = opt.output.get(i).map(|p| p.display().to_string()),
                filtered = opt.filtered_output.get(i).map(|p| p.display().to_string()),
                ambiguous = opt.ambiguous_output.get(i).map(|p| p.display().to_string()),
                "Multi-file output mode"
            );
            self.output_mode = OutputMode::MultiFile {
                output,
                filt,
                ambiguous,
            };
        }
        Ok(())
    }

    fn variant_store(&self) -> Option<Arc<dyn StoreTrait>> {
        self.variant_store.clone()
    }

    fn header(&self) -> &Header {
        &self.header
    }
}

#[cfg(test)]
pub(crate) mod tests;
