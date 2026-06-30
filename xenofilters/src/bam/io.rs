//! `src/bam/io.rs`
//!
//! BAM output helpers — single-threaded and multithreaded writers,
//! plus the merged-output variant.

use crate::Error;
use noodles::bam::io::Writer as BamWriter;
use noodles::bgzf::io::{
    MultithreadedWriter, Writer as BgzfSyncWriter,
    multithreaded_writer::Builder as MultiBuilder,
    writer::{Builder, CompressionLevel},
};
use noodles::sam::{
    Header,
    alignment::{io::Write as AlignmentWrite, record_buf::RecordBuf},
    header::record::value::{Map, map::program::tag},
};
use std::io::Stdout;
use std::{fs::File, num::NonZeroUsize, path::Path};

// -- @PG helper ----------------------------------------------------------------

fn add_pg_line(header: &mut Header) -> Result<(), Error> {
    let program = Map::builder()
        .insert(tag::NAME, "xenofilter")
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(
            tag::COMMAND_LINE,
            std::env::args().collect::<Vec<_>>().join(" "),
        )
        .build()
        .expect("Failed to build PG record");
    header.programs_mut().add("xenofilter", program)?;
    Ok(())
}

// -- BamOutput enum ------------------------------------------------------------

/// A BAM writer that is either single-threaded or multithreaded.
///
/// `threads == 1` → single-threaded [`bgzf::io::Writer`].
/// `threads  > 1` → [`bgzf::MultithreadedWriter`] with a crossbeam thread
///                  pool of `threads` compression workers.
pub(crate) enum BamOutput {
    Single(BamWriter<BgzfSyncWriter<File>>),
    Multi(BamWriter<MultithreadedWriter<File>>),
    Stdout(BamWriter<BgzfSyncWriter<Stdout>>),
}

impl BamOutput {
    pub(crate) fn write_alignment_record(
        &mut self,
        header: &Header,
        rec: &RecordBuf,
    ) -> std::io::Result<()> {
        match self {
            Self::Single(w) => w.write_alignment_record(header, rec),
            Self::Multi(w) => w.write_alignment_record(header, rec),
            Self::Stdout(w) => w.write_alignment_record(header, rec),
        }
    }
}

impl Default for BamOutput {
    fn default() -> Self {
        Self::Stdout(BamWriter::from(
            Builder::default().build_from_writer(std::io::stdout()),
        ))
    }
}

// -- MergedOutput --------------------------------------------------------------

/// A single BAM file that receives winners, discarded reads, and ambiguous reads
/// from all streams, distinguished by `RG:Z` tag suffixes.
///
/// Created when `--merged-output` is supplied.  The header it was opened with
/// already contains the expanded `@RG` lines (see [`crate::bam::expand_header`]).
///
/// The `Arc` wrapper lets `LineByLine` hold one `MergedOutput` and hand
/// references to the IO helpers that need to call `write_alignment_record`.
#[derive(Default)]
pub(crate) struct MergedOutput {
    writer: BamOutput,
    /// The expanded header (with `_xenofilt` / `_xenoambig` RG groups).
    header: Header,
}

impl MergedOutput {
    pub(crate) fn new(
        path: &Path,
        header: Header, // already expanded by expand_header()
        add_pg: bool,
        threads: usize,
    ) -> Result<Self, Error> {
        let mut hdr = header;
        if add_pg {
            add_pg_line(&mut hdr)?;
        }
        let writer = open_writer(path, &hdr, threads)?;
        Ok(Self {
            writer,
            header: hdr,
        })
    }

    pub(crate) fn write_alignment_record(&mut self, rec: &RecordBuf) -> std::io::Result<()> {
        self.writer.write_alignment_record(&self.header, rec)
    }
}

// -- Public constructors -------------------------------------------------------

pub(crate) fn path_unicode_ok<'a, P: 'a + AsRef<Path>>(path: P) -> Result<(), Error> {
    path.as_ref()
        .to_str()
        .ok_or_else(|| Error::InvalidPathUtf8 {
            path: path.as_ref().to_path_buf(),
        })?;
    Ok(())
}

/// Open a per-destination BAM output file (original multi-file mode).
pub(crate) fn out_from_file(
    f: &Path,
    header: &Header,
    add_pg: bool,
    threads: usize,
) -> Result<BamOutput, Error> {
    tracing::debug!(path = %f.display(), threads, "Opening BAM output file");

    let mut modified_header = header.clone();
    if add_pg {
        add_pg_line(&mut modified_header)?;
    }
    open_writer(f, &modified_header, threads)
}

/// Shared writer construction logic.
fn open_writer(f: &Path, header: &Header, threads: usize) -> Result<BamOutput, Error> {
    let file = File::create(f).map_err(|e| Error::CreateOutputFileFailed {
        path: f.to_path_buf(),
        source: e,
    })?;

    if threads <= 1 {
        let enc = Builder::default()
            .set_compression_level(CompressionLevel::FAST)
            .build_from_writer(file);
        let mut writer = BamWriter::from(enc);
        writer.write_header(header)?;
        Ok(BamOutput::Single(writer))
    } else {
        let worker_count = NonZeroUsize::new(threads).expect("threads > 1 guaranteed by branch");
        let enc = MultiBuilder::default()
            .set_worker_count(worker_count)
            .build_from_writer(file);
        let mut writer = BamWriter::from(enc);
        writer.write_header(header)?;
        Ok(BamOutput::Multi(writer))
    }
}

// -- Tests ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_path_unicode_ok() {
        assert!(path_unicode_ok("file.bam").is_ok());
        assert!(path_unicode_ok("file with spaces.bam").is_ok());
        assert!(path_unicode_ok("file_with_üñîçødé.bam").is_ok());
    }

    #[test]
    fn test_add_pg_line() -> Result<(), Error> {
        let mut header = Header::default();
        add_pg_line(&mut header)?;
        let mut roots = header.programs().roots();
        let (id, map) = roots
            .next()
            .map(|(id, map)| (id.to_string(), map))
            .expect("No PG record written");
        assert_eq!(&id, "xenofilter");
        let of = map.other_fields();
        let vn = of.get(&tag::VERSION).expect("VN tag missing").to_string();
        let vn_parts = vn.split('.').collect::<Vec<_>>();
        let env_vn_parts = env!("CARGO_PKG_VERSION").split('.').collect::<Vec<_>>();
        assert_eq!(
            vn_parts[0].parse::<u32>().unwrap(),
            env_vn_parts[0].parse::<u32>().unwrap()
        );
        assert_eq!(
            vn_parts[1].parse::<u32>().unwrap(),
            env_vn_parts[1].parse::<u32>().unwrap()
        );
        assert_eq!(roots.next(), None);
        Ok(())
    }
}
