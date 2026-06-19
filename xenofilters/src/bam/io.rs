//! BAM output helpers.
//!
//! The key type is [`BamOutput`], an enum that dispatches to either a
//! single-threaded bgzf writer (1 thread) or noodles'
//! [`bgzf::MultithreadedWriter`] (> 1 thread).  Both variants implement the
//! same `write_alignment_record` call so callers need no conditional logic.
//!
//! ## Why not the async writer?
//! `bgzf::r#async::writer::Builder` does expose `set_worker_count`, but the
//! async writer requires a tokio executor and returns `Future`s.  Wiring an
//! async sink into an otherwise fully synchronous main loop adds considerable
//! complexity for zero practical benefit — the crossbeam-based
//! `bgzf::MultithreadedWriter` achieves the same parallelism with plain
//! `std::thread` workers and a synchronous `Write` interface.

use anyhow::{anyhow, Result};
use noodles::bam::io::Writer as BamWriter;
use noodles::bgzf::io::MultithreadedWriter;
use noodles::bgzf::{
    self,
    io::{writer::CompressionLevel, Writer as BgzfSyncWriter},
};
use noodles::sam::{
    alignment::{io::Write as AlignmentWrite, record_buf::RecordBuf},
    header::record::value::{map::program::tag, Map},
    Header,
};
use std::{fs::File, num::NonZeroUsize, path::Path};

// -- @PG helper ----------------------------------------------------------------

fn add_pg_line(header: &mut Header) -> Result<()> {
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
/// `threads == 1` → [`bgzf::io::Writer`] (single sync writer, no extra threads).  
/// `threads  > 1` → [`bgzf::MultithreadedWriter`] (crossbeam thread pool with
///                  `threads` compression workers).
pub(crate) enum BamOutput {
    Single(BamWriter<BgzfSyncWriter<File>>),
    Multi(BamWriter<MultithreadedWriter<File>>),
}

impl BamOutput {
    /// Write one alignment record through whichever writer is active.
    pub(crate) fn write_alignment_record(
        &mut self,
        header: &Header,
        rec: &RecordBuf,
    ) -> std::io::Result<()> {
        match self {
            Self::Single(w) => w.write_alignment_record(header, rec),
            Self::Multi(w) => w.write_alignment_record(header, rec),
        }
    }
}

// -- Public constructor --------------------------------------------------------

pub(crate) fn path_unicode_ok<'a, P: 'a + AsRef<Path>>(path: P) -> Result<()> {
    path.as_ref()
        .to_str()
        .ok_or_else(|| anyhow!("Path '{}' is not valid UTF-8", path.as_ref().display()))?;
    Ok(())
}

/// Open a BAM output file with `threads` bgzf compression workers.
///
/// - `threads == 1` → single-threaded [`bgzf::io::Writer`], no extra OS
///   threads are created.
/// - `threads  > 1` → [`bgzf::MultithreadedWriter`] with a crossbeam-based
///   thread pool of `threads` workers; compression runs concurrently with the
///   main record-processing loop.
pub(crate) fn out_from_file(
    f: &Path,
    header: &Header,
    add_pg: bool,
    threads: usize,
) -> Result<BamOutput> {
    tracing::debug!(path = %f.display(), threads, "Opening BAM output file");

    let mut modified_header = header.clone();
    if add_pg {
        add_pg_line(&mut modified_header)?;
    }

    let file =
        File::create(f).map_err(|e| anyhow!("Cannot create output file '{}': {e}", f.display()))?;

    if threads <= 1 {
        // Single-threaded path: plain sync bgzf writer, no extra threads.
        let enc = bgzf::io::writer::Builder::default()
            .set_compression_level(CompressionLevel::FAST)
            .build_from_writer(file);
        let mut writer = BamWriter::from(enc);
        writer.write_header(&modified_header)?;
        Ok(BamOutput::Single(writer))
    } else {
        // Multithreaded path: bgzf::MultithreadedWriter uses crossbeam channels
        // and std::thread workers — no tokio required.
        let worker_count =
            NonZeroUsize::new(threads).expect("threads > 1 guaranteed by branch condition");
        let enc = bgzf::io::multithreaded_writer::Builder::default()
            .set_worker_count(worker_count)
            .build_from_writer(file);
        let mut writer = BamWriter::from(enc);
        writer.write_header(&modified_header)?;
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
    fn test_add_pg_line() -> Result<()> {
        let mut header = Header::default();
        add_pg_line(&mut header)?;
        let pg = header.programs();
        let mut roots = pg.roots();
        let (id, map): (&[u8], &_) = roots
            .next()
            .map(|(id, map)| (id.as_ref(), map))
            .expect("No PG record written");
        assert_eq!(id, b"xenofilter");

        let of = map.other_fields();
        let vn = of.get(&tag::VERSION).expect("VN tag missing").to_string();
        let vn_parts = vn.split('.').collect::<Vec<_>>();
        let env_vn_parts = env!("CARGO_PKG_VERSION").split('.').collect::<Vec<_>>();
        assert_eq!(
            vn_parts[0].parse::<u32>().unwrap(),
            env_vn_parts[0].parse::<u32>().unwrap(),
            "major version mismatch"
        );
        assert_eq!(
            vn_parts[1].parse::<u32>().unwrap(),
            env_vn_parts[1].parse::<u32>().unwrap(),
            "minor version mismatch"
        );
        let cl = of.get(&tag::COMMAND_LINE).expect("CL tag missing");
        assert_eq!(
            cl.to_string(),
            std::env::args().collect::<Vec<_>>().join(" ")
        );
        assert_eq!(roots.next(), None, "Unexpected second PG entry");
        Ok(())
    }
}
