//! `src/bam/reader.rs`
//!
//! [`BgzfBamReader`] -- wraps single-threaded (seekable) and multithreaded
//! bgzf BAM readers behind a uniform interface.

use noodles::bam::{io::Reader as BamReader, record::Record};
use noodles::bgzf::{self, VirtualPosition, io::MultithreadedReader};
use noodles::sam::Header;
use rayon::ThreadPool;
use std::fs::File;
use std::sync::Arc;

/// Wraps both bgzf reader variants in a zero-cost internal enum.
///
/// `Single` is seekable (required by `HashLookup` pass-2).
/// `Multi` parallelises bgzf block decompression via a dedicated rayon
/// thread pool but does not support seeking. The pool must be `install()`-ed
/// around every call that may trigger decompression (not just construction),
/// since `MultithreadedReader` spawns rayon tasks lazily on each read.
pub(crate) enum BgzfBamReader {
    Single(BamReader<bgzf::io::Reader<File>>),
    Multi {
        reader: BamReader<MultithreadedReader<File>>,
        pool: Arc<ThreadPool>,
    },
}

impl BgzfBamReader {
    pub(crate) fn read_header(&mut self) -> std::io::Result<Header> {
        match self {
            Self::Single(r) => r.read_header(),
            Self::Multi { reader, pool } => pool.install(|| reader.read_header()),
        }
    }

    /// Returns `Err(Unsupported)` for the multithreaded variant.
    /// `HashLookup` must always construct with `Single`.
    pub(crate) fn seek_vpos(&mut self, pos: VirtualPosition) -> std::io::Result<VirtualPosition> {
        match self {
            Self::Single(r) => r.get_mut().seek(pos),
            Self::Multi { .. } => Err(std::io::Error::new(
                std::io::ErrorKind::Unsupported,
                "seek_vpos unavailable on MultithreadedReader; \
                 use --matching-algorithm namesorted or collated with --threads",
            )),
        }
    }

    pub(crate) fn next_record(&mut self) -> Option<std::io::Result<Record>> {
        match self {
            Self::Single(r) => r.records().next(),
            Self::Multi { reader, pool } => pool.install(|| reader.records().next()),
        }
    }
}
