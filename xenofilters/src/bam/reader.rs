//! `src/bam/reader.rs`
//!
//! [`BgzfBamReader`] -- wraps single-threaded (seekable) and multithreaded
//! bgzf BAM readers behind a uniform interface.

use noodles::bam::{io::Reader as BamReader, record::Record};
use noodles::bgzf::{self, io::MultithreadedReader, VirtualPosition};
use noodles::sam::Header;
use std::fs::File;
use std::io::Read as ioRead;

/// Wraps both bgzf reader variants in a zero-cost internal enum.
///
/// `Single` is seekable (required by `HashLookup` pass-2).
/// `Multi` parallelises bgzf block decompression but does not support seeking.
pub(crate) enum BgzfBamReader {
    Single(BamReader<bgzf::io::Reader<File>>),
    Multi(BamReader<MultithreadedReader<File>>),
}

impl BgzfBamReader {
    pub(crate) fn read_header(&mut self) -> std::io::Result<Header> {
        match self {
            Self::Single(r) => r.read_header(),
            Self::Multi(r)  => r.read_header(),
        }
    }

    /// Returns `Err(Unsupported)` for the multithreaded variant.
    /// `HashLookup` must always construct with `Single`.
    pub(crate) fn seek_vpos(&mut self, pos: VirtualPosition) -> std::io::Result<VirtualPosition> {
        match self {
            Self::Single(r) => r.get_mut().seek(pos),
            Self::Multi(_)  => Err(std::io::Error::new(
                std::io::ErrorKind::Unsupported,
                "seek_vpos unavailable on MultithreadedReader; \
                 use --matching-algorithm namesorted or collated with --threads",
            )),
        }
    }

    pub(crate) fn next_record(&mut self) -> Option<std::io::Result<Record>> {
        match self {
            Self::Single(r) => r.records().next(),
            Self::Multi(r)  => r.records().next(),
        }
    }

    /// Re-read the raw SAM header bytes for sort-order validation.
    /// Called during `AlnStream::new()` before any records are consumed.
    pub(crate) fn read_raw_header_bytes(&mut self) -> std::io::Result<Vec<u8>> {
        fn read_raw<R: ioRead>(r: &mut BamReader<R>) -> std::io::Result<Vec<u8>> {
            let mut hr  = r.header_reader();
            hr.read_magic_number()?;
            let mut rhr = hr.raw_sam_header_reader()?;
            let mut buf = Vec::new();
            rhr.read_to_end(&mut buf)?;
            Ok(buf)
        }
        match self {
            Self::Single(r) => read_raw(r),
            Self::Multi(r)  => read_raw(r),
        }
    }
}
