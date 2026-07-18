// src/aln_stream/open.rs
use crate::file_spec::path_for_stream;
use crate::{
    aln_stream::{AlignmentStream, AlnStream},
    config::run_config::RunConfig,
    Error,
};
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{smallvec, SmallVec};
use std::num::NonZeroUsize;
use crate::file_spec::path_from_stream;

/// Open all streams for namesorted/collated (unified reader path).
/// `bgzf_threads` applies to BAM decompression; ignored for SAM/CRAM.
pub(crate) fn open_streams_unified(
    run: &mut RunConfig,
    bgzf_threads: usize,
) -> Result<SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>, Error> {
    let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];

    // MultithreadedReader parallelises bgzf block decompression.
    // Requires threads > 1 AND a non-seeking backend (namesorted / collated).
    // HashLookup pass-2 uses seek_vpos → must use Single.
    let threads = NonZeroUsize::new(bgzf_threads).unwrap_or(NonZeroUsize::MIN);

    let n = run.alignment.len();
    for i in 0..n {
        let path = &run.alignment[i];
        let path_str = path.to_string_lossy().to_string();
        tracing::debug!(stream = i, path_str, "Opening stream");
        let positive_regions = path_for_stream(specs, 0)
            .map(|s| 


        let stream = AlnStream::<RecordBuf>::new(run, i, threads.clone(), positive_regions)?;
        aln.push(Box::new(stream));
        if i > 0 {
            if aln[i].next_qname() != aln[0].next_qname() {
                return Err(Error::InvalidInput(format!(
                    "HashLookup requires all input streams to be namesorted/collated. \
                     Stream 0 next_qname: {:?}, stream {} next_qname: {:?}",
                    aln[0].next_qname(),
                    i,
                    aln[i].next_qname()
                )));
            }
        }
    }
    Ok(aln)
}

/// Open exactly 2 streams as raw seekable BGZF BAM readers (HashLookup pass-2
/// requires seek_vpos, which the unified/CRAM-capable reader cannot provide).
pub(crate) fn open_streams_raw_bam(
    run: &mut RunConfig,
) -> Result<SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>, Error> {
    if !run.alignment.iter().all(|p| p.ends_with(".bam")) {
        return Err(Error::InvalidInput(format!(
            "hashlookup requires BAM input (CRAM/SAM seeking not yet supported; \
             see ROADMAP: CRAM index seeking): {:?}",
            run.alignment
        )));
    }
    let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];
    for i in 0..2 {
        let stream = AlnStream::<RecordBuf>::new_raw_bam(i, run)?;
        aln.push(Box::new(stream));
    }
    Ok(aln)
}
