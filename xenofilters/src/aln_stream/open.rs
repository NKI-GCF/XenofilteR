// src/aln_stream/open.rs
use smallvec::{smallvec, SmallVec};
use noodles::sam::alignment::record_buf::RecordBuf;
use crate::{
    aln_stream::{AlignmentStream, AlnStream},
    config::run_config::RunConfig,
    Error,
};

/// Open all streams for namesorted/collated (unified reader path).
/// `bgzf_threads` applies to BAM decompression; ignored for SAM/CRAM.
pub(crate) fn open_streams_unified(
    run: &mut RunConfig,
    bgzf_threads: usize,
) -> Result<SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]>, Error> {
    let mut aln: SmallVec<[Box<dyn AlignmentStream<RecordBuf>>; 2]> = smallvec![];
    let n = run.alignment.len();
    for i in 0..n {
        let path = &run.alignment[i];
        tracing::debug!(stream = i, path, "Opening stream");
        let stream = AlnStream::<RecordBuf>::new_unified(
            path, run.io.reference.as_deref(), bgzf_threads, i, run,
        )?;
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
        let stream = AlnStream::<RecordBuf>::new_raw_bam(&run.alignment[i], i, run)?;
        aln.push(Box::new(stream));
    }
    Ok(aln)
}
