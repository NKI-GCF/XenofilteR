use rust_htslib::bam::record::{Cigar, Record, CigarStringView, Aux};
use anyhow::{Result, anyhow};

/// Our new, owned CIGAR enum
pub enum CigarOp {
    Match(u32),
    Ins(u32),
    Del(u32),
    RefSkip(u32),
    SoftClip(u32),
    HardClip(u32),
    Pad(u32),
    Equal(u32),
    Diff(u32),
    /// Our new op for supplementary jumps.
    Translocate {
        new_tid: i32,
        new_pos: i64,
        new_strand: bool,
        // We can add a score/penalty here if we want
        gap_score: f64,
    },
}

impl From<&Cigar> for CigarOp {
    fn from(cigar: &Cigar) -> Self {
        match *cigar {
            Cigar::Match(l) => CigarOp::Match(l),
            Cigar::Ins(l) => CigarOp::Ins(l),
            Cigar::Del(l) => CigarOp::Del(l),
            Cigar::RefSkip(l) => CigarOp::RefSkip(l),
            Cigar::SoftClip(l) => CigarOp::SoftClip(l),
            Cigar::HardClip(l) => CigarOp::HardClip(l),
            Cigar::Pad(l) => CigarOp::Pad(l),
            Cigar::Equal(l) => CigarOp::Equal(l),
            Cigar::Diff(l) => CigarOp::Diff(l),
        }
    }
}

fn get_read_span(cigar: CigarStringView) -> (i64, i64) {
    let mut read_start = 0;
    let mut read_end = 0;

    for op in cigar.iter() {
        match *op {
            Cigar::SoftClip(len) | Cigar::HardClip(len) => {
                if read_start == 0 {
                    read_start += len as i64;
                }
            }
            Cigar::Match(len)
            | Cigar::Ins(len)
            | Cigar::Equal(len)
            | Cigar::Diff(len) => {
                read_end += len as i64;
            }
            _ => {}
        }
    }
    (read_start, read_end)
}

fn get_md_tag(record: &Record) -> Result<String> {
    match record.aux(b"MD")? {
        Aux::String(md_bytes) => Ok(md_bytes.to_string()),
        _ => Err(anyhow!("MD tag not found or of unexpected type")),
    }
}

// A new struct to hold the combined alignment data
pub struct StitchedAlignment<'a> {
    pub seq: &'a [u8],
    pub qual: &'a [u8],
    pub tid: i32,
    pub pos: i64,
    pub strand: bool,
    pub stitched_cigar: Vec<CigarOp>, // A new, combined CIGAR
    pub stitched_md: String,        // A new, combined MD string
}

/// Pre-processes a primary alignment and its supplementary records
/// into a single, stitch-aware alignment.
pub fn stitch_alignment_segments(records: &[Record]) -> Result<StitchedAlignment<'_>> {
    let mut rec_iter = records.iter();
    let primary = rec_iter.next().ok_or(anyhow!("No primary alignment found among records"))?;

    // 2. Collect and sort all segments (primary + supplementary)
    //    by their unclipped start position *on the read*.
    let mut segments: Vec<_> = rec_iter.map(|r| {
        let (read_start, read_end) = get_read_span(r.cigar());
        (read_start, read_end, r)
    }).collect();
    segments.sort_by_key(|(start, _, _)| *start);

    // 3. Build the stitched CIGAR and MD
    let mut stitched_cigar: Vec<CigarOp> = Vec::new();
    let mut last_segment: Option<&Record> = None;

    for (read_start, read_end, record) in segments.iter() {
        // 1. Add a Translocate op if this isn't the first segment
        if let Some(last_rec) = last_segment {
            // Check for continuity
            if last_rec.tid() != record.tid() || record.pos() != *read_end {
                stitched_cigar.push(CigarOp::Translocate {
                    new_tid: record.tid(),
                    new_pos: record.pos(),
                    new_strand: !record.is_reverse(),
                    gap_score: -60.0, // Example static penalty
                });
            }
        }

        // 2. Add this segment's CIGAR ops (minus clips)
        stitched_cigar.extend(
            record.cigar().iter()
                .filter(|op| !matches!(op, Cigar::SoftClip(_) | Cigar::HardClip(_)))
                .map(CigarOp::from)
        );
        
        last_segment = Some(record);
    }
    
    let primary = last_segment.unwrap(); // Or find primary record
    // The MD string is now just the simple, unmodified MD
    let md_string = get_md_tag(primary)?; 

    Ok(StitchedAlignment {
        seq: primary.seq().encoded,
        qual: primary.qual(),
        tid: primary.tid(),
        pos: primary.pos(),
        strand: !primary.is_reverse(),
        stitched_cigar,  // <-- The new Vec<CigarOp>
        stitched_md: md_string, // <-- The simple MD
    })
}

