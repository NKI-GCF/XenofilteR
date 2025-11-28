use crate::alignment::{UnifiedOp, lift_alignment_ops};
use anyhow::{Result, anyhow};
use rust_htslib::bam::record::{Aux, Cigar, CigarStringView, Record};
use smallvec::{SmallVec, smallvec};

#[allow(dead_code)]
// A new struct to hold the combined alignment data
pub struct StitchedAlignment<'a> {
    pub seq: Box<dyn Iterator<Item = u8> + 'a>,
    pub qual: Box<dyn Iterator<Item = u8> + 'a>,
    pub tid: i32,
    pub pos: i64,
    pub strand: bool,
    pub stitched_ops: SmallVec<[UnifiedOp; 1]>, // A new, combined CIGAR
}

pub fn stitch_alignment_segments<'a>(
    records: &'a [Record],
    gap_penalty: f64,
) -> Result<(Option<StitchedAlignment<'a>>, Option<StitchedAlignment<'a>>)> {
    // 1. Partition records into R1, R2, and other.
    let mut groups: SmallVec<[SmallVec<[&'a Record; 1]>; 2]> = smallvec![SmallVec::new()];

    for rec in records {
        if rec.is_secondary() {
            continue;
        } // Skip secondaries
        if rec.is_first_in_template() {
            groups[0].push(rec); // R1 group
        } else if rec.is_last_in_template() {
            if groups.len() == 1 {
                groups.push(SmallVec::new());
            }
            groups[1].push(rec); // R2 group
        }
    }

    let r1_stitched = groups
        .first()
        .map(|r1_records| stitch_one_read(r1_records, gap_penalty))
        .transpose()?;

    let r2_stitched = groups
        .get(1)
        .map(|r2_records| stitch_one_read(r2_records, gap_penalty))
        .transpose()?;

    Ok((r1_stitched, r2_stitched))
}

const fn reverse_complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'N' => b'N',
        _ => base, // Non-standard bases are returned as-is
    }
}

fn get_unclipped_read_start(cigar: CigarStringView) -> u32 {
    let mut read_start = 0;
    for op in cigar.iter() {
        match *op {
            Cigar::SoftClip(len) | Cigar::HardClip(len) => read_start += len,
            _ => break,
        }
    }
    read_start
}

pub fn stitch_one_read<'a>(
    records: &[&'a Record],
    gap_penalty: f64,
) -> Result<StitchedAlignment<'a>> {
    let anchor = records
        .iter()
        .find(|r| !r.is_supplementary())
        .ok_or_else(|| anyhow!("No primary record found for stitching"))?;

    let anchor_is_reverse = anchor.is_reverse();
    let anchor_seq = anchor.seq();

    let mut segments: SmallVec<[_; 1]> = records
        .iter()
        .map(|rec| {
            let read_start = get_unclipped_read_start(rec.cigar());
            Ok((read_start, *rec))
        })
        .collect::<Result<SmallVec<[_; 1]>>>()?;

    segments.sort_by_key(|(start, _)| *start);

    let mut stitched_ops: SmallVec<[UnifiedOp; 1]> = SmallVec::new();

    for (i, (_, record)) in segments.iter().enumerate() {
        if i > 0 {
            stitched_ops.push(UnifiedOp::Translocate {
                new_tid: record.tid(),
                new_pos: record.pos(),
                new_strand: !record.is_reverse(),
                gap_score: gap_penalty,
            });
        }

        let md = match record.aux(b"MD")? {
            Aux::String(md_bytes) => md_bytes.to_string(),
            _ => return Err(anyhow!("MD tag not found or of unexpected type")),
        };

        let mut ops = lift_alignment_ops(record.cigar(), &md)?;

        if anchor_is_reverse != record.is_reverse() {
            ops.reverse();
        }

        stitched_ops.extend(ops);
    }

    let (seq, qual): (Box<dyn Iterator<Item = u8>>, Box<dyn Iterator<Item = u8>>) =
        if anchor_is_reverse {
            (
                Box::new(
                    anchor_seq
                        .encoded
                        .iter()
                        .map(|&b| reverse_complement_base(b)),
                ),
                Box::new(anchor.qual().iter().rev().copied()),
            )
        } else {
            (
                Box::new(anchor.seq().encoded.iter().copied()),
                Box::new(anchor.qual().iter().copied()),
            )
        };

    Ok(StitchedAlignment {
        seq,
        qual,
        tid: anchor.tid(),
        pos: anchor.pos(),
        strand: !anchor_is_reverse, // fastq-order strand
        stitched_ops,
    })
}
