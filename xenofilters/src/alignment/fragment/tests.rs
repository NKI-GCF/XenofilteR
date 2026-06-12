use crate::config::Config;
use crate::penalty::Penalty;
use crate::alignment::fragment::Fragment;
use crate::tests::create_record;
use anyhow::Result;
use noodles::sam::alignment::record_buf::RecordBuf;
use smallvec::{SmallVec, smallvec};
use crate::alignment::MdCigFlags;
use crate::alignment::fragment::READ_CT;

pub(crate) fn setup_penalties() -> Penalty {
    let c = Config::default();
    let mut p = c.to_penalties();
    p.log_likelihood_match = [0.0; 93]; // 0 log-like for match
    p.log_likelihood_mismatch = [-1.0; 93];
    p.gap_open = -2.0;
    p.gap_extend = -0.5;
    p
}

#[test]
fn test_stitched_fragment_creation() -> Result<()> {
    let record1 = create_record(b"read1", "5M3S", &[b'A'; 8], &[30; 8], "5", false)?;
    let record2 = create_record(b"read1", "4M4S", &[b'A'; 8], &[30; 8], "4", false)?;
    let flags1 = record1.flags();
    let flags2 = record2.flags();

    let records: SmallVec<[&RecordBuf; READ_CT]> = smallvec![&record1, &record2];
    let p = setup_penalties();
    let mut md_cig_flags: SmallVec<[MdCigFlags; READ_CT]> = SmallVec::new();
    md_cig_flags.push(MdCigFlags::try_from_record(&record1, &flags1)?);
    md_cig_flags.push(MdCigFlags::try_from_record(&record2, &flags2)?);
    let _stitched = Fragment::new(&p, records, md_cig_flags);
    Ok(())
}
