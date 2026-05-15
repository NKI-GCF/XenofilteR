use crate::Config;
use crate::Penalty;
use crate::alignment::fragment::Fragment;
use crate::tests::create_record;
use anyhow::Result;
use noodles::bam::record::Record;
use smallvec::{SmallVec, smallvec};

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

    let records: SmallVec<[&Record; 2]> = smallvec![&record1, &record2];
    let p = setup_penalties();
    let _stitched = Fragment::new(&p, records);
    Ok(())
}
