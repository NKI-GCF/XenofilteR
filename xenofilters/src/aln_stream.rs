use anyhow::{anyhow, ensure, Result};
use rust_htslib::bam::{self, Read, record::Record};

use crate::Config;
use crate::bam_format::{out_from_file, out_stdout};

pub struct AlnStream {
    pub bam: bam::Reader,
    next: Option<Record>,
    out: Option<bam::Writer>,
    ambiguous: Option<bam::Writer>,
    filt: Option<bam::Writer>,
}

impl AlnStream {
    pub fn new(opt: &Config, i: usize)
        -> Result<Self>
    {
        let bam_str = opt.alignment[i].as_str();
        let mut bam = bam::Reader::from_path(bam_str)?;

        let mut test_record = Record::new();
        bam.read(&mut test_record).transpose()?.ok_or_else(|| anyhow!("{bam_str} has no records"))?;

        let next = Some(test_record);
        let out = opt.output.get(i).map(|f| out_from_file(f, bam.header())).transpose()?;
        let filt = opt.filtered_output.get(i).map(|f| out_from_file(f, bam.header())).transpose()?;
        let ambiguous = opt.ambiguous_output.get(i).map(|f| out_from_file(f, bam.header())).transpose()?;

        let mut stream = AlnStream {
            bam,
            next,
            out,
            filt,
            ambiguous,
        };

        if i == 0 {
            stream.out = Some(out_stdout(stream.bam.header(), opt.stdout_format.into())?);
        }
        let parts: Vec<Vec<&[u8]>> = stream.bam.header()
            .as_bytes()
            .split(|&b| b == b'\n').map(|s| s.split(|&b| b == b'\t').collect())
            .collect();

        // Should be first line in header, but in some version does not comply to specs.
        for p in &parts {
            ensure!(p.len() < 3 || p[0] != b"@HD" || (p[2] != b"SO:coordinate" && p[2] != b"GO:reference"), "Coordinate sorted input, would require hashmap lookup.");
        }

        Ok(stream)
    }
    pub fn un_next(&mut self, rec: Record) {
        self.next = Some(rec);
    }

    pub fn next_rec(&mut self) -> Result<Option<Record>> {
        self.next.take().map(Ok).or_else(|| {
            let mut record = Record::new();
            match self.bam.read(&mut record) {
                Some(Err(e)) => Some(Err(anyhow!(e))),
                Some(Ok(())) => Some(Ok(record)),
                None => None,
            }
        }).transpose()
    }

    pub fn write_record(&mut self, rec: Record, is_best: Option<bool>) -> Result<()> {
        let out = match is_best {
            Some(true) => self.out.as_mut(),
            Some(false) => self.filt.as_mut(),
            None => self.ambiguous.as_mut(),
        };
        out.map(|out| out.write(&rec)).transpose()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn skip_filtered_write_first_bit_differs() {
        // in progress_state() this is used to skip unmapped, if in config
        assert_eq!(SKIP_FILTERED | 1, WRITE_FILTERED);
    }
}
