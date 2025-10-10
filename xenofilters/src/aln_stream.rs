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

/*
    pub fn process_record(&mut self, opt: &Config, rec: Record, fragment_state: &mut FragmentState) -> Result<()> {
        match fragment_state.update_state(opt, &rec, &self.header_view)? {
            FragmentState::Default(_) => panic!("Illegal state"),
            new_State => {
                self.first.push_back(rec);
                *fragment_state = new_State;
            }
        }
        Ok(())
    }

    pub fn process(&mut self, opt: &Config) -> Result<()> {
        let mut record = Record::new();
        if self.bam.read(&mut record).transpose()?.is_some() {
            self.add(opt, record);
            aln[i].process(&config)?;
        }
        
        while let Some(os) = self.algo.front().map(|r| std::str::from_utf8(r.qname())).transpose()?
                .and_then(|n| self.algo.get(n)) {
            match *os {
                FragmentState::Write => {
                    let r = self.algo.pop_front().unwrap();
                    self.out.as_mut().map(|out| out.write(&r)).transpose()?;
                },
                FragmentState::Filter => {
                    let r = self.algo.pop_front().unwrap();
                    self.filt.as_mut().map(|out| out.write(&r)).transpose()?;
                },
                _ => break, // not yet conclusive
            }
        }
        Ok(())
    }

    fn get_filter_status(&self, is_ok: bool) -> u32 {
        if is_ok {
            DO_WRITE
        } else if self.filt.is_some() {
            WRITE_FILTERED
        } else {
            SKIP_FILTERED
        }
    }

const PENDING: u32 = 0x8000_0000;
const SKIP_FILTERED: u32 = 0xffff_ffff;
const WRITE_FILTERED: u32 = 0xffff_fffe;
const DO_WRITE: u32 = 0xffff_fffd;
    fn get_mismatches(&self, opt: &Config, mut mismatches: i64) -> Result<u32> {
        if let Some(last) = self.m.back() {
            let cigar = last.cigar();
            let mut it = cigar.into_iter();

            if opt.ignore_clips_beyond_contig {
                let pos = last.pos();
                let readlen: i64 = last.seq_len().try_into()?;
                if pos < readlen {
                    // ignore start softclip if before contig
                    mismatches += match it.next() {
                        Some(Cigar::Ins(u)) => *u as i64,
                        Some(Cigar::SoftClip(u)) => pos - *u as i64,
                        _ => 0,
                    };
                } else {
                    let tid: u32 = last.tid().try_into()?;
                    let contig_len = self.bam.header().target_len(tid).ok_or_else(|| anyhow!("tid {tid} not in header"))?;

                    let beyond = pos + readlen - contig_len as i64;

                    if beyond > 0 {
                        mismatches -= beyond;
                    }
                }
            }
            for cig in it {
                match cig {
                    Cigar::SoftClip(u) | Cigar::Ins(u) => mismatches += *u as i64,
                    _ => {},
                }
            }
            Ok(mismatches.try_into()?)
        } else {
            Err(anyhow!("No read to get mismatches for"))
        }
    }

    fn is_pending(&self) -> bool {
        (self.state & PENDING) != 0 && self.state != DO_WRITE
    }

    fn get_state(&self) -> u32 {
        if USE_HASHING {
            let r = self.m.front().unwrap();
            let q_str = std::str::from_utf8(r.qname()).unwrap();
            self.h.get(q_str).map(|&u| u).unwrap_or(self.default_score)
        } else {
            self.state
        }
    }

    fn set_state(&mut self, new_state: u32) {
        if USE_HASHING {
            let r = self.m.front().unwrap();
            let q = String::from_utf8_lossy(r.qname()).to_string();
            self.h.insert(q, new_state);
        } else {
            self.state = new_state;
        }
    }

    fn process_singleton(&self, r: &Record) -> Result<bool> {
        if r.tid() != r.mtid() || r.is_unmapped() {
            // other contig, don't wait for it,
            Ok(true)
        } else {
            let record_isize = r.insert_size(); // 2nd in pair if < 0
            Ok(record_isize < 0 || r.pos() + record_isize < self.m[0].pos())
            //    - (v[9].len() as i64) + get_start_clipped::<i64>(v).unwrap_or(0)
            //    - (f[9].len() as i64) + get_start_clipped::<i64>(f).unwrap_or(0)
        }
    }*/
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
