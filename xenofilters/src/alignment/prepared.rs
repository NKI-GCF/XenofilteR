use super::PrepareError;
use crate::alignment::{UnifiedOp, UnifiedOpIterator};
use rust_htslib::bam::record::{Aux, Record};

#[allow(dead_code)]
pub fn print_req(i: usize, rec: &Record) {
    let qname = String::from_utf8_lossy(rec.qname());
    let cigar = rec.cigar().to_string();
    eprint!("{i}:{qname}\t{cigar}");
    if let Ok(Aux::String(md)) = rec.aux(b"MD") {
        eprint!("\tMD:Z:{md}");
    }
    eprintln!();
}

pub struct PreparedAlignmentPair<'a> {
    iter1: UnifiedOpIterator<'a>, // host
    iter2: UnifiedOpIterator<'a>, // graft
    seq_len: u32,
}

impl<'a> PreparedAlignmentPair<'a> {
    pub fn are_perfect_match(&'a mut self) -> (bool, bool) {
        let first = match self.iter1.peek() {
            Some(UnifiedOp::Match(len)) => self.seq_len == *len,
            _ => false,
        };
        let second = match self.iter2.peek() {
            Some(UnifiedOp::Match(len)) => self.seq_len == *len,
            _ => false,
        };
        (first, second)
    }
}

pub struct PreparedAlignmentPairIter<'a> {
    records1: std::slice::Iter<'a, Record>,
    records2: std::slice::Iter<'a, Record>,
}

impl<'a> PreparedAlignmentPairIter<'a> {
    #[must_use]
    pub fn new(alns1: &'a [Record], alns2: &'a [Record]) -> Self {
        Self {
            records1: alns1.iter(),
            records2: alns2.iter(),
        }
    }
}

impl<'a> Iterator for PreparedAlignmentPairIter<'a> {
    type Item = Result<PreparedAlignmentPair<'a>, PrepareError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut read_host: Option<&Record> = None;
        let mut read_graft: Option<&Record> = None;
        loop {
            if read_host.is_none() {
                read_host = match self.records1.next() {
                    Some(r) if r.is_secondary() => continue,
                    None => return None,
                    r => r,
                };
            }
            if read_graft.is_none() {
                read_graft = match self.records2.next() {
                    Some(r) if r.is_secondary() => continue,
                    None => return None,
                    r => r,
                };
            }
            break;
        }
        let read_host = read_host.take().unwrap();
        let seq_len = read_host.seq_len();
        let read_graft = read_graft.take().unwrap();
        assert_eq!(seq_len, read_graft.seq_len());

        let iter1 = match read_host.is_unmapped() {
            true => UnifiedOpIterator::empty(read_host.is_reverse()),
            false => match UnifiedOpIterator::new(read_host) {
                Ok(iter) => iter,
                Err(e) => return Some(Err(e)),
            },
        };
        let iter2 = match read_graft.is_unmapped() {
            true => UnifiedOpIterator::empty(read_graft.is_reverse()),
            false => match UnifiedOpIterator::new(read_graft) {
                Ok(iter) => iter,
                Err(e) => return Some(Err(e)),
            },
        };
        Some(Ok(PreparedAlignmentPair {
            iter1,
            iter2,
            seq_len: seq_len as u32,
        }))
    }
}
