use super::PrepareError;
use crate::alignment::UnifiedOpIterator;
use noodles::bam::record::Record;
use noodles::sam::alignment::record::data::field::Value;

pub(crate) fn stringify_record(rec: &Record) -> String {
    let qname = rec.name();
    let cigar = rec.cigar();
    let mut s = format!("{qname:?}\t{cigar:?}");
    if let Some(Ok(Value::String(md))) = rec.data().get(b"MD") {
        s.push_str(&format!("\tMD:Z:{md}"));
    }
    s.push_str(&format!("\treverse:{}", rec.flags().is_reverse_complemented()));
    //s.push_str(&format!("\tseq:{}", rec.seq().as_bytes().iter().map(|&b| b as char).collect::<String>()));
    //s.push_str(&format!("\tqual:{}", rec.qual().iter().map(|&q| (q + 33) as char).collect::<String>()));

    s
}

#[cfg_attr(test, derive(Debug))]
pub(crate) struct PreparedAlignmentPair<'a> {
    iter1: UnifiedOpIterator<'a>, // host
    iter2: UnifiedOpIterator<'a>, // graft
}

impl<'a> PreparedAlignmentPair<'a> {
    pub(crate) fn host_graft_are_perfect_match(&mut self) -> (bool, bool) {
        (self.iter1.is_single_match(), self.iter2.is_single_match())
    }
}

pub(crate) struct PreparedAlignmentPairIter<'a> {
    records1: std::slice::Iter<'a, Record>,
    records2: std::slice::Iter<'a, Record>,
}

impl<'a> PreparedAlignmentPairIter<'a> {
    #[must_use]
    pub(crate) fn new(alns1: &'a [Record], alns2: &'a [Record]) -> Self {
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
                    Some(r) if r.flags().is_secondary() => continue,
                    None => return None,
                    r => r,
                };
            }
            if read_graft.is_none() {
                read_graft = match self.records2.next() {
                    Some(r) if r.flags().is_secondary() => continue,
                    None => return None,
                    r => r,
                };
            }
            break;
        }
        let read_host = read_host.take().unwrap();
        // for supplementary alignments, seq_len may differ due to hard clipping
        let read_graft = read_graft.take().unwrap();

        let iter1 = match read_host.flags().is_unmapped() {
            true => UnifiedOpIterator::empty(),
            false => match UnifiedOpIterator::new(read_host) {
                Ok(iter) => iter,
                Err(e) => return Some(Err(e)),
            },
        };
        let iter2 = match read_graft.flags().is_unmapped() {
            true => UnifiedOpIterator::empty(),
            false => match UnifiedOpIterator::new(read_graft) {
                Ok(iter) => iter,
                Err(e) => return Some(Err(e)),
            },
        };
        Some(Ok(PreparedAlignmentPair { iter1, iter2 }))
    }
}

#[cfg(test)]
pub(crate) mod tests;
