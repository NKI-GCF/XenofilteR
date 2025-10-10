use std::cmp::Ordering;

use anyhow::{anyhow, Result};
use once_cell::sync::Lazy;
use rust_htslib::bam::{record::{Aux, Cigar, Record}};
use smallvec::{SmallVec, smallvec};

//const DEL_VALUE: usize = 30; // arbitrary penalty for deletions

static IPHRED: Lazy<[f64; 93]> = Lazy::new(|| {
    let mut arr = [0.0; 93];
    for (i, item) in arr.iter_mut().enumerate() {
        *item = 10.0f64.powf(i as f64 / 10.0);
    }
    arr
});

pub struct FragmentState {
    alns: SmallVec<[Record; 2]>,
}

impl FragmentState {
    pub fn from_record(r: Record) -> Self {
        FragmentState {
            alns: smallvec![r],
        }
    }
    pub fn first_qname(&self) -> &[u8] {
        self.alns.first().map_or(b"", |r| r.qname())
    }

    pub fn add_record(&mut self, r: Record) {
        self.alns.push(r);
    }
    pub fn drain(&mut self) -> SmallVec<[Record; 2]> {
        self.alns.drain(..).collect()
    }

    /*pub fn mismatch_score(&self) -> Result<f64> {
        let mut score = 0.0_f64;
        for r in self.alns.iter() {
            if r.is_secondary() {
                continue;
            }
            // FIXME: unmapped score
            if r.is_unmapped() {
                score += IPHRED[DEL_VALUE] * r.seq().len() as f64;
                continue;
            }

            // MD:Z parsing
            let md_str = match r.aux(b"MD")? {
                Aux::String(md) => md,
                _ => return Err(anyhow!("No MD tag found")),
            };
            let mut md_iter = parse_md(md_str).into_iter().peekable();
            let mut current = md_iter.next();

            let mut read_i = 0;
            
            let qual = r.qual();
            for cig in r.cigar().iter() {
                match cig {
                    Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                        for _ in 0..*len {
                            match current {
                                Some(MdOp::Match(n)) if n > 1 => current = Some(MdOp::Match(n - 1)),
                                Some(MdOp::Match(_)) => current = md_iter.next(),
                                Some(MdOp::Mismatch) => {
                                    score += IPHRED[qual[read_i] as usize];
                                    current = md_iter.next();
                                }
                                Some(MdOp::Deletion(_)) => {
                                    score += IPHRED[DEL_VALUE]; // some default penalty
                                    current = md_iter.next();
                                }
                                None => panic!("MD string ended before read"),
                            }
                            read_i += 1;
                        }
                    }
                    // Insertion relative to reference
                    Cigar::Ins(len) => {
                        for _ in 0..*len {
                            score += IPHRED[qual[read_i] as usize];
                            read_i += 1;
                        }
                    }
                    // Deletion relative to reference
                    Cigar::Del(len) => {
                        // consume MdOp deletions (if any)
                        let mut consumed = 0u64;
                        while consumed < *len as u64 {
                            match current {
                                Some(MdOp::Deletion(ref del)) => {
                                    for _ in del {
                                        score += IPHRED[DEL_VALUE];
                                        consumed += 1;
                                    }
                                    current = md_iter.next();
                                }
                                Some(MdOp::Match(n)) if n > 1 => current = Some(MdOp::Match(n - 1)),
                                Some(MdOp::Match(_)) | Some(MdOp::Mismatch) => {
                                    current = md_iter.next();
                                }
                                None => break,
                            }
                        }
                    }
                    Cigar::SoftClip(len) | Cigar::HardClip(len) => {
                        for _ in 0..*len {
                            if *cig == Cigar::SoftClip(*len) {
                                score += IPHRED[qual[read_i] as usize];
                                read_i += 1;
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
        Ok(score)
    }*/

    pub fn score(&self) -> Result<f64> {
        let mut score = 0.0_f64;
        for r in self.alns.iter() {
            if r.is_secondary() {
                continue;
            }
            if r.is_unmapped() {
                continue;
            }
            // MD:Z parsing
            let md_str = match r.aux(b"MD")? {
                Aux::String(md) => md,
                _ => return Err(anyhow!("No MD tag found")),
            };
            let mut md_iter = parse_md(md_str).into_iter().peekable();
            let mut current = md_iter.next();

            let mut read_i = 0;

            let qual = r.qual();
            for cig in r.cigar().iter() {
                match cig {
                    Cigar::Match(len) | Cigar::Equal(len) => {
                        for _ in 0..*len {
                            match current {
                                Some(MdOp::Match(n)) if n > 1 => {
                                    score += IPHRED[qual[read_i] as usize];
                                    current = Some(MdOp::Match(n - 1))
                                },
                                Some(MdOp::Match(_)) => {
                                    score += IPHRED[qual[read_i] as usize];
                                    current = md_iter.next()
                                },
                                Some(MdOp::Mismatch) => {
                                    current = md_iter.next();
                                }
                                Some(MdOp::Deletion(_)) => {
                                    current = md_iter.next();
                                }
                                None => panic!("MD string ended before read"),
                            }
                            read_i += 1;
                        }
                    }
                    Cigar::Diff(len) => {
                        for _ in 0..*len {
                            match current {
                                Some(MdOp::Match(n)) if n > 1 => current = Some(MdOp::Match(n - 1)),
                                Some(MdOp::Match(_)) => current = md_iter.next(),
                                Some(MdOp::Mismatch) => {
                                    current = md_iter.next();
                                }
                                Some(MdOp::Deletion(_)) => {
                                    current = md_iter.next();
                                }
                                None => panic!("MD string ended before read"),
                            }
                            read_i += 1;
                        }
                    }
                    // Deletion relative to reference
                    Cigar::Del(len) => {
                        // consume MdOp deletions (if any)
                        let mut consumed = 0u64;
                        while consumed < *len as u64 {
                            match current {
                                Some(MdOp::Deletion(ref del)) => {
                                    consumed += del.len() as u64;
                                    current = md_iter.next();
                                }
                                Some(MdOp::Match(n)) if n > 1 => current = Some(MdOp::Match(n - 1)),
                                Some(MdOp::Match(_)) | Some(MdOp::Mismatch) => {
                                    current = md_iter.next();
                                }
                                None => break,
                            }
                        }
                    }
                    Cigar::Ins(len) | Cigar::SoftClip(len) | Cigar::HardClip(len) => read_i += *len as usize,
                    Cigar::Pad(_) | Cigar::RefSkip(_) => {}, // do nothing
                }
            }
        }
        Ok(score)
    }
}

impl Ord for FragmentState {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self.score(), other.score()) {
            (Ok(s), Ok(os)) => s.partial_cmp(&os).unwrap_or(Ordering::Equal),
            (Err(e), _) | (_, Err(e)) => {
                eprintln!("Error computing score for self: {}", e);
                Ordering::Equal
            },
        }
    }
}
impl PartialEq for FragmentState {
    fn eq(&self, other: &Self) -> bool {
        match (self.score(), other.score()) {
            (Ok(s), Ok(os)) => s == os,
            (Err(e), _) | (_, Err(e)) => {
                eprintln!("Error computing score for self: {}", e);
                false
            },
        }
    }
}

impl PartialOrd for FragmentState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Eq for FragmentState {}


        
// A very simple MD parser
#[derive(Debug)]
enum MdOp {
    Match(u64),
    Mismatch,
    Deletion(Vec<u8>),
}

fn parse_md(md: &str) -> Vec<MdOp> {
    let mut ops = Vec::new();
    let mut num = String::new();
    let mut chars = md.chars().peekable();

    while let Some(c) = chars.next() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            if !num.is_empty() {
                ops.push(MdOp::Match(num.parse().unwrap()));
                num.clear();
            }
            if c == '^' {
                let mut del = Vec::new();
                while let Some(&d) = chars.peek() {
                    if d.is_ascii_alphabetic() {
                        del.push(d as u8);
                        chars.next();
                    } else {
                        break;
                    }
                }
                ops.push(MdOp::Deletion(del));
            } else {
                ops.push(MdOp::Mismatch);
            }
        }
    }
    if !num.is_empty() {
        ops.push(MdOp::Match(num.parse().unwrap()));
    }
    ops
}


