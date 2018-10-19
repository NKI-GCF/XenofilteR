#[macro_use]
extern crate nom;
extern crate clap;
use clap::{Arg, App};
use std::io::{self, BufReader, stderr};
use std::io::prelude::*;
use std::fs::File;
use nom::digit;
use std::result::Result::*;
use std::str::FromStr;
use std::collections::{HashMap,VecDeque};
//use nom::CompareResult::Error;
// first argument must be host alignment, 2nd must be graft. Arguments may be pipes.
// Both must be name sorted (or unsorted, raw bwa ouptut, which should be in fastq order).
// graft must contain all reads present in host, host may miss some (e.g. unmapped reads)


// could also consider template length in PE.

// cargo rustc -- -Z trace-macros
// cargo +nightly clippy

const PENDING: u32 = 0x4000_0000;
const DO_WRITE: u32 = 0xc000_0000;
const SKIP_WRITE: u32 = 0x8000_0000;
const SCORE_MASK: u32 = 0x3fff_ffff;

const SAMFLAG_PAIRED: u64 = 0x1;
const SAMFLAG_PRIMARY: u64 = 0x100;

named!(uint32 <&str, u32>, map_res!(digit, FromStr::from_str));

named!(cigar<&str, u32>, do_parse!(
        nr: complete!(uint32) >>
        tp: one_of!(&b"MIDNSHPX="[..]) >>
        ({
            //println!("{}, {}", nr, tp);
            match tp {
                'S' => nr,
                'I' => nr,
                _ => 0,
            }
        })
    )
);

named!(cigar_sum<&str, u32>, fold_many1!(cigar, 0, |acc, v| {
    //println!("cigar_sum: {}, {}", acc, v);
    acc + v
}));

// TODO:
// 1: scheid per readgroup (vereist readgroup ammendment)
// 2: behandel clippings in primary/non-primary beter
//

struct SamRec {
    n: String,
    s: BufReader<File>,
    v: Vec<String>,
    pub w: Option<std::io::Stdout>,
}

impl SamRec {
    pub fn new(n: &str) -> Self {
        SamRec {
            n: n.to_string(),
            s: BufReader::new(File::open(n).expect(n)),
            v: Vec::new(),
            w: None,
        }
    }
}

fn is_coord_sorted(v: &[String]) -> bool {
    v[0].starts_with("@HD") && v[2] == "SO:coordinate\n"
}

fn write_record(r: &mut SamRec, buf: &[u8]) {
    if let Some(ref mut w) = r.w {
        if let Err(e) = w.write(buf) {
            panic!("{}", e);
        }
    }
}

fn ensure_same_reads(record: &mut Vec<SamRec>) -> bool {
    let rn = record[0].v[0].clone();
    for rec in record {
        if rec.v[0] != rn {
            return true;
        }
    }
    false
}

fn handle_headers(record: &mut Vec<SamRec>) -> bool {
    let mut hashing_required = false;
    let mut line;
    for rec in record {
         while {
            line = String::new();
            rec.v = match rec.s.read_line(&mut line) {
                Err(e) => panic!("{}", e),
                Ok(0) => panic!("No reads after {} header?", rec.n),
                Ok(_) => line.split('\t').map(|x| x.to_string()).collect(),
            };
            rec.v[0].starts_with('@')
         } {
            if !hashing_required {
                hashing_required = is_coord_sorted(&rec.v)
            }
            let buf = &rec.v.join("\t").into_bytes();
            write_record(rec, buf);
         }
    }
    hashing_required
}

/// if read is mapped returns tuple: nr-mismatches, nr-clipped
fn get_mismatches(v: &[String], unmapped_penalty: u32) -> (u32, u32) {
    match cigar_sum(&v[5]) {
        Ok((_, t)) => {
            let nm = &v[v.iter().rposition(|ref x| x.starts_with("NM:i:")).unwrap()];
            (u32::from_str(nm.get(5..).unwrap()).unwrap(), t)
        },
        Err(e) => {
            if e != nom::Err::Error(nom::Context::Code("*", nom::ErrorKind::Many1)) {
                panic!("{}, {}", e, v[5]);
            }
            (unmapped_penalty, 0)
        }
    }
}

fn read_line(s: &mut BufReader<File>, line: &mut String) -> bool {
    match s.read_line(line) {
        Err(e) => panic!("{}", e),
        Ok(i) => i > 0,
    }
}

/// requires host and graft reads to be retrieved from alignments for all readnames in the same order.
fn line_by_line_filter(record: &mut Vec<SamRec>, mm_threshold: u32, unmapped_penalty: u32) -> bool {
    let graft = record.len() - 1;
    let mut i = 0;
    let mut score: u32 = 0;
    let mut graft_lines: Vec<u8> = vec![];
    let mut last = false;

    while !last {
        let flag = record[i].v[1].to_string().parse::<u64>().unwrap();

        if (flag & SAMFLAG_PRIMARY) == 0 {
            let mism = get_mismatches(&record[i].v, unmapped_penalty);
            // one added to score, to ensure it cannot be 0, which means skip.
            score = if i != graft {
                score + mism.0 + mism.1 + 1 // mismatches in host, higher score, means more likely graft.
            } else if mism.0 <= mm_threshold {
                if score > mism.0 + mism.1 + 1 {
                    graft_lines.extend_from_slice(&record[i].v.join("\t").into_bytes());
                    score - mism.0 - mism.1 - 1
                } else {0}
            } else {0};
        } else {
            // FIXME: also write non-primaries if graft in favor?
        }

        let mut line = String::new();
        let same_readname = if read_line(&mut record[i].s, &mut line) {
            record[i].v = line.split('\t').map(|x| x.to_string()).collect();
            record[i].v[0] == record[if i == 0 {graft} else {0}].v[0]
        } else {
            if i == graft {
                last = true;
            }
            false
        };
        if i != graft && !same_readname || same_readname && i == graft {
            if i != graft {
                i += 1;
            } else {
                if score > 0 {
                    write_record(&mut record[i], &graft_lines);
                }
                graft_lines.clear();
                score = 0;
                i = 0;
            };
        }
    }
    true
}

fn hashmap_filter(record: &mut Vec<SamRec>, mm_threshold: u32, unmapped_penalty: u32) {
    let mut readscore: HashMap<String, u32> = HashMap::new();
    let mut pending: VecDeque<Vec<String>> = VecDeque::new();
    let mut line;

    for rec in record {
        let is_host = rec.w.is_none();
        let mut got_line = true; // first read is already retrieved at end of header.

        while got_line {

            let flag = rec.v[1].to_string().parse::<u64>().unwrap();

            if (flag & SAMFLAG_PRIMARY) == 0 {

                let mism = get_mismatches(&rec.v, unmapped_penalty);
                let sum = mism.0 + mism.1;
                if is_host {
                    if !readscore.contains_key(&rec.v[0]) {
                        readscore.insert(rec.v[0].clone(), sum);
                        line = String::new();
                        got_line = read_line(&mut rec.s, &mut line);
                        rec.v = line.split('\t').map(|x| x.to_string()).collect();
                        continue;
                    }
                    if let Some(x) = readscore.get_mut(&rec.v[0]) {
                        *x += sum;
                    }
                } else if let Some(x) = readscore.get_mut(&rec.v[0]) {
                    let score = *x & SCORE_MASK;
                    if *x == SKIP_WRITE || mism.0 > mm_threshold || sum > score {
                        *x = SKIP_WRITE;
                    } else if (flag & SAMFLAG_PAIRED) == 0 || (*x & PENDING) != 0 {
                        *x = DO_WRITE; // single or 2nd in pair
                    } else {
                        *x = (score - sum) | PENDING;
                    }
                } else {
                    panic!("Read {} with flag {:x} in graft but not in host?", rec.v[0], flag);
                }
            }
            if !is_host {
                pending.push_back(rec.v.drain(..).collect());
                if pending.front().map_or(false, |f| f[0] == rec.v[0]) {

                    while let Some(&score) = pending.front().map_or(None,|f| readscore.get(&f[0]))
                        .filter(|&&score| score == SKIP_WRITE || score == DO_WRITE) {
                            let f = pending.pop_front().unwrap();
                            if score == DO_WRITE {
                                write_record(rec, &f.join("\t").into_bytes());
                            }
                    }
                }
            }
            line = String::new();
            got_line = read_line(&mut rec.s, &mut line);
            rec.v = line.split('\t').map(|x| x.to_string()).collect();
        }
        if let Some(f) = pending.pop_front() {
            panic!("pending reads, front with score {:x} is:\n{}", readscore.get(&f[0]).unwrap(), &f.join("\t"));
        }
        let _ = stderr().flush();
    }
}


fn main() {
    let matches = App::new("xenofilters").version("0.01").author("Roel Kluin <r.kluin@nki.nl>")
        .about("Filter host reads from xenografts")
        .arg(
                Arg::with_name("alignments")
                    .help("Alignments, reads must be in same order, and all accounted for, all identical readnames consecutive (e.g. name sorted). Only the reads mapping best in the last alignment are written to stdout.")
                    .required(true)
                    .multiple(true)
                    .number_of_values(1),
            )
        .get_matches();

    let mut record: Vec<SamRec> = vec![];
    for f in matches.values_of("alignments").unwrap().collect::<Vec<_>>() {
        record.push(SamRec::new(&f));
    }
    let graft = record.len() - 1;
    record[graft].w = Some(io::stdout());

    let mm_threshold = 4; // applies to MM_I_human
    let unmapped_penalty = 8; // applies to paired-end unampped mate

    let mut use_hashmap = handle_headers(&mut record) && ensure_same_reads(&mut record);
    if !use_hashmap {
        eprintln!("Attempting to use line-by-line processing. Requires readnames to be retrieved consecutive and in the same order for host and graft.");
        use_hashmap = !line_by_line_filter(&mut record, mm_threshold, unmapped_penalty);
    }
    if use_hashmap {
        eprintln!("Coordinate sorted input or reads not in same order, falling back to hashmap lookup for read names.");
        let _ = stderr().flush();
        hashmap_filter(&mut record, mm_threshold, unmapped_penalty)
    }
}

