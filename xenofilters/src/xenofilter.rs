#[macro_use]
extern crate nom;
extern crate clap;
use clap::{Arg, App};
use std::io::{self, BufReader};
use std::io::prelude::*;
use std::fs::File;
use nom::digit;
use std::result::Result::*;
use std::str::FromStr;
//use nom::CompareResult::Error;
// first argument must be host alignment, 2nd must be graft. Arguments may be pipes.
// Both must be name sorted (or unsorted, raw bwa ouptut, which should be in fastq order).
// graft must contain all reads present in host, host may miss some (e.g. unmapped reads)


// could also consider template length in PE.

// cargo rustc -- -Z trace-macros
// cargo +nightly clippy



named!(int32 <&str, i32>, map_res!(digit, FromStr::from_str));

named!(cigar<&str, i32>, do_parse!(
        nr: complete!(int32) >>
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

named!(cigar_sum<&str, i32>, fold_many1!(cigar, 0, |acc, v| {
    //println!("cigar_sum: {}, {}", acc, v);
    acc + v
}));

struct SamRec {
    n: String,
    spec: BufReader<File>,
    v: Vec<String>,
}

impl SamRec {
    pub fn new(n: &str) -> Self {
        SamRec {
            n: n.to_string(),
            spec: BufReader::new(File::open(n).expect(n)),
            v: Vec::new(),
        }
    }
}

fn is_coord_sorted(v: &[String]) -> bool {
    v[0].starts_with("@HD") && v[2] == "SO:coordinate"
}

fn main() -> io::Result<()> {
    let matches = App::new("xenofilters").version("0.01").author("Roel Kluin <r.kluin@nki.nl>")
        .about("Filter host reads from xenografts")
        .arg(
                Arg::with_name("alignments")
                    .help("Alignments, reads must be in same order, and all accounted for, all identical readnames consecutive (e.g. name sorted). Only the reads attributed to the last alignment are written to stdout.")
                    .required(true)
                    .multiple(true)
                    .number_of_values(1),
            )
        .get_matches();

    let mut spec: Vec<SamRec> = vec![];
    for f in matches.values_of("alignments").unwrap().collect::<Vec<_>>() {
        spec.push(SamRec::new(&f));
    }

    let mm_threshold = 4; // applies to MM_I_human
    let unmapped_penalty = 8; // applies to paired-end unampped mate

    let mut i = 0;
    loop {
        let mut line = String::new();
        spec[i].v = match spec[i].spec.read_line(&mut line) {
            Err(e) => panic!("{}", e),
            Ok(0) => panic!("No reads after {} header?", spec[i].n),
            Ok(_) => line.split('\t').map(|x| x.to_string()).collect(),
        };
        if spec[i].v[0].starts_with('@') {
            assert_ne!(is_coord_sorted(&spec[i].v), true,
                "Namesorted or fastq order required. {} is coordinate sorted.", spec[i].n);
            if i != 0 {
                if let Err(e) = io::stdout().write(&spec[i].v.join("\t").into_bytes()) {
                    panic!("{}", e);
                }
            }
        } else if i != spec.len() - 1 {
            i += 1;
        } else {
            //println!("{}", line);
            break;
        }
    }
    i = 0;
    let mut score: i32 = 0;
    let mut graftlines: Vec<u8> = vec![];
    let mut last = false;
    while !last {
        let is_graft = i == spec.len() - 1;
        let flag = spec[i].v[1].to_string().parse::<u64>().unwrap();
        if (flag & 256) == 0 { //primary alignments only
            let mism = match cigar_sum(&spec[i].v[5]) {
                Ok((_, v)) => {
                    let nm = &spec[i].v[spec[i].v.iter().rposition(|ref x| x.starts_with("NM:i:")).unwrap()];
                    v + i32::from_str(nm.get(5..).unwrap()).unwrap()
                },
                Err(e) => {
                    if e != nom::Err::Error(nom::Context::Code("*", nom::ErrorKind::Many1)) {
                        panic!("{}, {}", e, spec[i].v[5]);
                    }
                    unmapped_penalty
                }
            };
            if !is_graft {
                score += mism; // more mismatches in host is better.
            } else if mism <= mm_threshold {
                score -= mism;
                if score > 0 {
                    graftlines.extend_from_slice(&spec[i].v.join("\t").into_bytes());
                }
            } else {
                score = 0;
            }
        }

        let mut line = String::new();
        let same_readname = match spec[i].spec.read_line(&mut line) {
            Err(e) => panic!("{}", e),
            Ok(0) => {
                if is_graft {
                    last = true;
                }
                false
            },
            Ok(_) => {
                spec[i].v = line.split('\t').map(|x| x.to_string()).collect();
                spec[0].v[0] == spec[1].v[0]
            },
        };
        if (!is_graft && !same_readname) || (same_readname && is_graft) {
            i = if !is_graft {
                i + 1
            } else {
                if score > 0 {
                    if let Err(e) = io::stdout().write(&graftlines) {
                        panic!("{}", e);
                    }
                }
                graftlines.clear();
                score = 0;
                0
            };
        }
    }
    Ok(())
}

