#[macro_use]
extern crate nom;
#[macro_use]
extern crate derive_new;
extern crate clap;
extern crate core;
use clap::{Arg, App};
use std::{io::{self,BufReader,prelude::*},collections::{HashMap,VecDeque},path::Path,fs::File,result::Result::*,str::FromStr};
use nom::digit;
//use nom::CompareResult::Error;
// first argument must be host alignment, 2nd must be graft. Arguments may be pipes.
// Both must be name sorted (or collate or unsorted, raw bwa ouptut, which should be in fastq order).
// graft must contain all reads present in host, host may miss some (e.g. unmapped reads)


// could also consider template length in PE.

// cargo rustc -- -Z trace-macros
// cargo +nightly clippy

const PENDING: u32 =    0x4000_0000;
const DO_WRITE: u32 =   0xc000_0000;
const SKIP_WRITE: u32 = 0x8000_0000;
const SKIP_FILTERED: u32 = 0x2000_0000;
const SCORE_MASK: u32 = 0x1fff_ffff;

const SAMFLAG_PAIRED: u32 = 1;
const SAMFLAG_UNMAPPED: u32 = 4;
const SAMFLAG_MATE_UNMAPPED: u32 = 8;
const SAMFLAG_FIRST_IN_PAIR: u32 = 64;
const SAMFLAG_NON_PRIMARY: u32 = 256;

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
    pub w: Option<Box<Write>>,
    excluded_out: Option<Box<Write>>
}

#[derive(new)]
struct XFOpt {
    mm_threshold: u32,
    unmapped_penalty: u32,
    graft_weight: u32,
    as_single: bool,
    skip_non_primary: bool
}

impl SamRec {
    pub fn new(n: &str) -> Self {
        SamRec {
            n: n.to_owned(),
            s: BufReader::new(File::open(n).unwrap_or_else(|e| panic!(e))),
            v: Vec::new(),
            w: None,
            excluded_out: None,
        }
    }
}

fn is_coord_sorted(v: &[String]) -> bool {
    v[0].starts_with("@HD") && v[2] == "SO:coordinate\n"
}

fn write_record(rw: &mut Option<Box<Write>>, buf: &[u8]) {
    if let Some(ref mut w) = rw {
        if let Err(e) = w.write(buf) {
            panic!("{}", e);
        }
    }
}

fn all_same_readname(record: &mut Vec<SamRec>) -> bool {
    let rn = record[0].v[0].clone();
    for rec in record {
        if rec.v[0] != rn {
            return false;
        }
    }
    true
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
                Ok(_) => line.split('\t').map(|x| x.to_owned()).collect(),
            };
            rec.v[0].starts_with('@')
         } {
            if !hashing_required {
                hashing_required = is_coord_sorted(&rec.v)
            }
            let buf = &rec.v.join("\t").into_bytes();
            write_record(&mut rec.w, buf);
         }
    }
    hashing_required
}

fn progress_mark(r: &mut SamRec, xf_opt: &XFOpt, mark: u32) -> u32 {

    let flag = r.v[1].to_owned().parse::<u32>().unwrap();
    let mate_unmapped = (flag & SAMFLAG_MATE_UNMAPPED) != 0;
    if (mark == SKIP_WRITE && !xf_opt.as_single) || (flag & SAMFLAG_NON_PRIMARY) != 0 {
        mark
    } else if (flag & SAMFLAG_UNMAPPED) != 0 {
        if mate_unmapped || xf_opt.as_single {
            // (hint) writeable if host, do skip if graft (pair, unlesss single-end or xf_opt.as_single).
            if r.w.is_none() {DO_WRITE} else {SKIP_WRITE}
        } else {
            mark // Unchanged: penalty for an unmapped read is already evaluated when the mapped mate occurrs.
        }
    } else {
        let nm = r.v[r.v.iter().rposition(|ref x| x.starts_with("NM:i:")).unwrap()]
            .get(5..).map(u32::from_str).unwrap().unwrap();

        // if read is mapped: nr-mismatches, nr-clipped
        match cigar_sum(&r.v[5]) {
            Ok((_, clips_indels)) => {
                let sum = nm + clips_indels + if mate_unmapped {xf_opt.unmapped_penalty} else {0};
                // one added to score, to ensure it cannot be 0, which means skip.
                if r.w.is_none() { // host
                    mark + sum // mismatches in host, higher score, means more likely graft.
                } else if (nm > xf_opt.mm_threshold) || (mark != DO_WRITE && sum > (mark & SCORE_MASK)) {
                    if r.excluded_out.is_some() {SKIP_FILTERED} else {SKIP_WRITE}
                } else if mark == DO_WRITE || (flag & SAMFLAG_PAIRED) == 0 || xf_opt.as_single || (mark & PENDING) != 0 || mate_unmapped {
                    // if mate_unmapped and here we can still write it, don't wait for mate.
                    DO_WRITE
                } else {
                    (mark - sum) | PENDING
                }
            },
            Err(e) => panic!("{}, {}", e, r.v[5]),
        }
    }
}

fn get_key(v: &[String], xf_opt: &XFOpt) -> String {
    if xf_opt.as_single {
        let flag = v[1].to_owned().parse::<u32>().unwrap();
        v[0].clone() + if (flag & SAMFLAG_FIRST_IN_PAIR) != 0 {"/1"} else {"/2"}
    } else {
        v[0].clone()
    }
}

/// return whether to continue iterating
fn lookup_and_eval_write(rw: &mut Option<Box<Write>>, filtered: &mut Option<Box<Write>>, v: &[String], opt: Option<&u32>) -> bool {
    match opt { 
        Some(&DO_WRITE) => {
            write_record(rw, &v.join("\t").into_bytes());
            true
        },
        Some(&SKIP_FILTERED) => {
            write_record(filtered, &v.join("\t").into_bytes());
            true
        }
        Some(&SKIP_WRITE) => true,
        Some(_) => false,
        _ => panic!("first time observed!\n{}", v.join("\t"))
    }
}


fn hashmap_filter(record: &mut Vec<SamRec>, xf_opt: &XFOpt) {
    let mut readscore: HashMap<String, u32> = HashMap::new();
    let mut pending: VecDeque<Vec<String>> = VecDeque::new();
    let mut line;
    let mut front: Option<String> = None;

    for mut rec in record {
        while rec.v[0] != "" {
            let flag = rec.v[1].to_owned().parse::<u32>().unwrap();
            if (flag & SAMFLAG_NON_PRIMARY) == 0 || !xf_opt.skip_non_primary {
                let k = get_key(&rec.v, xf_opt);
                {
                    let x = readscore.entry(k.clone()).or_insert(xf_opt.graft_weight);
                    if rec.w.is_some() {
                        let mut tot = (*x as u64) << 32;
                        *x = progress_mark(&mut rec, xf_opt, *x);
                        tot |= *x as u64;
                    } else {
                        *x = progress_mark(&mut rec, xf_opt, *x);
                    }
                }
                if rec.w.is_some() {
                    if xf_opt.as_single {
                        let _ = lookup_and_eval_write(&mut rec.w, &mut rec.excluded_out, &rec.v, readscore.get(&k));

                    } else {
                        pending.push_back(rec.v.drain(..).collect());
                        if front.is_none() || front == Some(k) {
                            while pending.front().map_or(false, |v| {
                                let newk = get_key(&v, xf_opt);
                                lookup_and_eval_write(&mut rec.w, &mut rec.excluded_out, &v, readscore.get(&newk))
                            }) {
                                let _ = pending.pop_front().unwrap();
                            }
                            front = pending.front().map(|f| get_key(f, xf_opt));
                        }
                    }
                }
            }
            line = String::new();
            let _ = rec.s.read_line(&mut line).map_err(|e| panic!("{}", e)).ok().unwrap();
            rec.v = line.split('\t').map(|x| x.to_owned()).collect();
        }
    }
    eprintln!("Observed {} unique readnames..", readscore.len());
    if let Some(v) = pending.pop_front() {
        let k = get_key(&v, xf_opt);
        panic!("pending reads, e.g. {} with score {:x}:\n{}", k, &readscore[&k], &v.join("\t"));
    }
}

/// requires host and graft reads to be retrieved from alignments for all readnames in the same order.
fn line_by_line_filter(record: &mut Vec<SamRec>, xf_opt: &XFOpt) {
    let mut i = 0;
    let mut mark: u32 = xf_opt.graft_weight;
    let mut graft_lines: Vec<u8> = vec![];
    let graft = record.len() - 1;

    while i <= graft {
        mark = progress_mark(&mut record[i], xf_opt, mark);
        let flag = record[i].v[1].to_owned().parse::<u32>().unwrap();
        if i == graft && mark != SKIP_WRITE && ((flag & SAMFLAG_NON_PRIMARY) == 0 || !xf_opt.skip_non_primary) {
            graft_lines.extend_from_slice(&record[i].v.join("\t").into_bytes());
        }
        let mut line = String::new();
        let nr_read = record[i].s.read_line(&mut line).map_err(|e| panic!("{}", e)).unwrap();
        if nr_read != 0 {
            record[i].v = line.split('\t').map(|x| x.to_owned()).collect();
            if record[i].v[0] == record[if i == graft {0} else {graft}].v[0] {
                if i == graft {
                    if mark == DO_WRITE {
                        write_record(&mut record[i].w, &graft_lines);
                    }
                    graft_lines.clear();
                    i = 0;
                    mark = xf_opt.graft_weight;
                }
                continue;
            } else if i == graft {
                continue;
            }
        } else {
            // the last record is tricky, we insert an empty record marked as primary to ensure
            // record[host].v[0] is not the same as record[graft].v[0], and it won't be counted
            // either.
            record[i].v = vec!["".to_owned(), "256".to_owned()];
        }
        i += 1;
    }
    if mark == DO_WRITE {
        write_record(&mut record[graft].w, &graft_lines);
    }
}


fn main() {
    let app = App::new("xenofilters");
    let matches = app.version("0.01").author("Roel Kluin <r.kluin@nki.nl>")
        .about("Filter reads n multiple alignments. Only reads mapping best in an alignment that is written remain after filtering.")
        .arg(
                Arg::with_name("mm_threshold")
                    .help("Number of mismatches allowed in the second alignment.")
                    .long("mismatch-threshold")
                    .short("m")
                    .value_name("INT")
                    .required(false)
                    .takes_value(true),
            )
        .arg(
                Arg::with_name("unmapped_penalty")
                    .help("Penalty given to unmapped reads in favor of the alternative alignment.")
                    .short("u")
                    .value_name("INT")
                    .long("unmapped-penalty")
                    .required(false)
                    .takes_value(true),
            )
        .arg(
                Arg::with_name("favor_last")
                    .help("On equal alignment score, consider read to be the last alignment")
                    .long("favor-last-alignment")
                    .short("l")
                    .required(false)
            )
        .arg(
                Arg::with_name("as_single")
                    .help("Evaluate each read in a pair on its own. assumes -p. NB: breaks pairs and no attempt to fix.")
                    .long("process-as-single")
                    .short("S")
                    .required(false)
            )
        .arg(
                Arg::with_name("skip_non_primary")
                    .help("skip non primary mappings even if the primary mapping is written.")
                    .long("skip-non-primary")
                    .short("p")
                    .required(false)
            )

        .arg(
                Arg::with_name("output")
                    .help("Write reads kept from last alignment to this file (unless specified, default is stdout).")
                    .required(false)
                    .short("o")
                    .value_name("FILE")
                    .number_of_values(1),
            )
        .arg(
                Arg::with_name("filtered_reads")
                    .help("Write mapping reads, excluded from last alignment to this file.")
                    .required(false)
                    .short("f")
                    .value_name("FILE")
                    .number_of_values(1),
            )

        .arg(
                Arg::with_name("alignments")
                    .help("Alignments, 2 at least required. If reads - readnames of alignments - are consecutive and in the same order for all alignment inputs, a low memory non-hashing strategy is adopted.")
                    .required(true)
                    .multiple(true)
                    .number_of_values(1),
            )

        .get_matches();

    let mut record: Vec<SamRec> = vec![];
    for f in matches.values_of("alignments").unwrap().collect::<Vec<_>>() {
        record.push(SamRec::new(&f));
    }
    if record.len() > 1 {
        let graft = record.len() - 1;
        record[graft].w = match matches.value_of("mm_threshold") {
            Some(ref f) => Some(Box::new(File::create(&Path::new(f)).unwrap())),
            None => Some(Box::new(io::stdout())),
        };
        if let Some(ref f) = matches.value_of("mm_threshold") {
            record[graft].excluded_out = Some(Box::new(File::create(&Path::new(f)).unwrap()));
        }

        // mm_threshold and unmapped_penalty: TODO parse from commandline
        let xf_opt = XFOpt::new(matches.value_of("mm_threshold").map_or(4, |s| s.parse::<u32>().unwrap()),
            matches.value_of("unmapped_penalty").map_or(8, |s| s.parse::<u32>().unwrap()),
            if matches.is_present("favor_last_alignment") {1} else {0},
            if matches.is_present("as_single") {true} else {false},
            if matches.is_present("skip_non_primary") || matches.is_present("as_single") {true} else {false});

        let use_hashmap = handle_headers(&mut record) || !all_same_readname(&mut record) || xf_opt.as_single;
        if use_hashmap {
            eprintln!("Coordinate sorted input or reads not in same order, falling back to hashmap lookup for read names.");
            hashmap_filter(&mut record, &xf_opt);
        } else {
            eprintln!("Attempting to use line-by-line processing. Requires readnames to be retrieved consecutive and in the same order for host and graft.");
            line_by_line_filter(&mut record, &xf_opt);
        }
    } else {
        eprintln!("At least two alignments required as input");
    }
}

