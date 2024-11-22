#[macro_use]
extern crate nom;
#[macro_use]
extern crate derive_new;
extern crate clap;
extern crate core;
use clap::{App, Arg};
use nom::digit;
use std::{
    cmp::min,
    collections::{HashMap, VecDeque},
    fs::File,
    io::{self, prelude::*, BufReader},
    path::Path,
    result::Result::*,
    str::FromStr,
};
//use nom::CompareResult::Error;
// first argument must be host alignment, 2nd must be graft. Arguments may be pipes.
// Both must be name sorted (or collate or unsorted, raw bwa ouptut, which should be in fastq order).
// graft must contain all reads present in host, host may miss some (e.g. unmapped reads)

// could also consider template length in PE.

// cargo rustc -- -Z trace-macros
// cargo +nightly clippy

const PENDING: u32 = 0x8000_0000;
const SKIP_WRITE: u32 = 0xffff_ffff;
const SKIP_FILTERED: u32 = 0xffff_fffe;
const DO_WRITE: u32 = 0xffff_fffd;

const SAMFLAG_PAIRED: u32 = 1;
const SAMFLAG_UNMAPPED: u32 = 4;
const SAMFLAG_MATE_UNMAPPED: u32 = 8;
const SAMFLAG_1ST_IN_PAIR: u32 = 64;
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

#[derive(new)]
struct XFOpt {
    mm_threshold: u32,
    unmapped_penalty: u32,
    graft_weight: u32,
    skip_non_primary: bool,
    ignore_clips_beyond_contig: bool,
}

struct SamRec {
    contig_len: HashMap<String, u32>,
    n: String,
    s: BufReader<File>,
    v: Vec<String>,
    pub w: Option<Box<dyn Write>>,
    coord_sorted: bool,
    excluded_out: Option<Box<dyn Write>>,
}

impl SamRec {
    pub fn new(n: &str) -> Result<Self, std::io::Error> {
        Ok(SamRec {
            contig_len: HashMap::new(),
            n: n.to_owned(),
            s: File::open(n).map(BufReader::new)?,
            v: Vec::new(),
            w: None,
            coord_sorted: false,
            excluded_out: None,
        })
    }
}

/// given a tab split header line, return true if it indicates the bam is coordinate sorted.
fn is_coord_sorted(v: &[String]) -> bool {
    v[0].starts_with("@HD") && v[2] == "SO:coordinate\n"
}

/// writes record to a stream if so configured in command line arguments
fn write_record(rw: &mut Option<Box<dyn Write>>, buf: &[u8]) {
    if let Some(ref mut w) = rw {
        if let Err(e) = w.write(buf) {
            panic!("{}", e);
        }
    }
}

/// if both alignment input streams do not start with the same readname, hashing is required
fn all_same_readname(record: &mut Vec<SamRec>) -> bool {
    let rn = record[0].v[0].clone();
    for rec in record {
        if rec.v[0] != rn {
            return false;
        }
    }
    true
}

/// proceses in both input streams the sam header, evaluates whether hashing is required.
/// first sam read record for both alignment is tab split after this function.
fn handle_headers(record: &mut Vec<SamRec>, xf_opt: &XFOpt) -> bool {
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
            if !rec.coord_sorted && is_coord_sorted(&rec.v) {
                hashing_required = true;
                rec.coord_sorted = true;
            }
            if xf_opt.ignore_clips_beyond_contig && rec.v[0] == "@SQ" {
                rec.contig_len.insert(
                    rec.v[1].chars().skip(3).collect::<String>(),
                    rec.v[2]
                        .trim()
                        .chars()
                        .skip(3)
                        .collect::<String>()
                        .parse::<u32>()
                        .unwrap(),
                );
            }
            let buf = &rec.v.join("\t").into_bytes();
            write_record(&mut rec.w, buf);
            write_record(&mut rec.excluded_out, buf);
        }
    }
    hashing_required
}

fn tot(x: u32) -> u32 {
    ((x & 0xffff_0000) >> 16) | (x & 0xffff)
}

fn sub(score: u32, sum: u32) -> u32 {
    let mut fst = score & 0xffff;
    let sub_fst = sum & 0xffff;
    let mut snd = (score & 0xffff_0000) >> 16;
    let sub_snd = (sum & 0xffff_0000) >> 16;

    if fst < sub_fst {
        fst = 0;
        snd -= min(snd, sub_snd + sub_fst - fst);
    } else if snd < sub_snd {
        snd = 0;
        fst -= min(fst, sub_fst + sub_snd - snd);
    } else {
        fst -= sub_fst;
        snd -= sub_snd;
    }
    (snd << 16) | fst
}

fn get_start_clipped<T: FromStr>(v: &[String]) -> Option<T> {
    v[5].find('S')
        .and_then(|n| v[5].split_at(n).0.parse::<T>().ok())
}

fn is_pending(score: u32) -> bool {
    (score & PENDING) != 0 && score != DO_WRITE
}

fn progress_score(r: &mut SamRec, xf_opt: &XFOpt, score: u32) -> u32 {
    let flag = r.v[1].to_owned().parse::<u32>().unwrap();
    let paired = (flag & SAMFLAG_PAIRED) != 0;
    let mate_unmapped = (flag & SAMFLAG_MATE_UNMAPPED) != 0;
    if score >= SKIP_FILTERED || (flag & SAMFLAG_NON_PRIMARY) != 0 {
        score
    } else if (flag & SAMFLAG_UNMAPPED) != 0 {
        if !paired || mate_unmapped {
            // do skip if graft unmapped pair.
            if r.w.is_none() {
                DO_WRITE
            } else {
                SKIP_WRITE
            }
        } else {
            score // Unchanged: penalty for an unmapped read is already evaluated when the mapped mate occurrs.
        }
    } else {
        let is_first_in_pair = (flag & SAMFLAG_1ST_IN_PAIR) != 0;
        let nm = r.v[r.v.iter().rposition(|x| x.starts_with("NM:i:")).unwrap()]
            .get(5..)
            .map(u32::from_str)
            .unwrap()
            .unwrap();

        // if read is mapped: nr-mismatches, nr-clipped
        match cigar_sum(&r.v[5]) {
            Ok((_, clips_indels)) => {
                // preserve per read info, sum later.
                let mut mismatches = nm + clips_indels;
                if xf_opt.ignore_clips_beyond_contig {
                    let pos = r.v[3].to_owned().parse::<u32>().unwrap() - 1;
                    if let Some(p) = get_start_clipped::<u32>(&r.v).filter(|&p| p > pos) {
                        mismatches -= p - pos;
                    } else if r.v[5].ends_with('S') {
                        let readlen = r.v[9].len() as u32;
                        let beyond_contig = r.contig_len[&r.v[2]];
                        if pos + readlen > beyond_contig {
                            mismatches -= pos + readlen - beyond_contig;
                        }
                    }
                }
                if paired {
                    if is_first_in_pair {
                        // mismatches are counted in lower half, for mate in higher
                        if mate_unmapped {
                            mismatches += xf_opt.unmapped_penalty << 16;
                        }
                    } else {
                        // mismatches in higher half, for mate in lower.
                        mismatches <<= 16;
                        if mate_unmapped {
                            mismatches += xf_opt.unmapped_penalty;
                        }
                    }
                }

                // one added to score, to ensure it cannot be 0, which means skip.
                if r.w.is_none() {
                    // host
                    score + mismatches // mismatches in host, higher score, means more likely graft.
                } else if score >= SKIP_FILTERED {
                    score
                } else if !paired || is_pending(score) || mate_unmapped {
                    if nm > xf_opt.mm_threshold
                        || tot(mismatches) + xf_opt.graft_weight > tot(score & !PENDING)
                    {
                        if r.excluded_out.is_some() {
                            SKIP_FILTERED
                        } else {
                            SKIP_WRITE
                        }
                    } else {
                        DO_WRITE
                    }
                }
                /*else if r.v[8].to_owned().parse::<i64>().unwrap() + (r.v[9].len() as i64) < 0 {
                }*/
                else if nm > xf_opt.mm_threshold {
                    if r.excluded_out.is_some() {
                        SKIP_FILTERED
                    } else {
                        SKIP_WRITE
                    }
                } else {
                    sub(score, mismatches) | PENDING
                }
            }
            Err(e) => panic!("{}, {}", e, r.v[5]),
        }
    }
}

// FIXME: allow writing to both host and graft.

/// requires host and graft reads to be retrieved from alignments for all readnames in the same order.
fn line_by_line_filter(record: &mut [SamRec], xf_opt: &XFOpt) {
    let mut i = 0;
    let default_insert = xf_opt.graft_weight | (xf_opt.graft_weight << 16);
    let mut score: u32 = default_insert;
    let mut graft_lines: Vec<u8> = vec![];
    let graft = record.len() - 1;

    while i <= graft {
        score = progress_score(&mut record[i], xf_opt, score);
        let flag = record[i].v[1].to_owned().parse::<u32>().unwrap();
        if i == graft
            && score != SKIP_WRITE
            && ((flag & SAMFLAG_NON_PRIMARY) == 0 || !xf_opt.skip_non_primary)
        {
            graft_lines.extend_from_slice(&record[i].v.join("\t").into_bytes());
        }
        let mut line = String::new();
        let nr_read = record[i]
            .s
            .read_line(&mut line)
            .map_err(|e| panic!("{}", e))
            .unwrap();
        if nr_read != 0 {
            record[i].v = line.split('\t').map(|x| x.to_owned()).collect();
            let alt = if i == graft { 0 } else { graft };
            if record[i].v[0] == record[alt].v[0] {
                if i == graft {
                    if score == DO_WRITE {
                        write_record(&mut record[i].w, &graft_lines);
                    } else if record[i].excluded_out.is_some() && score == SKIP_FILTERED {
                        write_record(&mut record[i].excluded_out, &graft_lines);
                    }
                    graft_lines.clear();
                    i = 0;
                    score = default_insert;
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
    if score == DO_WRITE {
        write_record(&mut record[graft].w, &graft_lines);
    } else if record[graft].excluded_out.is_some() && score == SKIP_FILTERED {
        write_record(&mut record[graft].excluded_out, &graft_lines);
    }
}

/// return whether to continue iterating
fn lookup_and_eval_write(rec: &mut SamRec, f: &[String], val: u32, original_score: u32) -> bool {
    match val {
        DO_WRITE => {
            write_record(&mut rec.w, &f.join("\t").into_bytes());
            true
        }
        SKIP_FILTERED => {
            write_record(&mut rec.excluded_out, f.join("\t").trim().as_bytes());
            write_record(
                &mut rec.excluded_out,
                &format!("\tXf:H:{:08X}\n", original_score).into_bytes(),
            );
            true
        }
        SKIP_WRITE => true,
        _ => false,
    }
}

fn is_surpassed(f: &[String], tid: &str, pos: i64) -> bool {
    let tlen = f[8].to_owned().parse::<i64>().unwrap();
    //    - (v[9].len() as i64) + get_start_clipped::<i64>(v).unwrap_or(0)
    //    - (f[9].len() as i64) + get_start_clipped::<i64>(f).unwrap_or(0);

    if f[2] != tid || tlen <= 0 {
        //next contig, unmapped or 2nd in pair.
        true
    } else {
        // test whether past pso + tlen
        f[3].to_owned().parse::<i64>().unwrap() + tlen < pos
    }
}

/// requires memory, but not that host and graft reads occur in the same order.
fn hashmap_filter(record: &mut [SamRec], xf_opt: &XFOpt) {
    let mut readscore: HashMap<String, u32> = HashMap::new();
    let mut pending: VecDeque<Vec<String>> = VecDeque::new();
    let mut line;
    let default_insert = xf_opt.graft_weight | (xf_opt.graft_weight << 16);

    for rec in record.iter_mut() {
        while !rec.v[0].is_empty() {
            let flag = rec.v[1].to_owned().parse::<u32>().unwrap();
            if (flag & SAMFLAG_NON_PRIMARY) == 0 || !xf_opt.skip_non_primary {
                {
                    let def = default_insert;
                    let x = readscore.entry(rec.v[0].clone()).or_insert(def);
                    *x = progress_score(rec, xf_opt, *x);
                }
                if rec.w.is_some() {
                    let mut is_partial = false;
                    let tid = &rec.v[2].to_string();
                    let pos = rec.v[3].to_owned().parse::<i64>().unwrap();
                    let do_process = if let Some(f) = pending.front() {
                        if f[0] == rec.v[0] {
                            true
                        } else if rec.coord_sorted && is_surpassed(f, tid, pos) {
                            //assert!(false, "{}\n{}", f.join("\t"), rec.v.join("\t"));
                            is_partial = true;
                            true
                        } else {
                            false
                        }
                    } else {
                        true
                    };
                    pending.push_back(rec.v.drain(..).collect());
                    if do_process {
                        while pending.front().map_or(false, |f| {
                            let mut score = readscore[&f[0]];
                            let originalscore = score;
                            if score < SKIP_FILTERED && is_partial {
                                // mismatches were already subtracted in progress_score()
                                let f_flag = f[1].to_owned().parse::<u32>().unwrap();
                                let is_first_in_pair = (f_flag & SAMFLAG_1ST_IN_PAIR) != 0;
                                if score == DO_WRITE
                                    || (score
                                        & if is_first_in_pair {
                                            0xffff
                                        } else {
                                            0x7fff_0000
                                        })
                                        > 0
                                {
                                    score = DO_WRITE
                                } else {
                                    score = if rec.excluded_out.is_some() {
                                        SKIP_FILTERED
                                    } else {
                                        SKIP_WRITE
                                    }
                                }
                            }
                            lookup_and_eval_write(rec, f, score, originalscore)
                        }) {
                            let _ = pending.pop_front().unwrap();
                            if let Some(f) = pending.front() {
                                is_partial = rec.coord_sorted && is_surpassed(f, tid, pos);
                            }
                        }
                    }
                }
            }
            line = String::new();
            let _ = rec
                .s
                .read_line(&mut line)
                .map_err(|e| panic!("{}", e))
                .ok()
                .unwrap();
            rec.v = line.split('\t').map(|x| x.to_owned()).collect();
        }
    }
    if !pending.is_empty() {
        eprintln!("processing {} pending reads..", pending.len());
        let graft = record.len() - 1;
        while pending.front().map_or(false, |f| {
            let mut score = readscore[&f[0]];
            let originalscore = score;
            if is_pending(score) {
                // mismatches were already subtracted in progress_score()
                let f_flag = f[1].to_owned().parse::<u32>().unwrap();
                let is_first_in_pair = (f_flag & SAMFLAG_1ST_IN_PAIR) != 0;
                if (score
                    & if is_first_in_pair {
                        0xffff
                    } else {
                        0x7fff_0000
                    })
                    > 0
                {
                    score = DO_WRITE
                } else {
                    score = if record[graft].excluded_out.is_some() {
                        SKIP_FILTERED
                    } else {
                        SKIP_WRITE
                    }
                }
            }
            lookup_and_eval_write(&mut record[graft], f, score, originalscore)
        }) {
            let _ = pending.pop_front().unwrap();
        }
    }
    eprintln!("Observed {} unique readnames..", readscore.len());
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let app = App::new("xenofilters");
    let matches = app.version("0.01").author("Roel Kluin <r.kluin@nki.nl>")
        .about("Filter reads n multiple alignments. Reads mapping best in an alignment are piped as configured.")
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
                    .help("Penalty given to unmapped reads in favor of the alternative alignment. Set to 0 to disable this.")
                    .short("u")
                    .value_name("INT")
                    .long("unmapped-penalty")
                    .required(false)
                    .takes_value(true),
            )
        .arg(
                Arg::with_name("favor_last")
                    .help("On equal alignment score, consider read to be the last alignment.")
                    .long("favor-last-alignment")
                    .short("l")
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
                    .long("output")
                    .short("o")
                    .value_name("FILE")
                    .number_of_values(1),
            )
        .arg(
                Arg::with_name("filtered_reads")
                    .help("Write mapping reads, excluded from last alignment to this file.")
                    .required(false)
                    .long("filtered-reads")
                    .short("f")
                    .value_name("FILE")
                    .number_of_values(1),
            )
        .arg(
                Arg::with_name("use_hashing")
                    .help("Enforce hashing algorithm.")
                    .required(false)
                    .long("use-hashing")
                    .short("H")
            )
        .arg(
                Arg::with_name("ignore_clips_beyond_contig")
                    .help("Ignore clips in reads before start of contig.")
                    .required(false)
                    .long("ignore-clips-beyond-contig")
                    .short("C")
            )
        .arg(
                Arg::with_name("alignments")
                    .help("Alignments, 2 at least required. If reads - readnames of alignments - are consecutive and in the same order for all alignment inputs, a low memory non-hashing strategy is adopted.")
                    .required(true)
                    .multiple(true)
                    .number_of_values(1),
            )

        .get_matches();
    // hla: --favor-last-alignment --ignore-clips-beyond-contig -m 5 -u 0

    let mut record: Vec<SamRec> = vec![];
    for f in matches.values_of("alignments").unwrap().collect::<Vec<_>>() {
        record.push(SamRec::new(f)?);
    }
    if record.len() < 2 {
        return Err("At least two alignments required as input".into());
    }
    let graft = record.len() - 1;
    record[graft].w = match matches.value_of("output") {
        Some(ref f) => File::create(Path::new(f))
            .map(|f| Box::new(f) as Box<dyn Write>)
            .ok(),
        None => Some(Box::new(io::stdout())),
    };
    if let Some(ref f) = matches.value_of("filtered_reads") {
        record[graft].excluded_out = File::create(Path::new(f))
            .map(|f| Box::new(f) as Box<dyn Write>)
            .ok();
    }

    let xf_opt = XFOpt::new(
        matches
            .value_of("mm_threshold")
            .map_or(4, |s| s.parse::<u32>().unwrap()),
        matches
            .value_of("unmapped_penalty")
            .map_or(8, |s| s.parse::<u32>().unwrap()),
        if matches.is_present("favor_last") {
            1
        } else {
            0
        },
        matches.is_present("skip_non_primary"),
        matches.is_present("ignore_clips_beyond_contig"),
    );

    let use_hashmap = handle_headers(&mut record, &xf_opt)
        || !all_same_readname(&mut record)
        || matches.is_present("use_hashing");
    if use_hashmap {
        eprintln!("Coordinate sorted input or reads not in same order, falling back to hashmap lookup for read names.");
        hashmap_filter(&mut record, &xf_opt);
    } else {
        eprintln!("Attempting to use line-by-line processing. Requires readnames to be retrieved consecutive and in the same order for host and graft.");
        line_by_line_filter(&mut record, &xf_opt);
    }
    Ok(())
}
