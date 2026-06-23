# xenofilters — Roadmap

This document tracks planned improvements in rough priority order.
Items marked ✓ are complete; ○ are planned; ◑ are in progress.

---

## v0.2 — Stability & coverage

- ○ **Integration test suite** — programmatic BAM generation for all edge cases
  (better-graft, better-host, tie, unmapped, paired-end mixed quality, orphan reads).
  Tracked in `tests/integration_test.rs`.
- ○ **Variant rescue — SNV complete** — single-nucleotide variants fully exercised by
  tests with synthetic BAM + VCF fixtures showing score-delta change drives a
  different decision.
- ○ **Variant rescue — indels** — extend `align_alt_to_read` NW path to handle
  multi-base insertions and deletions in the alt allele.
- ○ **Multiple ALT alleles** — `parse_sample_record` and `parse_population_record` both
  carry a `FIXME` noting that only the first ALT allele is handled. Extend to iterate
  all ALT alleles and pick the highest-scoring one per read.
- ○ **Error-type unification** — `AlignmentError`, `anyhow::Error`, and `std::io::Error`
  are all in play. Consolidate into a single `xenofilters::Error` enum via `thiserror`
  so callers can pattern-match without downcasting.

---

## v0.3 — Performance

- ○ **Progress report** — Print stats oneliner, including completion read from the BAM
  (obtainable from `File::metadata().len()` vs. current file position).  Must not
  stall the main loop; update at most once per 100 k fragments.
- ○ **bgzf read parallelism** — noodles ≥ 0.65 exposes `.set_worker_count()` on the
  bgzf reader builder.  Wire `--threads` into the reader (currently only wired into
  writers).
- ○ **Reduce `RecordBuf` allocations** — the `FromBamRecord` path in `AlnStream` clones
  every record into a heap-allocated `RecordBuf` before comparisons.  Profile whether
  keeping `noodles::bam::Record` lazy (zero-copy slice into the bgzf block) through the
  comparison phase saves meaningful wall time.
- ○ **Table-driven unit tests** — reduce duplication in `fragment_state/tests.rs` and
  `ops/tests.rs` via parameterised / table-driven patterns.

---

## v0.4 — Features

- ○ **CRAM support** — noodles has a CRAM reader; wire it behind a `--input-format` flag.
  Requires a reference FASTA path for decoding.
- ○ **Streaming stdin/stdout pipeline** — primary use case: pipe directly from `bwa mem`
  without writing intermediate BAM files to disk.  Currently blocked by the requirement
  that both streams start at the same read name: need a ring-buffer resync strategy or
  guaranteed-identical order contract from the aligner.
- ○ **Memory-mapped I/O** — for NVMe storage, `mmap` over the bgzf blocks may reduce
  syscall overhead.  Benchmark before committing; kernel page cache often wins.
- ○ **Summary statistics JSON** — emit a machine-readable JSON summary (fragment counts,
  % assigned per stream, % ambiguous) to `--stats-output`. Parseable by MultiQC.
- ○ **Multi-species (> 2 streams)** — the data structures already support N streams
  (`SmallVec<[…; 2]>`, `ARG_MAX = 4`), but the fast-path ordering logic in
  `ordering.rs` only compares `best[0]` vs `best[last]`.  Extend to an N-way tournament.
- ○ **`--help` grouping** — group CLI flags into sections (Input, Output, Scoring,
  Variants, Advanced) using `clap`'s `ArgGroup` / `next_help_heading`.

---

## v1.0 — Production

- ○ **Publish to crates.io** — once variant rescue and the test suite are complete.
- ○ **Benchmark vs XenofilteR (R) and xenome** — reproducible benchmark pipeline using
  a public PDX dataset (e.g. SRP072062) with known ground truth.
- ○ **Conda / Bioconda packaging** — recipe for `bioconda` channel.
- ○ **GitHub Actions CI** — build matrix (Linux x86_64, Linux aarch64, macOS); run
  `cargo test`, `cargo clippy -- -D warnings`, `cargo fmt --check`.

---

## Known limitations (not on the roadmap)

- **Hash-map mode for non-name-sorted input** is explicitly out of scope. Use
  `samtools sort -n` to pre-sort.
- **CRAM writing** is not planned in the near term; bgzf BAM is the standard output
  for this class of tool.
