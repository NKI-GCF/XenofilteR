# xenofilters — Roadmap

This document tracks planned improvements in rough priority order.
Items marked ✓ are complete; ○ are planned; ◑ are in progress.

---

## v0.2 — Stability & coverage

- [✓] **Three backends wired** — `namesorted`, `hashlookup`, `collated` all dispatch from
  `main.rs` with correct BED/VCF region wiring.
- [✓] **`EarlyKind` refactor** — `StreamKind::Early` carries `AllUnmapped` or `AllPerfect`
  enabling a correct 2×2 decision matrix in `resolve_fragment` without scoring.
- [✓] **Supplementary read scoring** — structural penalty
  `gap_open + non_clipped_bases × gap_extend` applied in sequential, parallel, and
  `HashLookup::score_records` paths.
- [✓] **CIGAR/MD subsumption (Tier 2.5)** — `alignment_sig` + `subsumes` check in
  `LineByLine::resolve` bypasses NW DP when one stream structurally dominates.
- [✓] **Parallel scoring pipeline** — `crossbeam_channel` bounded channels; `score_threads`
  workers each own `Scratch`; IO thread holds all writers; `Arc<dyn StoreTrait>` clones
  give workers zero-cost variant access.
- [✓] **Merged output** — `OutputMode::Merged`, `expand_header()`, `rewrite_rg()`,
  `_xenofilt`/`_xenoambig` RG suffixes.
- [✓] **`AmbiguousVcfEvaluator` trait stub** — defined in `src/filter_algorithm/ambiguous_vcf.rs`.
- [✓] **`VariantRescued` decision tag** — `XR:C:<phred>` written when winning margin came
  from `scratch.last_variant_delta > 0`.
- ○ **Integration test suite** — programmatic BAM generation for all edge cases
  (better-graft, better-host, tie, unmapped, paired-end, supplementary, orphan reads).
  Tracked in `tests/integration_test.rs`.
- ○ **Variant rescue — SNV complete** — single-nucleotide variants fully exercised by
  synthetic BAM + VCF fixtures showing score-delta change drives a different decision.
- ○ **Variant rescue — indels** — extend `align_alt_to_read` NW path to handle
  multi-base insertions and deletions.
- ○ **Multiple ALT alleles** — `parse_sample_record` and `parse_population_record` currently
  handle only the first ALT. Extend to iterate all ALTs and pick the highest-scoring one.
- ○ **Error-type unification** — consolidate `AlignmentError`, `anyhow::Error`, and
  `std::io::Error` into a single `xenofilters::Error` via `thiserror`.

---

## v0.3 — Performance

- ○ **Progress report** — one-line stats updated every 100k fragments using
  `File::metadata().len()` vs current file position. Must not stall the main loop.
- ○ **bgzf read parallelism** — wire `--threads` into the bgzf reader via
  `.set_worker_count()` (currently only wired into writers).
- ○ **Reduce `RecordBuf` allocations** — profile whether keeping `noodles::bam::Record`
  lazy (zero-copy slice into bgzf block) through the comparison phase saves wall time.
- ○ **HashLookup seek-IO thread** — offload `fetch_by_virtual_offset` calls to a dedicated
  seek thread via `crossbeam_channel` to overlap seeks with scoring work.
  Stub documented in `src/filter_algorithm/hash_lookup/stage.rs`.
- ○ **Collated worker pool** — dispatch `score_pair` to a parallel thread pool.
  Output order is not guaranteed (acceptable for Collated). Stub in `collated.rs`.
- ○ **Table-driven unit tests** — reduce duplication in `fragment_state/tests.rs` and
  `ops/tests.rs` via parameterised patterns.

---

## v0.4 — Features

- ○ **`AmbiguousVcfEvaluator` implementation** — fill in `evaluate_ambiguous_variants` to
  bypass Tier-3 NW scoring when population-specific variant evidence is conclusive.
  Requires `p_variant > 0.5` invariant enforcement.
- ○ **`--variant-rescue-bonus f64`** — config flag to scale the variant rescue contribution
  independently of the mismatch penalty.
- ○ **CRAM support** — wire a CRAM reader behind `--input-format` (requires reference FASTA).
- ○ **Streaming stdin/stdout pipeline** — pipe directly from `bwa mem` without intermediate
  BAM files. Currently blocked by the requirement that both streams share a common read
  order.
- ○ **Summary statistics JSON** — emit machine-readable JSON
  (fragment counts, % assigned, % ambiguous) to `--stats-output`. Parseable by MultiQC.
- ○ **N-stream support (> 2)** — `namesorted` and `collated` are architecturally ready
  (`SmallVec<[…; 2]>` already accepts N). `HashLookup` is flagged high-memory risk for
  N > 2 (NameTable scales as O(in-flight × streams)).
- ○ **`--help` grouping** — group CLI flags into sections (Input, Output, Scoring,
  Variants, Advanced) using clap's `next_help_heading`.
- ○ **N-way round-robin tournament** — `namesorted` currently uses a sequential
  champion-vs-challenger tournament (O(N) comparisons, correct winner guaranteed
  for strict orderings). When best[0] = best[1] but best[2] is strictly better,
  all three are declared ambiguous rather than best[2] winning. A full round-robin
  over all pairs, or sorting by NW score, would handle this correctly. Only
  affects fragments where two or more streams produce identical alignment quality.

---

## v1.0 — Production

- ○ **Publish to crates.io** — once variant rescue and the test suite are complete.
- ○ **Benchmark vs XenofilteR (R) and xenome** — reproducible benchmark using a public PDX
  dataset (e.g. SRP072062) with known ground truth.
- ○ **Conda / Bioconda packaging** — recipe for the `bioconda` channel.
- ○ **GitHub Actions CI** — build matrix (Linux x86_64, Linux aarch64, macOS);
  `cargo test`, `cargo clippy -- -D warnings`, `cargo fmt --check`.

---

## Known limitations (not on the roadmap)

- **Hash-map mode for non-name-sorted input** is explicitly out of scope for `namesorted`.
  Use `samtools sort -n` to pre-sort, or use `hashlookup`.
- **CRAM writing** is not planned in the near term; bgzf BAM is the standard output format.
- **`p_variant > 0.5` requirement** is a structural constraint of the rescue scoring formula,
  not a bug. Below 0.5 the delta is always non-positive.
