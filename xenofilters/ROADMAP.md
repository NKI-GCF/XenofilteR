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
- [✓] **Error-type unification** — consolidate `AlignmentError`, `anyhow::Error`, and
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
- [✓] **N-way round-robin tournament** (`namesorted`) — single scan pass per
  stream builds MCFs (borrowing from `best`), derives Tier 2/2.5/3 metrics into
  fixed-size stack arrays, then eliminates losers in O(N) backward sweeps.  No
  heap allocations; correctness no longer depends on stream ordering within `best`.

---
## v0.5 — Performance

### memchr-accelerated MD parsing
`memchr::memchr3(b'A', b'C', b'T', md)` + `memchr::memchr(b'^', md)` can replace the
hand-rolled loop in `md_mismatches`. Benefit is SIMD-width dependent and may be
negligible on typical 20–30 byte MD strings. Requires criterion benchmark showing
>5% improvement before adoption.

### jemalloc global allocator
`tikv-jemallocator` as `#[global_allocator]` benefits HashLookup's
`Box<MappedRecord>` churn. Unknown magnitude without profiling. Benchmark with
`heaptrack` against a real PDX dataset (>10 M fragments) before adopting.

### crossbeam work-stealing for parallel scoring
Replace bounded `crossbeam_channel` in `parallel_io_loop` with
`crossbeam::deque::{Injector, Worker, Stealer}`. Eliminates head-of-line blocking
when one fragment triggers deep variant rescue. Requires restructuring
`FragmentBundle` dispatch and result collection. Significant refactor.

### compact_str for short read names
Replace `Box<[u8]>` keys in `FragmentTable` with `compact_str::CompactString`
(inline up to 24 bytes). Real Illumina names are 50+ bytes; benefit likely zero
on production data. Profile first.

## v0.6 — Internal Architecture

### ScoreCtx: shared scoring context struct

**Goal:** Eliminate `penalties`, `scratch`, `add_decision_tag`, and
`ambiguous_log_threshold` being independently defined in `CollatedMatcher`,
`HashLookup`, and `LineByLine`.

**Prerequisite:** Benchmark all three backends on a representative PDX dataset
before and after, using `criterion` or wall-clock on ≥10 M fragment inputs.
The layout change moves `Scratch` (containing `SmallVec` NW buffers) behind
a struct boundary, which may affect cache line alignment and compiler inlining.

**Implementation steps:**

1. Define in `src/filter_algorithm/shared.rs`:
```rust
   pub(crate) struct ScoreCtx {
       pub(crate) penalties:               Penalty,
       pub(crate) scratch:                 Scratch,
       pub(crate) add_decision_tag:        bool,
       pub(crate) ambiguous_log_threshold: f64,
   }
```
2. Replace fields in all three backends with `ctx: ScoreCtx`.
3. Update all call sites (`score_candidate`, `nw_score_fragment`, etc.)
   to pass `&mut self.ctx`.
4. Validate: `cargo test --all` must pass; benchmark must show < 2% regression.

**Deferred because:** Field access through an extra struct indirection (`self.ctx.scratch`)
may inhibit `#[inline]` propagation across module boundaries in debug builds and
could affect LTO decisions in release. Must measure before committing.

---

### SimpleRec → correctly designed record abstraction

**Goal:** Consolidate the ad-hoc `quality_at` / `ref_seq_id` / `as_record_buf`
implementations across `bam::Record` and `RecordBuf`.

**Current blocker:** The proposed `as_record_buf(&self) -> &[u8]` in the review
notes is incorrect. The correct signature requires a `&Header` parameter and
returns `Result<RecordBuf, io::Error>`. The existing `SimpleRec` trait already
captures this correctly; what is missing is a home file, not a redesign.

**Implementation steps:**

1. Move `SimpleRec` trait + two `impl` blocks to `src/alignment/fragment/record.rs`
   (already done in this PR as a simple refactor).
2. Evaluate whether `CramStream` records require a third `impl` once CRAM support
   lands. If so, extend then; do not pre-abstract.

**Deferred because:** CRAM support is in-flight; premature abstraction before the
third implementor exists risks designing the wrong interface.

---

### Error hierarchy

**Goal:** Categorize the 40+ flat error variants for navigability.

**Rejected approach:** Sub-enums with `From` conversions — adds per-`?` overhead
in tight loops.

**Accepted approach (deferred):** Use `#[non_exhaustive]` category marker enums
that are NOT used as wrapping variants, only for `match` grouping in error
display logic. Implement a `Display` on the flat `Error` that prefixes the
category: `"[bam] {msg}"`. Zero runtime cost; improves `--verbose` output
readability.

**Implementation steps:**
1. Add `fn category(&self) -> &'static str` to `Error` via a match arm per
   variant.
2. Update `Display` / `Debug` formatting to include the category prefix.
3. No `From` conversions; no new enum variants.

## v0.7 — Architecture

### Config Composition via Clap Flatten

**Goal:** Reduce the 40-field flat `Config` struct into composed domain configurations,
improving cohesion and reducing `validate_and_init` line count.

**Prerequisite:** All existing tests pass; no open refactoring PRs.

**Implementation steps:**

1. Define three sub-structs using `#[derive(Parser)]` with clap's `#[command(flatten)]`:

```rust
   #[derive(Parser, Clone, Debug, Default)]
   pub struct ScoringConfig {
       #[arg(short, long, default_value = "4")]  pub mismatch_penalty: f64,
       #[arg(short, long, default_value = "6")]  pub gap_open: f64,
       #[arg(short = 'e', long, default_value = "1")] pub gap_extend: f64,
       #[arg(short = 'c', long, default_value = "5")] pub clipping_penalty: f64,
       #[arg(long, default_value = "illumina")] pub error_model: ErrorModel,
       #[arg(long)]                             pub bisulfite: bool,
       #[arg(long, default_value_t = u32::MAX)] pub ambiguous_threshold: u32,
       #[arg(short = 'J', long, default_value = "20")] pub chimeric_junction_bases: u32,
   }

   #[derive(Parser, Clone, Debug, Default)]
   pub struct ParallelConfig {
       #[arg(short = 't', long, default_value = "4")] pub threads: usize,
       #[arg(short = 'S', long, default_value = "1")] pub score_threads: usize,
   }

   #[derive(Parser, Clone, Debug, Default)]
   pub struct VariantConfig {
       #[arg(short, long, num_args = 0..4)] pub sample_variants: Vec<String>,
       #[arg(short, long, num_args = 0..4)] pub population_variants: Vec<String>,
       #[arg(long)]                         pub expand_indels: bool,
       #[arg(long, requires = "reference")] pub reference: Option<PathBuf>,
       // ... etc.
   }
```

2. In the root `Config`, replace individual fields with:
```rust
   #[command(flatten)]  pub scoring:  ScoringConfig,
   #[command(flatten)]  pub parallel: ParallelConfig,
   #[command(flatten)]  pub variants: VariantConfig,
```

3. Update `validate_and_init` to delegate to sub-struct validators:
   `self.scoring.validate()`, `self.parallel.validate()`, etc.

4. Update all call sites (grep for `config.mismatch_penalty`, etc.) to use
   `config.scoring.mismatch_penalty`. Expected ~150 mechanical substitutions.

5. Validate: `cargo test --all` must pass; `cargo clippy -- -D warnings`.

**Deferred because:** Step 4 touches every file that reads config fields.
The mechanical substitution is high-volume and carries merge-conflict risk.
Doing it in one focused PR is safer than mixing with feature work.

---

### AlignmentStream Trait Extraction

**Goal:** Move the `AlignmentStream` trait out of `src/aln_stream.rs` into
`src/stream/traits.rs` so that `aln_stream.rs` contains only the concrete
`AlnStream` struct, not the interface contract.

**Prerequisite:** Config composition complete (reduces cross-module coupling first).

**Implementation steps:**

1. Create `src/stream/mod.rs` and `src/stream/traits.rs`.
2. Move `pub(crate) trait AlignmentStream<R: SimpleRec>` and `FromBamRecord`
   into `src/stream/traits.rs`.
3. Add `pub(crate) mod stream;` to `src/lib.rs`.
4. Update all `use crate::aln_stream::AlignmentStream` imports.
5. In `src/aln_stream.rs`, `use crate::stream::traits::{AlignmentStream, FromBamRecord};`.
6. Run `cargo check` to enumerate remaining import failures; fix mechanically.

**Deferred because:** `AlignmentStream<R>` is a generic trait with `R: SimpleRec`
bounded receivers. The current test infrastructure (`MockStream` in
`src/aln_stream/tests.rs`) closely couples to the same file. Moving the trait
requires moving or re-exporting the test mock as well, expanding the scope.

---

### Workspace Crate Extraction

**Goal:** Extract `src/bam/`, `src/alignment/ops.rs`, and `src/variant/store.rs`
into separate workspace crates (`crates/xf-bam`, `crates/xf-align-math`,
`crates/xf-variant-store`) to enforce strict dependency boundaries and improve
incremental compilation times.

**Prerequisite:** Both Config composition and AlignmentStream extraction complete.

**Estimated impact:** ~30% reduction in full-rebuild time for the main crate.
Requires `cargo workspaces` tooling for version management.

**Deferred because:** This restructuring changes every internal import path and
requires deciding on stable inter-crate API surfaces before `v1.0` publication.
Profile first: run `cargo build --timings` on a real PDX dataset rebuild to
confirm the compile-time bottleneck is in these specific modules before investing.

## v0.8 — Runtime wiring completion

### noodles-util unified reader integration
`open_streams_unified` currently assumes hand-written per-extension dispatch
(SAM/BAM/CRAM) rather than `noodles_util::alignment::io::reader::Builder`.
Needs verification of exact builder API at the pinned noodles version;
swap in once confirmed. Fallback dispatch is functionally complete but
duplicates format-detection logic noodles-util would centralize.

### TabixScored — positive-score tabix-indexed BED reader
`open_tabix_scored` references a `TabixScored` type that does not exist yet.
Implement analogous to `TabixBed`/`TabixVcf`: tabix query → parse score
(column 5) and strand (column 6) → return `ScoredRegion` list for the
queried interval. Required for `collated --positive-regions`.

### Wire positive_regions into score_state_nw call sites
`HashLookup`/`CollatedMatcher` now carry a `positive: [Option<(Regions, ScoreFn)>; 2]`
field but the NW scoring calls (`nw_score_records`, `nw_score_fragment`) do not
yet pass it through to `score_state_nw`'s `positive_regions` parameter.
Three call-site edits; deferred to keep this change reviewable in isolation.

### chimeric_junction_bases as a CLI flag
Currently hardcoded to 20 in `ScoringArgs::to_penalty()`. Promote to a
`--chimeric-junction-bases` flag on `ScoringArgs` (shared across all
subcommands) once the default value has been empirically validated.

### CRAM index seeking for HashLookup
`open_streams_raw_bam` rejects CRAM input with a clear error. Implementing
CRAM support requires `.crai` index seeking in place of BGZF virtual offsets —
a distinct code path in `fetch_by_virtual_offset`. Significant scope; tracked
separately from this wiring pass.

---

## v1.0 — Production

- ○ **Publish to crates.io** — once variant rescue and the test suite are complete.
- ○ **Benchmark vs XenofilteR (R) and xenome** — reproducible benchmark using a public PDX
  dataset (e.g. SRP072062) with known ground truth.
- ○ **Conda / Bioconda packaging** — recipe for the `bioconda` channel.
- ○ **GitHub Actions CI** — build matrix (Linux x86_64, Linux aarch64, macOS);
  `cargo test`, `cargo clippy -- -D warnings`, `cargo fmt --check`.

### crossbeam work-stealing for parallel fragment scoring

**Problem:** Bounded `crossbeam_channel` in `parallel_io_loop` causes head-of-line
blocking when a fragment with N overlapping variants takes O(N log N) to score.
Cheap fragments behind it in the channel wait despite idle workers.

**Affected configuration:** `--score-threads > 1` AND `--sample-variants` or
`--population-variants` with many overlapping sites. Negligible impact otherwise
(NW cost is O(L) and uniform without variants).

**Design:**

```
IO thread
  └- injector.push(bundle)          // non-blocking; unbounded Injector
  └- in_flight.fetch_add(1)
  └- throttle: spin while in_flight > MAX_IN_FLIGHT
                                     // MAX_IN_FLIGHT = score_threads * 4

N workers (each owns a Worker<FragmentBundle> + knows all Stealers)
  loop {
    bundle = local.pop()
        or steal_from_siblings()
        or injector.steal()
        or park_until_woken();
    score_bundle(bundle);
    result_tx.send(scored);           // bounded channel back to IO thread
    in_flight.fetch_sub(1);
  }
```

Back-pressure: `in_flight` atomic counter + IO thread spin. Prevents injector
from growing unboundedly when scoring is slower than reading. `MAX_IN_FLIGHT =
score_threads × 4` is a starting point; needs tuning against real PDX datasets.

**Changes required:**
- Replace `(work_tx, work_rx): (Sender, Receiver)` with
  `crossbeam::deque::Injector<FragmentBundle>` + per-worker
  `Worker<FragmentBundle>` + `Vec<Stealer<FragmentBundle>>`.
- Worker thread startup: create `Worker::new_fifo()`, share `Arc<[Stealer]>`.
- IO thread wakeup: replace channel blocking with `Parker`/`Unparker` from
  `crossbeam_utils::sync`.
- Result collection: keep existing bounded `result_tx`/`result_rx` channel
  (already non-blocking for IO thread).

**Prerequisite:** Criterion benchmark showing measurable head-of-line blocking
on a dataset with variant rescue enabled (`--sample-variants`) and high variant
density (>1 variant per 10 bp). Without this evidence, the simpler bounded
channel is correct.

**Risk:** injector is not FIFO across stealing. Output order already
non-deterministic in parallel mode; this is acceptable.
```

---

## Known limitations (not on the roadmap)

- **Hash-map mode for non-name-sorted input** is explicitly out of scope for `namesorted`.
  Use `samtools sort -n` to pre-sort, or use `hashlookup`.
- **CRAM writing** is not planned in the near term; bgzf BAM is the standard output format.
- **`p_variant > 0.5` requirement** is a structural constraint of the rescue scoring formula,
  not a bug. Below 0.5 the delta is always non-positive.
