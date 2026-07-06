# xenofilter — Rust implementation

Fast, streaming read classifier for xenograft / PDX data.

Compares multiple aligned BAM streams (host vs graft vs …) per fragment and decides:

- alignment is better than in other species → `-o` / stdout
- alignment is equally good / too close → `--ambiguous-output`
- alignment is worse than in other species → `--discarded-output`

Use cases: PDX host/graft separation, aligner comparison, viral integration, dual-genome
cell-line work, multi-strain within-species disambiguation.

---

## Architecture

Three backend matchers are fully wired and dispatched from `main.rs`:

**`LineByLine` (`--matching-algorithm namesorted`).**
Streaming merge of name-sorted BAM streams. Reads advance in lockstep; a fragment is
complete when a new query name arrives. Supports `N > 2` streams via `SmallVec`. Parallel
scoring via `crossbeam_channel` bounded channels: IO thread reads and writes; `N` worker
threads each own a `Scratch` (NW DP tables) and call `score_bundle`. `Arc<dyn StoreTrait>`
clones give workers zero-cost access to variant stores.

**`HashLookup` (`--matching-algorithm hashlookup`).**
Two-pass design for position-sorted BAMs. Pass 1 reads lightweight `ScoringRecord`s
(no sequences) and inserts into `NameTable` (`HashMap<Box<[u8]>, PendingFragment>`).
At fragment completion each stream is classified as `StreamKind::Early` (with
`EarlyKind::AllUnmapped` or `EarlyKind::AllPerfect`) or `StreamKind::Scoring`.
Early streams drop their `ScoringRecord`s immediately. The `resolve_fragment` 2×2
matrix handles all `EarlyKind` combinations without scoring. Full NW runs only for
`(Scoring, Scoring)` pairs. Pass 2 seeks via BGZF virtual offsets using
`fetch_by_virtual_offset`. `StagedOutput` (a `BTreeMap<u64, ScoredFragment>`) preserves
driving-stream order.

**`CollatedMatcher` (`--matching-algorithm collated`).**
Each stream is individually collated. `CollatedReader` groups consecutive same-name
records. `waiting_a` / `waiting_b` `HashMap`s buffer unmatched fragments. Tier-1
unmapped fast-path runs before BED/VCF I/O. BED and VCF files are queried via tabix
random access (`TabixBed`, `TabixVcf`).

**Chimeric fragment detection** (namesorted only, paired-end).
When `--chimeric-pairs A:B` is configured, fragments where mates split
across two streams (one mate mapped in stream A, complementary mate in stream B)
are classified as chimeric.  Both streams' records are written to their assigned
outputs with `XC:Z:<other_stream_label>`.  Useful for viral integration studies
(e.g. HPV integration into human: `--chimeric-pairs 0:1`).

Three-stream example: human (0) + HPV (1) + mouse (2) xenograft.
`--chimeric-pairs 0:1` lets human and HPV form chimeric pairs while
the mouse stream competes normally in the tournament.

---

## Scoring cascade

All three backends share the same cascade:

| Tier | Condition | Cost |
|------|-----------|------|
| 1 | Unmapped vs mapped | O(1) flag check |
| 2 | Perfect vs imperfect (CIGAR len=1, MD all-digits) | O(records) |
| 2.5 | CIGAR/MD structural subsumption | O(records × ops) |
| 3 | Per-base log-likelihood + affine gap | O(read length) |
| 3+ | Variant-aware NW rescue + WIS scheduling | O(read × variants) |

Supplementary alignments contribute a structural penalty `gap_open + non_clipped_bases × gap_extend`
instead of per-base scoring.

---

## Decision tags

```
XF:C:<phred>   Score-based confidence (phred-scaled log-likelihood ratio)
XR:C:<phred>   Same, but the margin came from variant rescue (p_variant > 0.5)
```

---

## Features compared to original XenofilteR (R)

- Three matching algorithms (namesorted / hashlookup / collated)
- Parallel fragment scoring (namesorted + crossbeam workers)
- Supplementary alignment structural penalties
- CIGAR/MD subsumption pre-check (bypasses NW when one stream dominates)
- Optional variant-aware rescue (SNV + indel) via NW DP + weighted interval scheduling
- `EarlyKind`-aware fast-path matrix in HashLookup (AllUnmapped vs AllPerfect)
- Merged output mode (`--merged-output`) with `_xenofilt` / `_xenoambig` RG suffixes
- Optional `--ambiguous-regions` BED and `--diagnostic-variants` VCF forcing full scoring
- Phred-scaled decision tag `XF:C` / `XR:C` (variant-rescued)
- @PG header line
- Rust performance, safety, and zero-dependency binary

---

## Quick build

```bash
cargo build --release
./target/release/xenofilters host.bam graft.bam \
  -o best.bam -f filtered.bam -a ambiguous.bam \
  --add-decision-tag
```

Requires Rust ≥ 1.88 (edition 2024, `let`-chain syntax).

## Verification
cargo check                          # binary builds
cargo check --lib                    # library builds
cargo test --lib                     # lib unit tests
cargo bench --no-run --features bench-internals   # benches compile
cargo test --test integration        # integration tests compile
