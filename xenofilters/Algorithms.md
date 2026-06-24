# xenofilters — Algorithm Reference

## Shared Scoring Cascade

All three algorithms apply the same cascade once a fragment's records can be evaluated.
`HashLookup` and `Collated` can fast-path single streams before the second stream
arrives; `Namesorted` always waits for records from every stream before deciding.

### Tier 1 — Unmapped fast-path

`FragmentState::partial_cmp` checks whether all primaries in a fragment are unmapped
(flags: `is_unmapped && (!is_segmented || is_mate_unmapped)`).

- One stream unmapped, other mapped → mapped stream wins immediately; no scoring.
- Both unmapped → tie (ambiguous or discard).

### Tier 2 — Perfect-match fast-path

`FragmentState::cmp_perfect` builds `MdCigFlags` for every record and checks:

```
is_perfect ≡ cigar.len() == 1  &&  md.bytes().all(is_ascii_digit)
```

**Supplementary** records (non-overlapping chimeric segments) suppress the `is_perfect`
flag for their stream. **Secondary** records (overlapping alternative alignments) do not.

- One stream perfect, other not → perfect wins; no scoring.
- Both perfect → tie.
- Neither perfect → Tier 2.5.

BED/VCF region overlap (`Collated` and `HashLookup` only) forces fall-through to Tier 3
even when Tier 2 would otherwise resolve the fragment.

### Tier 2.5 — CIGAR/MD structural subsumption

Before running full NW DP, `alignment_sig` computes aggregate mismatch count, soft-clip
count, and indel base count for each stream's `MdCigFlags`. `subsumes` then checks
whether one stream's signature is dominated on all three axes simultaneously:

```
a dominates b ≡ a.mismatches ≤ b.mismatches
              ∧ a.soft_clips  ≤ b.soft_clips
              ∧ a.indel_bases ≤ b.indel_bases
```

If one stream strictly dominates, it wins without any per-base scoring. If both
dominate equally (tie), the result is `Ordering::Equal` → ambiguous. If neither
dominates, Tier 3 runs.

Subsumption is only applied when both sides have the same number of primary segments;
mismatched split structures fall through to full NW unconditionally.

### Tier 3 — Per-base log-likelihood scoring

`Fragment::score` iterates `ScoreOpIter` (CIGAR + MD joint iterator) and accumulates
per-base log-likelihood penalties:

| Op type            | Penalty                                        |
|--------------------|------------------------------------------------|
| Match              | `log_likelihood_match[q]`                      |
| Mismatch           | `log_likelihood_mismatch[q]`                   |
| Soft-clip          | `log_likelihood_mismatch[q]` per clipped base  |
| Insertion/Deletion | affine: `gap_open + len × gap_extend`          |
| RefSkip            | zero                                           |

Quality scores index the penalty arrays (Phred-capped at `MAX_Q = 93`).

**Supplementary alignment structural penalty.** Supplementary records are excluded from
the per-base NW segment. Instead, each supplementary contributes:

```
penalty = gap_open + non_clipped_bases × gap_extend
```

where `non_clipped_bases = read_length − clip_count` (soft + hard clips from the
supplementary record's own CIGAR). This accounts for the chimeric junction without
double-penalising the primary mapping's soft clips.

The total fragment score is the sum of per-base NW scores over primary records plus
structural penalties over supplementary records.

### Tier 3 — Variant-aware rescue (optional)

When a VCF is supplied per stream, `Fragment::score` additionally runs
`score_variants_in_window` for every scoring window.

For each overlapping variant:

1. **Weighted reference baseline** — `weighted_ref_score` computes the expected
   log-likelihood of the window assuming mixture prior `p_variant`:

   ```
   score = Σ_j (1 - p) × lm[q_j] + p × lmm[q_j]
   ```

2. **Alt alignment** — `align_alt_to_read` runs Needleman–Wunsch DP of the alt allele
   against the read window, weighted by `p_variant`.

3. **Delta** stored in `Eval`: `alt_score − incurred`. Positive only when `p_variant > 0.5`.
   This is a structural invariant of the scoring formula; rescue cannot occur below this
   threshold.

Variants fully processed within a window are moved to a `finished` accumulator to
prevent double-counting across windows.

**Weighted Interval Scheduling** — `maximize_delta` runs an O(n log n) DP over all
`Eval` entries with positive delta, sorting by `end()` coordinate and using
`partition_point` to find the latest non-overlapping predecessor. The result is the
maximum non-overlapping sum of variant rescue deltas.

The rescue delta is added to the fragment score. A non-zero
`scratch.last_variant_delta` causes the decision tag to be written as `XR`
(variant-rescued) instead of `XF`.

### Decision

```
delta = score(stream_0) − score(stream_1)
```

- `delta >  ambiguous_threshold` → stream 0 wins
- `delta < -ambiguous_threshold` → stream 1 wins
- `|delta| ≤ threshold` → ambiguous

With `--add-decision-tag`:
- `XF:C:<phred>` when the margin came from per-base log-likelihood
- `XR:C:<phred>` when `scratch.last_variant_delta > 0` (variant rescue tipped the balance)

`phred = round(10 × |delta| / ln 10)`, capped at 255.

---

## Algorithm 1 — `Namesorted` (`LineByLine`)

**Input requirement:** all BAMs must be in identical query-name order.

### Data flow

```
stream[0]  stream[1]  …
    │           │
    └───────────┴──► handle_record_is_fragment_finished → AlnBuffer
                                    │
                         resolve() ─┤
                                    ├─ Ordered(ord) → handle_ordering
                                    └─ Scored(delta, s1_vd, s2_vd) → apply_delta
                                              │
                                         write_record
```

The IO thread walks streams in round-robin order. Records are appended to
`AlnBuffer<FragmentState>` until a new qname is seen; the new record is pushed back via
`un_next`. When `best.len() > 1` the cascade runs.

`resolve()` applies the cascade in order: `partial_cmp` → `cmp_perfect` →
subsumption → `score_delta`. Each tier short-circuits the ones below it.

### Sequential mode (`score_threads == 1`)

`process_sequential` — one thread; output order is deterministic.

### Parallel mode (`score_threads > 1`)

`process_parallel` — IO thread assembles `FragmentBundle`s and sends them on a bounded
`crossbeam_channel`. `N` worker threads each own a `Scratch` (NW DP tables) and call
`score_bundle`. Results drain back to the IO thread which owns all writers. Output order
is **nondeterministic** across fragments.

`Arc<dyn StoreTrait>` clones (O(1) atomic increment) are distributed to workers at
pipeline start; no locks during scoring.

### N-stream note

`namesorted` supports N > 2 streams via `SmallVec`. The cascade always compares
`best[0]` vs `best[last]`; a full N-way tournament is on the roadmap.

### Memory

O(reads-per-fragment × streams). Essentially streaming.

---

## Algorithm 2 — `Hashlookup` (`HashLookup`)

**Input requirement:** position-sorted, indexed BAMs. No ordering guarantee between
streams. Two-pass. Single-threaded only.

### EarlyKind fast-path

`StreamBuf::classify` evaluates each stream at fragment-completion time:

```
AllUnmapped  ← all primaries have is_unmapped flag set
AllPerfect   ← all primaries: is_perfect() == true  AND  no BED/VCF overlap
Scoring      ← otherwise (retain ScoringRecord bodies for NW)
```

`resolve_fragment` dispatches on the 2×2 matrix of `(driving_kind, lookup_kind)`:

| Driving         | Lookup          | Decision                             |
|-----------------|-----------------|--------------------------------------|
| AllPerfect      | AllPerfect      | Ambiguous tie                        |
| AllUnmapped     | AllUnmapped     | Ambiguous tie                        |
| AllPerfect      | AllUnmapped     | Driving wins (stream 0)              |
| AllUnmapped     | AllPerfect      | Lookup wins (stream 1)               |
| AllPerfect      | Scoring         | Driving wins (stream 0)              |
| AllUnmapped     | Scoring         | Lookup wins (stream 1)               |
| Scoring         | AllPerfect      | Lookup wins (stream 1)               |
| Scoring         | AllUnmapped     | Driving wins (stream 0)              |
| Scoring         | Scoring         | Full NW via `score_and_build`        |

### Pass 1 — Lightweight scan

`read_scoring_record` reads flags, ref_id, pos, CIGAR bytes, MD, qualities, and BGZF
virtual offset. No sequence. Records insert into `NameTable`
(`HashMap<Box<[u8]>, PendingFragment>`).

`PendingFragment::push` tracks `is_paired` from the first record's `is_segmented` flag.
`expected_primaries` returns 1 (single-end) or 2 (paired). Both streams must reach
`primary_count >= expected_primaries` before the fragment is complete and classifies.

Supplementary records never count as primaries; their virtual offsets are stored
separately in `supplementary_offsets`.

### Pass 2 — Selective seek

`StagedOutput` holds a `BTreeMap<u64, ScoredFragment>` keyed by driving-stream sequence
number. `flush` emits only contiguous entries from `next_emit`, preserving order.
`fetch_by_virtual_offset` calls `bam.seek_vpos(vpos)` on the inner BGZF reader.

Concurrency stub: see `src/filter_algorithm/hash_lookup/stage.rs` for a channel-based
seek-IO thread design that can overlap disk seeks with scoring work.

### Memory

O(in-flight fragments). Worst case: all fragments of stream 0 before any of stream 1.
N > 2 streams are flagged as high-memory risk.

---

## Algorithm 3 — `Collated` (`CollatedMatcher`)

**Input requirement:** each BAM must be collated (all records sharing a name are
contiguous) but the two streams may present fragments in arbitrary relative order.

### Data flow

```
CollatedReader[0]   CollatedReader[1]
      │                    │
next_fragment()      next_fragment()
      │                    │
      └──────► handle_fragment ◄──┘
                     │
          waiting_a / waiting_b  (HashMap buffers)
                     │
              score_pair
```

`CollatedReader::next_fragment` groups consecutive same-name records into one
`FragmentState` using a one-record peek buffer.

### Cascade in `score_pair`

1. **Tier 1 unmapped fast-path** — `partial_cmp` before any BED/VCF I/O (avoids
   tabix overhead for entirely unmapped fragments).
2. **BED/VCF region check** — `fragment_overlaps_region` queries `TabixBed` /
   `TabixVcf` for each mapped primary. A hit forces full NW even if `cmp_perfect`
   would resolve the pair.
3. **Tier 2 perfect-match fast-path** — `cmp_perfect`; bypassed if region overlap.
4. **Tier 3 full NW** — `score_one` on both fragments.

### Region overlap (Collated-specific)

BED/VCF files must be bgzf-compressed and tabix-indexed. `HashLookup` uses in-memory
`AmbiguousRegions` / `DiagnosticVariants` loaded at startup instead.

### Output order

Not guaranteed. Matched pairs write concurrently as they clear evaluation; unmatched
structural singletons flush at stream exhaustion.

Concurrency stub: see `src/filter_algorithm/collated.rs` for a thread-pool design
where `score_pair` can be dispatched to parallel workers.

### Memory

O(name-order skew between streams). Scales cleanly to N > 2 via N waiting maps.

---

## Algorithm Comparison

| Property                | Namesorted                                    | HashLookup                    | Collated                              |
|-------------------------|-----------------------------------------------|-------------------------------|---------------------------------------|
| BAM sort order required | identical name order, all streams             | position-sorted               | collated (name-grouped) per stream    |
| Inter-stream order      | must match                                    | arbitrary                     | arbitrary                             |
| Output order            | deterministic (sequential) / non- (parallel)  | driving-stream order          | unordered                             |
| Parallelism             | IO + N scoring workers (crossbeam)            | single-threaded (seek-IO stub)| IO + N thread pool (stub)             |
| Memory                  | O(fragment)                                   | O(in-flight fragments)        | O(name-order skew)                    |
| Region files            | not supported                                 | in-memory BED/VCF             | tabix-indexed BED/VCF                 |
| Pass count              | 1                                             | 2 (scan + seek)               | 1                                     |
| Variant rescue          | yes                                           | yes                           | yes                                   |
| BAM index needed        | no                                            | yes (BGZF virtual offsets)    | no (tabix for region files only)      |
| EarlyKind fast-path     | Tier 1 + Tier 2 (no explicit EarlyKind)       | AllUnmapped / AllPerfect      | Tier 1 unmapped before tabix I/O      |
| Supplementary penalty   | yes (sequential + parallel)                   | yes (score_records)           | yes (score_one)                       |
| Subsumption (Tier 2.5)  | yes (LineByLine::resolve)                     | no (uses EarlyKind instead)   | no (uses TabixBed/Vcf instead)        |
| N > 2 streams           | architecturally ready                         | high-memory risk              | architecturally ready                 |
