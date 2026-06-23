# xenofilters — Algorithm Reference

## Shared Scoring Cascade

All three algorithms apply the same three-tier cascade once a fragment's
records from all streams have been assembled.

### Tier 1 — Unmapped fast-path
`FragmentState::partial_cmp` checks whether all primaries in a fragment are
unmapped (flags: `is_unmapped && (!is_segmented || is_mate_unmapped)`).

- One stream unmapped, other mapped → mapped wins immediately; no scoring.
- Both unmapped → tie (ambiguous or discard).

### Tier 2 — Perfect-match fast-path
`FragmentState::cmp_perfect` builds `MdCigFlags` for every record and checks:

```
is_perfect ≡ cigar.len() == 1  &&  md.bytes().all(is_ascii_digit)
```

Secondary records suppress the `is_perfect` flag for their stream.

- One stream perfect, other not → perfect wins; no scoring.
- Both perfect → tie.
- Neither perfect → fall through to Tier 3.

BED/VCF region overlap (Collated and HashLookup only) can force Tier 3 even
when Tier 2 would resolve the fragment.

### Tier 3 — Per-base log-likelihood scoring
`Fragment::score` iterates `ScoreOpIter` (CIGAR + MD joint iterator) and
accumulates per-base log-likelihood penalties:

| Op type  | Penalty                          |
|----------|----------------------------------|
| Match    | `log_likelihood_match[q]`        |
| Mismatch | `log_likelihood_mismatch[q]`     |
| Soft-clip | `log_likelihood_mismatch[q]` per clipped base |
| Insertion/Deletion | affine: `gap_open + len × gap_extend` |
| RefSkip  | zero                             |

Quality scores index the penalty arrays (Phred-capped at `MAX_Q = 93`).

The fragment score is the sum over all primary reads in the fragment.

### Tier 3 — Variant-aware rescue (optional)
When a VCF is supplied per stream, `Fragment::score` additionally runs
`score_variants_in_window` for every scoring window.

For each overlapping variant:

1. **Weighted reference baseline** — `weighted_ref_score` computes the expected
   log-likelihood of the window assuming mixture prior `p_variant`:

   ```
   score = Σ_j (1 - p) × lm[q_j]  +  p × lmm[q_j]
   ```

2. **Alt alignment** — `align_alt_to_read` runs Needleman–Wunsch DP of the
   alt allele against the read window, weighted by `p_variant`.

3. **Delta** stored in `Eval`: `alt_score - incurred`.  Positive only when
   `p_variant > 0.5` (known structural property of the formula).

Variants fully processed within a window are moved to a `finished` accumulator
to prevent double-counting across windows; they are merged back before the
final scheduling step.

**Weighted Interval Scheduling** — `maximize_delta` runs a O(n log n) DP over
all `Eval` entries with positive delta, sorting by `end()` coordinate and using
`partition_point` to find the latest non-overlapping predecessor.  The result is
the maximum non-overlapping sum of variant rescue deltas.

The rescue delta is added to the fragment's log-likelihood score.  A non-zero
`scratch.last_variant_delta` causes the decision tag to be written as
`XR` (variant-rescued) instead of `XF`.

### Decision
```
delta = score(stream_0) - score(stream_1)
```

- `delta >  ambiguous_threshold` → stream 0 wins; records written to `--output[0]`, losers to `--filtered-output[1]`
- `delta < -ambiguous_threshold` → stream 1 wins
- `|delta| ≤ threshold` → ambiguous; both go to `--ambiguous-output`

With `--add-decision-tag`, winners receive `XF:C:<phred>` where
`phred = round(10 × |delta| / ln10)`, capped at 255.

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
                                    └─ Scored(delta) → apply_delta
                                              │
                                         write_record
```

The IO thread walks all streams in round-robin order.  Records are appended to
`AlnBuffer<FragmentState>` until a new qname is seen; the new record is pushed
back via `un_next`.  When `best.len() > 1` the cascade runs.

### Sequential mode (`score_threads == 1`)
`process_sequential` — one thread; output order is deterministic.

### Parallel mode (`score_threads > 1`)
`process_parallel` — IO thread assembles `FragmentBundle`s and sends them on a
bounded `crossbeam_channel`.  N worker threads each own a `Scratch` (NW DP
tables) and call `score_bundle`.  Results drain back to the IO thread which owns
all writers.  Output order is **nondeterministic** across fragments.

`Arc<dyn StoreTrait>` clones (O(1) atomic increment) are distributed to workers
at pipeline start; no locks during scoring.

### Memory
O(reads-per-fragment × streams).  Essentially streaming.

---

## Algorithm 2 — `Hashlookup` (`HashLookup`)

**Input requirement:** position-sorted, indexed BAMs.  No ordering guarantee
between streams.  Two-pass.

### Pass 1 — Lightweight scan

`read_scoring_record` reads each BAM record and extracts a `ScoringRecord`
(flags, ref\_id, pos, cigar bytes, MD, qualities, BGZF virtual offset).
Sequences are **not** stored — qualities + MD + CIGAR suffice for NW scoring.

Records are inserted into a `NameTable` (`HashMap<Box<[u8]>, PendingFragment>`).

`PendingFragment::push` detects fragment completion:

- `expected_primaries` = 1 (single-end) or 2 (paired, detected from first record's `is_segmented` flag).
- Both streams must have `primary_count >= expected_primaries` before `is_complete_inner` returns true.
- At classification time, each stream's `StreamBuf` is evaluated:
  - All primaries perfect **and** no BED/VCF overlap → `StreamKind::Early` (virtual offsets only; `ScoringRecord`s dropped).
  - Otherwise → `StreamKind::Scoring` (records retained in a `Box<SmallVec>`).

Supplementary records never count as primaries; their virtual offsets are stored
separately.

### Scoring

```
(Early, Early)     → ambiguous (no scoring)
(Early, Scoring)   → Early wins
(Scoring, Early)   → Early wins
(Scoring, Scoring) → full Tier 3 scoring on retained ScoringRecord fields
```

For `(Scoring, Scoring)`, `score_and_build` reconstructs temporary `RecordBuf`
values from the stored fields, runs the cascade, and produces a `ScoredFragment`.

### Pass 2 — Selective seek

`emit_scored` calls `fetch_by_virtual_offset` → `bam.seek_vpos(vpos)` on the
inner `bgzf::io::Reader` to seek directly to each record for output.

`StagedOutput` holds a `BTreeMap<u64, ScoredFragment>` keyed by driving-stream
insertion sequence number.  `flush` emits only contiguous entries from
`next_emit`, preserving driving-stream order.

### Memory
O(in-flight fragments).  Worst case: all fragments of one stream before any of
the other — the entire first stream stays in the hash table.

---

## Algorithm 3 — `Collated` (`CollatedMatcher`)

**Input requirement:** each BAM must be *collated* (all records sharing a name
are contiguous) but the two streams may present fragments in different orders.

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

`CollatedReader::next_fragment` groups consecutive records sharing the same
canonical name (after optional `/1`/`/2` stripping) into one `FragmentState`.
A one-record peek buffer handles the look-ahead.

`handle_fragment` checks the opposite stream's waiting map.  On a match,
`score_pair` runs; otherwise the fragment is buffered.

### Region overlap (Collated-specific)
`fragment_overlaps_region` queries `TabixBed` / `TabixVcf` via tabix random
access for each mapped primary record.  A hit forces Tier 3 scoring even when
`cmp_perfect` would have resolved the pair.

The BED/VCF files must be bgzf-compressed and tabix-indexed.  HashLookup uses
in-memory `AmbiguousRegions`/`DiagnosticVariants` instead.

### Output order
Not guaranteed.  Matched pairs are emitted immediately; unmatched fragments
(one stream only) are emitted as winners at stream exhaustion.

### Memory
O(name-order skew between the two streams).

---

## Algorithm Comparison

| Property | Namesorted | HashLookup | Collated |
|---|---|---|---|
| BAM sort order required | identical name order, both streams | position-sorted | collated (name-grouped) per stream |
| Inter-stream order | must match | arbitrary | arbitrary |
| Output order | deterministic (sequential) / nondeterministic (parallel) | driving-stream order | unordered |
| Parallelism | IO + N scoring workers | single-threaded | single-threaded |
| Memory | O(fragment) | O(in-flight fragments) | O(name-order skew) |
| Region files | not supported | in-memory BED/VCF | tabix-indexed BED/VCF |
| Pass count | 1 | 2 (scan + seek) | 1 |
| Variant rescue | yes | yes | yes |
| BAM index needed | no | yes (BGZF virtual offsets) | no (tabix for region files only) |
