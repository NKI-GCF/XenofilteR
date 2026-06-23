# xenofilters — Algorithm Reference

## Shared Scoring Cascade

All three algorithms apply the same scoring cascade once a fragment's records can be evaluated. For `Namesorted`, this
occurs when records from all streams are aggregated. For `HashLookup` and `Collated`, an early decision can be made
using data from **only a single stream** if all gathered reads for that stream are entirely unmapped or are primary-only
perfect matches outside of ambiguous regions. In such scenarios, only the read-name key and the decision value are
stored, drastically reducing hot-path memory allocations by discarding the corresponding reads of the losing stream
immediately upon arrival.

### Tier 1 — Unmapped fast-path

`FragmentState::partial_cmp` checks whether all primaries in a fragment are unmapped (flags: `is_unmapped && (!is_segmented || is_mate_unmapped)`).

* One stream unmapped, other mapped → mapped wins immediately; no scoring. If all reads for stream 0 are unmapped, they
  are discarded; stream 1 records are written out if mapped.
* Both unmapped → tie (ambiguous or discard).

### Tier 2 — Perfect-match fast-path

`FragmentState::cmp_perfect` builds `MdCigFlags` for every record and checks:

```
is_perfect ≡ cigar.len() == 1  &&  md.bytes().all(is_ascii_digit)

```

**Supplementary** records (non-overlapping chimeric or structural biological segments) suppress the `is_perfect` flag
for their stream. **Secondary** records (overlapping alternative alignments) do not suppress the `is_perfect` flag; they
tag along and do not influence early sorting.

* One stream perfect, other not → perfect wins; no scoring. Loser stream is routed to `--discarded-output`.
* Both perfect → tie.
* Neither perfect → fall through to Tier 3.

BED/VCF region overlap (`Collated` and `HashLookup` only) forces a fall-through to Tier 3 even when Tier 2 would
otherwise resolve the fragment. The `Namesorted` algorithm explicitly rejects ambiguous BED/VCF files as arguments
during initialization.

### Tier 3 — Per-base log-likelihood scoring

Before running full dynamic programming, a **CIGAR/MD Subsumption** check is executed: if all records for all streams
are accounted for and the CIGAR/MD alignment profile of one stream is structurally a subset of (and inferior to) the
other, the superior stream wins, bypassing alignment scoring.

Otherwise, `Fragment::score` iterates `ScoreOpIter` (CIGAR + MD joint iterator) and accumulates per-base log-likelihood
penalties:

| Op type            | Penalty                                       |
|--------------------|-----------------------------------------------|
| Match              | `log_likelihood_match[q]`                     |
| Mismatch           | `log_likelihood_mismatch[q]`                  |
| Soft-clip          | `log_likelihood_mismatch[q]` per clipped base |
| Insertion/Deletion | affine: `gap_open + len × gap_extend`         |
| RefSkip            | zero                                          |

Quality scores index the penalty arrays (Phred-capped at `MAX_Q = 93`).

The total fragment score is calculated across all **primary and supplementary** reads. Supplementary reads incur an
structural configuration penalty computed as:


$$\text{Penalty} = \text{gap\_open} + (\text{clipped\_bases\_len} \times \text{gap\_extend})$$


where `clipped_bases_len` corresponds to the clipped bases in the primary mapping hosting the non-overlapping part of
the supplementary read.

### Tier 3 — Variant-aware rescue (optional)

When a VCF is supplied per stream, `Fragment::score` additionally runs `score_variants_in_window` for every scoring
window. Population, sample, and strain-diagnostic variants are evaluated here. Conversely, *ambiguous VCFs* can be
leveraged via fast lookups against CIGAR/MD strings to explicitly accelerate early assignment to the alternate strain.

For each overlapping variant:

1. **Weighted reference baseline** — `weighted_ref_score` computes the expected log-likelihood of the window assuming
   mixture prior `p_variant`:
$$score = \sum_j (1 - p) \times lm[q_j] + p \times lmm[q_j]$$


2. **Alt alignment** — `align_alt_to_read` runs Needleman–Wunsch DP of the alt allele against the read window, weighted
   by `p_variant`.
3. **Delta** stored in `Eval`: `alt_score - incurred`. Positive only when $p_{\text{variant}} > 0.5$.

Variants fully processed within a window are moved to a `finished` accumulator to prevent double-counting across
windows; they are merged back before the final scheduling step.

**Weighted Interval Scheduling** — `maximize_delta` runs an $O(n \log n)$ DP over all `Eval` entries with positive
delta, sorting by `end()` coordinate and using `partition_point` to find the latest non-overlapping predecessor. The
result is the maximum non-overlapping sum of variant rescue deltas.

The rescue delta is added to the fragment's log-likelihood score. A non-zero `scratch.last_variant_delta` causes the
decision tag to be written as `XR` (variant-rescued) instead of `XF`.

### Decision

```
delta = score(stream_0) - score(stream_1)

```

* `delta >  ambiguous_threshold` → stream 0 wins; records written to `--output[0]`, losers to `--discarded-output[1]`
* `delta < -ambiguous_threshold` → stream 1 wins
* `|delta| ≤ threshold` → ambiguous; both go to `--ambiguous-output`

With `--add-decision-tag`, winners receive `XF:C:<phred>` where $\text{phred} = \text{round}(10 \times |\delta| / \ln 10)$, capped at 255.

---

## Algorithm 1 — `Namesorted` (`LineByLine`)

**Input requirement:** all BAMs must be in identical query-name order. Supports scalability to $N > 2$ low-memory
alignment streams simultaneously.

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

The IO thread walks all streams in round-robin order. Records are appended to `AlnBuffer<FragmentState>` until a new
qname is seen; the new record is pushed back via `un_next`. When `best.len() > 1` the cascade runs.

### Sequential mode (`score_threads == 1`)

`process_sequential` — one thread; output order is deterministic.

### Parallel mode (`score_threads > 1`)

`process_parallel` — IO thread assembles `FragmentBundle`s and sends them on a bounded `crossbeam_channel`. $N$ worker
threads each own a `Scratch` (NW DP tables) and call `score_bundle`. Results drain back to the IO thread which owns all
writers. Output order is **nondeterministic** across fragments.

`Arc<dyn StoreTrait>` clones ($O(1)$ atomic increment) are distributed to workers at pipeline start; no locks during
scoring.

### Memory

$O(\text{reads-per-fragment} \times \text{streams})$. Essentially streaming.

---

## Algorithm 2 — `Hashlookup` (`HashLookup`)

**Input requirement:** position-sorted, indexed BAMs. No ordering guarantee between streams. Two-pass. Multi-stream
processing ($N > 2$) is bounded by strict heap memory limits due to cross-stream sorting skews.

### Pass 1 — Lightweight scan

`read_scoring_record` reads each BAM record and extracts a `ScoringRecord` (flags, ref_id, pos, cigar bytes, MD,
qualities, BGZF virtual offset). Sequences are **not** stored — qualities + MD + CIGAR suffice for NW scoring.

Records are inserted into a `NameTable` (`HashMap<Box<[u8]>, PendingFragment>`).

`PendingFragment::push` detects fragment completion:

* `expected_primaries` = 1 (single-end) or 2 (paired, detected from first record's `is_segmented` flag).
* Both streams must have `primary_count >= expected_primaries` before complete validation flags flip true.
* At classification time, each stream's `StreamBuf` is evaluated:
* All primaries perfect **and** no BED/VCF overlap → `StreamKind::Early` (virtual offsets only; `ScoringRecord`s dropped
  immediately).
* Otherwise → `StreamKind::Scoring` (records retained in a `Box<SmallVec>`).



Supplementary records never count as primaries; their virtual offsets are stored separately.

### Scoring

```
(Early, Early)     → ambiguous (no scoring)
(Early, Scoring)   → Early wins
(Scoring, Early)   → Early wins
(Scoring, Scoring) → full Tier 3 scoring on retained ScoringRecord fields

```

For `(Scoring, Scoring)`, `score_and_build` reconstructs temporary `RecordBuf` values from the stored fields, runs the
cascade, and produces a `ScoredFragment`.

### Pass 2 — Selective seek

`emit_scored` hands off virtual offsets to a dedicated **Seek/IO Thread** running `fetch_by_virtual_offset` $\to$
`bam.seek_vpos(vpos)` on the inner reader to isolate disk bounds from computation.

`StagedOutput` holds a `BTreeMap<u64, ScoredFragment>` keyed by driving-stream insertion sequence number. `flush` emits
only contiguous entries from `next_emit`, preserving driving-stream order.

### Memory

$O(\text{in-flight fragments})$. Worst case: all fragments of one stream before any of the other — the entire first
stream stays in the hash table.

---

## Algorithm 3 — `Collated` (`CollatedMatcher`)

**Input requirement:** each BAM must be *collated* (all records sharing a name are contiguous) but different streams may
present fragments in completely arbitrary orders. Highly scalable to $N > 2$ low-memory streams.

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
                score_pair (Dispatched to Worker Pool)

```

`CollatedReader::next_fragment` groups consecutive records sharing the same canonical name (after optional `/1`/`/2`
stripping) into one `FragmentState`. A one-record peek buffer handles look-ahead.

`handle_fragment` checks the opposite stream's waiting map. On a match, `score_pair` tasks are dispatched out to a
parallel thread worker pool; otherwise, the minimal footprint fragment configuration is buffered.

### Region overlap (Collated-specific)

`fragment_overlaps_region` queries `TabixBed` / `TabixVcf` via tabix random access for each mapped primary record. A hit
forces Tier 3 scoring even when `cmp_perfect` would have resolved the pair.

The BED/VCF files must be bgzf-compressed and tabix-indexed. HashLookup uses in-memory alternatives instead.

### Output order

Not guaranteed. Matched pairs are written concurrently as they clear thread evaluation; unmatched structural singletons
are flushed at stream exhaustion.

### Memory

$O(\text{name-order skew between streams})$.

---

## Algorithm Comparison

| Property                | Namesorted                                   | HashLookup                    | Collated                             |
|-------------------------|----------------------------------------------|-------------------------------|--------------------------------------|
| BAM sort order required | identical name order, both streams           | position-sorted               | collated (name-grouped) per stream   |
| Inter-stream order      | must match                                   | arbitrary                     | arbitrary                            |
| Output order            | deterministic (sequential) / non- (parallel) | driving-stream order          | unordered (parallelized)             |
| Parallelism             | IO + N scoring workers                       | Async Seek IO Thread Pipeline | IO + N thread pool workers           |
| Memory                  | O(fragment)                                  | O(in-flight fragments)        | O(name-order skew)                   |
| Region files            | not supported                                | in-memory BED/VCF             | tabix-indexed BED/VCF                |
| Pass count              | 1                                            | 2 (scan + seek)               | 1                                    |
| Variant rescue          | yes                                          | yes                           | yes                                  |
| BAM index needed        | no                                           | yes (BGZF virtual offsets)    | no (tabix for region files only)     |
