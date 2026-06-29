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
is_perfect ≡ supp_count == 0        (SA:Z tag absent — no chimeric segments pending)
           ∧ cigar.len() == 1       (single alignment operation, end-to-end)
           ∧ md.bytes().all(digit)  (zero mismatches against the reference)
```

The `SA:Z` check is required because any primary record with a non-empty `SA:Z`
tag will be followed by supplementary records carrying chimeric-junction penalties.
A fragment with pending supplementaries cannot be declared perfect even if its
primary CIGAR and MD are flawless.

The third condition is critical. The SA:Z: tag (stored in `MdCigFlags` as
`has_supplementary`, read from `Tag::OTHER_ALIGNMENTS` during construction) encodes
the alignments of any supplementary records for that fragment. Its presence guarantees
supplementary records will follow in the stream; supplementaries always add a structural
penalty to the total score, so the fragment can never be "perfect" for fast-path
purposes. Only primary records need to be checked: primaries always precede their
supplementaries in the BAM stream, so checking SA:Z: on primaries is sufficient and
correct.

Supplementary records themselves do **not** have their own SA tag checked; the
per-primary SA:Z: count (`b.iter().filter(|&&c| c == b';').count()`) can be stored
alongside `MdCigFlags` in a `SmallVec<[usize; READ_CT]>` for downstream use during
scoring.

Decision outcomes:
- One stream perfect (and SA:Z:-free), other not → perfect wins; no scoring.
- Both perfect → tie.
- Neither perfect → Tier 2.5.

BED/VCF region overlap (`Collated` and `HashLookup` only) forces fall-through to
Tier 3 even when Tier 2 would otherwise resolve the fragment.

### Tier 2.5 — Effective-match-count domination

Before running full NW DP, the pre-assessment computes the **effective match count**
for each stream from the `ReadProfile` already built during the single CIGAR+MD walk
(no additional parsing). An effective match is a read-coordinate position classified
as `ReadOp::Match` — a CIGAR M base confirmed by the MD tag as correctly aligned.
For a fragment that includes supplementary alignments, the structural penalty for each
supplementary reduces the stream's effective match advantage:

```
effective_matches(stream) =
    Σ ReadOp::Match positions over all primary records
  − Σ structural_penalty_units(supplementary records)
```

where the structural penalty unit for one supplementary is
`|gap_open| + non_clipped_bases × |gap_extend|`, expressed in the same per-base cost
space as a mismatch. This means a stream with supplementaries must overcome their
penalty to dominate.

Domination rules derived from `compare_fragment_profiles`:
- If all per-mate positions where A and B differ favour A, and A's deletion counts
  are not worse than B's, A wins without NW (EarlyDecision Greater).
- Symmetrically for B.
- If both streams have positions that favour them (PartialScoring), or if either
  stream has insertions of differing length, fall through to full NW.
- If both streams are identical position-by-position, deletion counts break the tie;
  if those are also equal, the result is Equal → ambiguous.

This formulation is strictly more powerful than a three-axis aggregate check
(mismatches ≤, clips ≤, indels ≤) because it operates at per-position granularity,
handles paired-end fragments across both mates, and correctly accounts for the case
where one stream has more mismatches but fewer clips while the other has fewer
mismatches but more clips — situations the three-axis aggregate leaves as incomparable
but where a per-position comparison often resolves cleanly.

For `HashLookup`, which operates on raw `ScoringRecord` bytes before any `RecordBuf`
construction, a lighter aggregate check (mismatch count from raw MD + clip/indel count
from raw CIGAR bytes) is used as a pre-filter. If it fires, the fragment is resolved
without NW; otherwise `score_records` runs and builds full profiles internally.
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

**Supplementary alignment scoring.** Supplementary records (chimeric read segments)
contribute to scoring in two additive components:

1. **Per-base NW score** — the supplementary's aligned bases are scored through the
   standard log-likelihood pipeline, using the sequence and quality scores from the
   primary record. The primary is guaranteed to be soft-clipped (not hard-clipped), so
   the full read sequence is always available.
2. **Structural penalty** — an additional cost of `gap_open + non_clipped_bases × gap_extend`
   is added to account for the chimeric junction. `non_clipped_bases` is the read
   length minus the count of soft- and hard-clipped bases in the supplementary's own
   CIGAR.

The two components combine so that even a supplementary with perfectly-aligning bases
always produces a lower total score than the same bases aligned without a chimeric
break, and always scores lower than a single-base indel with equal matched content.
This ensures supplementaries penalise their stream appropriately without being ignored.
Secondary alignments (overlapping alternatives, flag 0x100) tag along with their
stream's decision without contributing to scoring.

The total fragment score is the sum of per-base NW scores over all primary and
supplementary records, plus structural penalties over supplementary records.
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
 No --threads flag currently

### Chimeric fragment routing (optional, namesorted + paired-end)

When `--chimeric-pairs A:B` is supplied, the cascade is preceded by a
complementary-mapping check:

1. Look up `FragmentState` for stream A and stream B in the current `FragmentBuffer`.
2. Collect *mapped segment identifiers* (0x40 = read 1, 0x80 = read 2) for each
   stream, considering only primary, mapped, non-supplementary records.
3. **Chimeric event**: both sets non-empty AND disjoint (no mate maps well in both
   streams simultaneously).

On a chimeric event:
- Stream A records → stream A assigned output + `XC:Z:<label_B>` tag.
- Stream B records → stream B assigned output + `XC:Z:<label_A>` tag.
- Remaining streams (outside the pair) → filtered (discarded) output.
- `routing_counters[COUNTER_CHIMERIC_BASE + i]` incremented for each chimeric stream.

The normal tournament cascade is **skipped** for detected chimeric fragments.
Non-chimeric fragments proceed through the full Tier 1–3 cascade unchanged.

Condition 2 (disjointness) ensures that genuinely ambiguous reads — where a
single mate aligns with high score to both species — are not misclassified
as chimeric integration events.


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

`resolve()` applies the cascade in order: `partial_cmp` → `cmp_perfect` (including
SA:Z: check) → `pre_assess_alignments` (Tier 2.5, single CIGAR+MD walk) → `score_delta`
(full NW). Each tier short-circuits the ones below it.
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
               (is_perfect includes the SA:Z:-absent check)
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
separately in `supplementary_offsets` for pass-2 retrieval and scoring.
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
3. **Tier 2 perfect-match fast-path** — `cmp_perfect` (SA:Z:-aware); bypassed if
   region overlap.
4. **Tier 2.5 effective-match-count domination** — `pre_assess_alignments` builds
   `ReadProfile`s in a single CIGAR+MD walk; supplementary penalty units are
   subtracted from effective match counts before comparison.
5. **Tier 3 full NW** — `score_one` on both fragments.

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

| Property                | Namesorted                                     | HashLookup                       | Collated                              |
|-------------------------|------------------------------------------------|----------------------------------|---------------------------------------|
| BAM sort order required | identical name order, all streams              | position-sorted                  | collated (name-grouped) per stream    |
| Inter-stream order      | must match                                     | arbitrary                        | arbitrary                             |
| Output order            | deterministic (sequential) / non- (parallel)   | driving-stream order             | unordered                             |
| Parallelism             | IO + N scoring workers (crossbeam)             | single-threaded (seek-IO stub)   | IO + N thread pool (stub)             |
| Memory                  | O(fragment)                                    | O(in-flight fragments)           | O(name-order skew)                    |
| Region files            | not supported                                  | in-memory BED/VCF                | tabix-indexed BED/VCF                 |
| Pass count              | 1                                              | 2 (scan + seek)                  | 1                                     |
| Variant rescue          | yes                                            | yes                              | yes                                   |
| BAM index needed        | no                                             | yes (BGZF virtual offsets)       | no (tabix for region files only)      |
| EarlyKind fast-path     | Tier 1 + Tier 2 (SA:Z:-aware)                  | AllUnmapped / AllPerfect         | Tier 1 unmapped before tabix I/O      |
| Supplementary scoring   | per-base NW + structural penalty               | per-base NW + structural penalty | per-base NW + structural penalty      |
| Subsumption (Tier 2.5)  | effective match count (ReadProfile)            | aggregate raw-byte pre-check     | effective match count (ReadProfile)   |
| N > 2 streams           | architecturally ready                          | high-memory risk                 | architecturally ready                 |
