# Architectural Trade-off: Indel Equivalence for Un-realigned Xenograft BAMs

## Problem Statement

Standard aligners produce CIGAR/MD strings that represent the same
biological insertion or deletion at different positions inside homopolymers
and tandem repeats, depending on the aligner's heuristics and whether
`samtools calmd` / `GATK IndelRealigner` has been applied.  When a VCF
variant records an indel at position P and the read expresses the same
indel at position P+k, the existing per-base NW scoring window misses the
match: the window is opened at the VCF position, and the read bases at that
window represent a different alignment.

xenofilters' existing NW path (`align_alt_to_read`) correctly scores the
alt allele *given the right window*.  The problem is entirely in the
**lookup / window-selection step**, not in the scoring arithmetic.

---

## Approach A — On-the-Fly Local Realignment (Smith-Waterman / banded NW)

### Mechanism

For every read that overlaps a variant locus, reconstruct a small haplotype
sequence (ref with alt spliced in) and perform SW or banded NW alignment of
the read against it.  Feed the resulting alignment into the existing scoring
framework.

### Pros

- Biologically exact: handles any misrepresentation including MNPs,
  complex alleles, and novel allele combinations.
- No reference FASTA access required at startup.
- Clean conceptual extension of the existing NW framework.

### Cons

| Concern | Detail |
|---------|--------|
| **Hot-path cost** | O(R × H) DP per read per overlapping variant, where R ≈ read length (150–20 000 bp) and H ≈ haplotype window (100–500 bp). For ONT/HiFi reads with hundreds of overlapping variants this is catastrophic. |
| **Redundant work** | The same computation is repeated for every read at a locus; no caching across reads. |
| **Integration complexity** | `align_alt_to_read` would need extension to handle the full haplotype graph; the existing `Scratch` buffers would need resizing per call. |
| **Cache pressure** | Constructing haplotypes per read thrashes the data cache. |
| **Scoring consistency** | A second DP inside the first breaks the clean per-base log-likelihood model; delta interpretation becomes ambiguous. |

### Verdict: **rejected**

The per-read DP cost is prohibitive at scale.  The existing NW is
deliberately *variant-scoped* (short window, O(read_len × variant_len));
a haplotype-level realignment is an order of magnitude more expensive and
is not needed once the lookup step is fixed.

---

## Approach B — Indel Equivalence Enumerator (recommended)

### Mechanism

At tool startup, for each indel in the input VCFs:

1. Fetch a reference context window (≤ 200 bp) via indexed FASTA.
2. Left-normalize the (POS, REF, ALT) tuple.
3. Enumerate every right-shift through the local repeat region that
   produces a mathematically equivalent representation.
4. Insert each equivalent tuple as an independent `Variant` entry in the
   `Store<V>`, carrying the same scoring parameters (p_variant / GQ) as
   the canonical form.

During BAM scanning, the existing `store.overlapping_multi()` query finds
whichever representation the aligner chose, and `align_alt_to_read` scores
it correctly because the stored alleles match the window the read presents.

### Pros

| Concern | Detail |
|---------|--------|
| **Hot-path cost** | O(1) lookup per locus — equivalent to current SNP handling. |
| **Pre-computation** | Done once at startup; amortized over millions of reads. |
| **Integration** | Transparent to all downstream code; `Store`, `Eval`, `wis_max_rescue_delta` unchanged. |
| **Memory** | Bounded: a k-bp repeat generates at most k equivalent entries. Typical indels in coding regions produce 1–10 equivalents. |
| **Cache locality** | Sorted Vec of equivalent positions allows binary-search range queries with good locality. |
| **Correctness** | WIS naturally deduplicates overlapping equivalent entries (two entries for the same biological event cannot both be non-overlapping). |

### Cons

| Concern | Detail |
|---------|--------|
| **Reference FASTA required** | `--reference` must be supplied when `--sample-variants` or `--population-variants` contains indels.  Already required for CRAM; this extends the requirement to any indel-containing VCF. |
| **Complex alleles** | MNPs and multi-base complex alleles are not enumerated (they cannot slide); stored as-is with a warning. |
| **Startup latency** | FASTA fetches per variant add ~1–10 ms total for typical PDX VCFs (< 10 k indels); negligible. |

### Verdict: **adopted**

The bottleneck is the per-read scoring loop, not the per-variant startup
cost.  Approach B leaves the hot path unchanged and fixes the lookup step
precisely where the problem lies.

---

## Integration with Existing Architecture

```
VCF record
    │
    ▼  IndelEquivalenceExpander::expand()
    │  (one reference FASTA fetch per indel; ≤ MAX_SHIFT=100 new entries)
    │
    ├─► canonical Variant entry  ─────────────┐
    ├─► shifted Variant (pos+1, ...)  ────────┤
    ├─► shifted Variant (pos+2, ...)  ────────┤  Store<V>::insert_all()
    └─► ...                           ────────┘
                                              │
                                              ▼
                              store.overlapping_multi(tid, start, end)
                                              │
                                              ▼
                                 align_alt_to_read()   ← unchanged
                                              │
                                              ▼
                              wis_max_rescue_delta()   ← unchanged
```

The `Store<V>` receives expanded entries via a new
`insert_expanded(variants)` method.  Coordinate handling: VCF is 1-based;
all internal positions stored 0-based (`vcf_pos - 1`).
