# xenofilters

Fast, streaming read classifier for xenograft / PDX sequencing data — written in Rust.

---

## What it does

Patient-derived xenograft (PDX) models graft human tumour tissue into mouse hosts.
Sequencing a PDX sample therefore captures reads from **both** species.
`xenofilters` takes two BAM files — one aligned to the **graft** reference
(e.g. human GRCh38) and one to the **host** (e.g. mouse GRCm39) — and assigns each
read fragment to the species it aligns to best.

The same approach applies wherever you need to separate reads from two co-sequenced
genomes: viral integration studies, transposon experiments, dual-genome cell-line work,
or comparison of two different aligners on the same sample.

---

## Algorithm

For every fragment (one or two reads sharing the same name) `xenofilters` applies a
cascading decision pipeline:

**Tier 1 — Unmapped fast-path.** If all primaries in one stream are unmapped (flags:
`UNMAP` and either `!PAIRED` or `MUNMAP`), the mapped stream wins immediately without
scoring. If both are unmapped the result is a tie.

**Tier 2 — Perfect-match fast-path.** If one stream has every primary alignment with a
single-op CIGAR and an all-digit MD tag (no mismatches), and the other does not, the
perfect alignment wins. Supplementary alignments suppress the perfect flag for their
stream. Secondary alignments do not.

**Tier 2.5 — CIGAR/MD structural subsumption.** Before full DP scoring, a quick
comparison of each stream's aggregate mismatch, soft-clip, and indel counts is performed.
If one stream's alignment profile is a strict subset of the other's across all axes
(mismatches ≤, clips ≤, indels ≤), the structurally superior stream wins without
Needleman–Wunsch scoring.

**Tier 3 — Per-base log-likelihood scoring.** A per-base log-likelihood score is computed
from CIGAR + MD tag + base quality scores. Mismatches and soft-clips are penalised;
indels use affine gap costs. Supplementary alignments contribute a structural configuration
penalty (`gap_open + non_clipped_bases × gap_extend`) rather than per-base scoring.
The decision score is the sum over all primary reads in the fragment.

**Tier 3 — Variant-aware rescue** (optional, `-s` / `-p`). If a VCF/BCF is supplied,
known variants are scored via Needleman–Wunsch DP of the alt allele against the read.
A positive score delta (alt allele better supported than reference) is added to the
fragment score. Non-overlapping variants are combined with weighted interval scheduling
(`maximize_delta`). `p_variant > 0.5` is required for a rescue delta to be positive.

**Decision.** The stream with the higher score wins. If the Phred-scaled score difference
is below `--ambiguous-threshold` (default 0) the fragment is written to
`--ambiguous-output` instead.

The optional `XF:C:<phred>` aux tag encodes the Phred-scaled confidence. When the margin
came from variant rescue, `XR:C:<phred>` is written instead.

---

## Three matching algorithms

`xenofilters` supports three backend algorithms selectable via `--matching-algorithm`:

**`namesorted` (default).** Both BAMs must be in identical query-name order. Streaming,
lowest memory, fastest. Supports parallel scoring via `--score-threads`. Ambiguous BED/VCF
region files are not supported in this mode.

**`hashlookup`.** Accepts position-sorted BAMs with no ordering guarantee between streams.
Two-pass design: pass 1 reads lightweight `ScoringRecord`s and classifies streams as
`EarlyKind::AllUnmapped`, `EarlyKind::AllPerfect`, or `Scoring`. Pass 2 seeks via BGZF
virtual offsets to fetch full records for output. In-memory BED and diagnostic VCF files
force fragments through full scoring. Single-threaded only. Memory proportional to
in-flight fragments.

**`collated`.** Each BAM must be internally collated (all records for a name contiguous)
but the two streams may present fragments in any relative order. A per-stream
`HashMap` buffers unmatched fragments. Tabix-indexed BED and VCF files are queried
per fragment. Output order is not guaranteed. Memory proportional to name-order skew.

---

## Installation

### From source

```bash
git clone https://github.com/NKI-GCF/XenofilteR
cd XenofilteR
cargo build --release
# Binary is at ./target/release/xenofilters
```

**Requirements:** Rust 1.88+ (edition 2024, requires `let`-chain syntax).
No external C libraries required.

---

## Quick start

```bash
# Name-sorted mode (default) — align and sort by name first
bwa mem -M -t 8 human_ref.fa reads_R1.fq.gz reads_R2.fq.gz \
  | samtools sort -n -@ 4 -o human.bam

bwa mem -M -t 8 mouse_ref.fa reads_R1.fq.gz reads_R2.fq.gz \
  | samtools sort -n -@ 4 -o mouse.bam

xenofilters human.bam mouse.bam \
  --output human_best.bam \
  --output mouse_best.bam \
  --discarded-output human_discarded.bam \
  --discarded-output mouse_discarded.bam \
  --ambiguous-output ambiguous.bam \
  --threads 4 \
  --score-threads 4 \
  --verbose

# Position-sorted mode — no name-sort required
xenofilters human_coord.bam mouse_coord.bam \
  --matching-algorithm hashlookup \
  --ambiguous-regions human_mask.bed mouse_mask.bed \
  --diagnostic-variants human_diag.bcf mouse_diag.bcf \
  --output human_best.bam --output mouse_best.bam
```

---

## CLI reference

```
Usage: xenofilters [OPTIONS] <ALIGNMENT>...

Arguments:
  <ALIGNMENT>...
      Input BAMs. Two required except in --single-alignment-mode.

Options:
  -o, --output <PATH>...
      Output file for winning reads (one per alignment stream).

  -f, --discarded-output <PATH>...
      Output file for reads that lost. Default: discard.

  -a, --ambiguous-output <PATH>...
      Output file for reads with equally good alignments. Default: discard.

      --merged-output <PATH>
      Write all reads (winners, discarded, ambiguous) to a single BAM.
      Read groups are suffixed with _xenofilt or _xenoambig for non-winners.
      Mutually exclusive with --output, --discarded-output, --ambiguous-output.

  -O, --stdout-format <FORMAT>
      Output format for stdout [sam|bam|cram]. [default: sam]

  -t, --threads <N>
      bgzf (de)compression worker threads. [default: 4]

  -S, --score-threads <N>
      Parallel scoring worker threads (namesorted only).
      0 = use all available logical CPUs (up to 16). [default: 1]

      --matching-algorithm <ALG>
      namesorted | hashlookup | collated [default: namesorted]

  -v, --verbose
      Increase verbosity (-v = INFO, -vv = DEBUG). Respects RUST_LOG.

  -U, --discard-unmapped
      Discard fragments unmapped in every stream, even from discarded output.

  -m, --mismatch-penalty <F>     [default: 4]
  -g, --gap-open <F>             [default: 6]
  -e, --gap-extend <F>           [default: 1]
  -c, --clipping-penalty <F>     [default: 5]

  -a, --ambiguous-threshold <PHRED>
      Minimum score difference (Phred) to call a winner. [default: 0]

  -s, --sample-variants <[IDX:]FILE>...
      Sample-specific VCF/BCF (FORMAT/GT + FORMAT/GQ).

  -p, --population-variants <[IDX:]FILE>...
      Population VCF/BCF (INFO/AF).

  -A, --add-decision-tag
      Write XF:C or XR:C aux tag to winning records.

  -P, --no-program-line
      Suppress @PG header line.

  -R, --strip-read-suffix <MODE>
      Handle /1 /2 suffixes [auto|true|false|variable]. [default: auto]

      --ambiguous-regions <FILE>...
      BED file(s) of regions that force full NW scoring (one per stream).
      hashlookup: plain BED, loaded into memory.
      collated: bgzf-compressed + tabix-indexed (.bed.gz + .tbi).

      --diagnostic-variants <FILE>...
      VCF/BCF of species-diagnostic positions (one per stream).
      Same format rules as --ambiguous-regions.

      --single-alignment-mode
      Allow a single input BAM (requires two variant profiles, one per strain).

      --matching-algorithm <ALG>
      namesorted | hashlookup | collated [default: namesorted]

  -h, --help     Print help
  -V, --version  Print version
```

---

## Output files

| Stream | Reads written |
|--------|--------------|
| `--output` | Reads that aligned **better** to this reference |
| `--discarded-output` | Reads that aligned **worse** (won by another stream) |
| `--ambiguous-output` | Reads where no stream was clearly better |

`--merged-output` combines all three into one file using `RG:Z` tag suffixes.

---

## Variant-aware rescue

```bash
xenofilters strain_a.bam strain_b.bam \
  -s 0:strain_a_variants.vcf \
  -s 1:strain_b_variants.vcf \
  -o best_a.bam -o best_b.bam \
  --add-decision-tag
```

VCF requirements: BCF preferred for performance. Must be coordinate-sorted.
Sample VCF: `GT` and `GQ` FORMAT fields. Population VCF: `AF` INFO field.
`p_variant > 0.5` is required for rescue to produce a positive delta.

---

## Decision tags

When `--add-decision-tag` is set:

| Tag | Meaning |
|-----|---------|
| `XF:C:<phred>` | Score-based confidence; phred = round(10 × \|delta\| / ln 10), capped at 255 |
| `XR:C:<phred>` | Same, but the winning margin came from variant rescue |

---

## Environment variables

| Variable | Effect |
|----------|--------|
| `RUST_LOG` | Override log level (`xenofilters=debug`) |

---

## Performance notes

- **Memory:** `namesorted` O(reads-per-fragment × streams). `hashlookup` O(in-flight
  fragments). `collated` O(name-order skew).
- **Threads:** bgzf (de)compression via `-t`. Parallel fragment scoring via `-S`
  (namesorted only; IO thread is always single-threaded).
- **Throughput:** ~300–600 M read-pairs/hour on a 16-core server with fast storage
  in `namesorted` sequential mode.

---

## License

GPL-3.0 — see [LICENSE](LICENSE).
