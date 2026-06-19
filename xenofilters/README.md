# xenofilters

Fast, streaming read classifier for xenograft / PDX sequencing data — written in Rust.

---

## What it does

Patient-derived xenograft (PDX) models graft human tumour tissue into mouse hosts.
Sequencing a PDX sample therefore captures reads from **both** species.
`xenofilters` takes two name-sorted BAM files — one aligned to the **graft** reference
(e.g. human GRCh38) and one to the **host** (e.g. mouse GRCm39) — and assigns each
read fragment to the species it aligns to best.

The same approach applies wherever you need to separate reads from two co-sequenced
genomes: viral integration studies, transposon experiments, dual-genome cell-line work,
or comparison of two different aligners on the same sample.

---

## Algorithm

For every fragment (one or two reads sharing the same name) `xenofilters`:

1. **Fast-path ordering** — if one alignment is mapped and the other is not, the mapped
   alignment wins immediately without any scoring.

2. **Perfect-match fast path** — if one alignment has a single CIGAR operation with no
   mismatches in the MD tag (a "perfect" alignment) and the other does not, the perfect
   alignment wins without scoring.

3. **Log-likelihood scoring** — when both alignments are imperfect, a per-base
   log-likelihood score is computed from the CIGAR + MD tag + base quality scores.
   Mismatches and soft-clips are penalised; indels are penalised with affine gap costs.
   The fragment score is the sum over all reads in the fragment.

4. **Variant-aware rescue** (optional, `-s` / `-p`) — if a VCF/BCF is supplied, known
   variants are scored via a Needleman–Wunsch alignment of the alt allele against the
   read.  A positive score delta (alt allele better supported than reference) is added to
   the fragment score, potentially rescuing alignments that would otherwise be filtered.
   Non-overlapping variants are combined with a weighted interval scheduling DP.

5. **Decision** — the stream with the higher score wins.  If the Phred-scaled score
   difference is below `--ambiguous-threshold` (default 0, i.e. disabled), the fragment
   is written to `--ambiguous-output` instead.

The optional `XF:C:<phred>` aux tag encodes the Phred-scaled confidence:
`XF = round(10 × |score_delta| / ln 10)`, capped at 255.

---

## Installation

### From source

```bash
git clone https://github.com/NKI-GCF/XenofilteR
cd XenofilteR
cargo build --release
# Binary is at ./target/release/xenofilters
```

### Via cargo install (once published)

```bash
cargo install xenofilters
```

**Requirements:** Rust 1.85+ (edition 2024). No external C libraries required.

---

## Quick start

```bash
# Align reads to both references first (name-sorted output required)
bwa mem -M -t 8 human_ref.fa reads_R1.fq.gz reads_R2.fq.gz \
  | samtools sort -n -@ 4 -o human.bam

bwa mem -M -t 8 mouse_ref.fa reads_R1.fq.gz reads_R2.fq.gz \
  | samtools sort -n -@ 4 -o mouse.bam

# Or use the provided bwa_mem.sh wrapper that does this in one step:
# ./bwa_mem.sh human_ref.fa mouse_ref.fa output.bam reads_R1.fq reads_R2.fq

# Classify
xenofilters human.bam mouse.bam \
  --output human_best.bam \
  --output mouse_best.bam \
  --filtered-output human_filtered.bam \
  --filtered-output mouse_filtered.bam \
  --ambiguous-output ambiguous.bam \
  --threads 8 \
  --verbose
```

---

## CLI reference

```
Usage: xenofilters [OPTIONS] <ALIGNMENT>...

Arguments:
  <ALIGNMENT>...
      Input alignments (name-sorted BAM). At least two required.

Options:
  -o, --output <PATH>...
      Output file for winning reads (one per alignment stream).
      Defaults to stdout for the first stream when omitted.

  -f, --filtered-output <PATH>...
      Output file for reads that lost (one per stream). Default: discard.

  -a, --ambiguous-output <PATH>...
      Output file for reads with equally good alignments. Default: discard.

  -O, --stdout-format <FORMAT>
      Output format for stdout [sam | bam | cram]. [default: sam]

  -t, --threads <N>
      bgzf (de)compression worker threads. [default: 4]

  -v, --verbose
      Increase verbosity (-v = INFO, -vv = DEBUG). Respects RUST_LOG.

  -U, --discard-unmapped
      Discard fragments that are unmapped in every stream, even from
      --filtered-output.

  -m, --mismatch-penalty <F>
      Mismatch penalty (Phred-scaled, positive). [default: 4]

  -g, --gap-open <F>
      Gap-open penalty (positive). [default: 6]

  -e, --gap-extend <F>
      Gap-extend penalty per base (positive). [default: 1]

  -c, --clipping-penalty <F>
      Soft-clip penalty per clipped base. [default: 5]

  -a, --ambiguous-threshold <PHRED>
      Minimum score difference (Phred) to call a winner. Pairs below this
      go to --ambiguous-output. 0 = disabled. [default: 0]

  -s, --sample-variants <[IDX:]FILE>...
      Sample-specific VCF/BCF for variant-aware rescue. Prefix with stream
      index (e.g. 0:graft.vcf) or assign positionally.

  -p, --population-variants <[IDX:]FILE>...
      Population-frequency VCF/BCF (INFO/AF field used). Same syntax.

  -A, --add-decision-tag
      Attach XF:C:<phred> aux tag to winning records.

  -P, --no-program-line
      Suppress the @PG header line in output BAMs.

  -R, --strip-read-suffix <MODE>
      Handle /1 /2 read-name suffixes [auto|true|false|variable].
      [default: auto]

      --single-alignment-mode
      Allow a single alignment stream (requires two variant profiles).

  -h, --help     Print help (see --help for full details)
  -V, --version  Print version
```

---

## Output files

| Stream | Reads written |
|--------|--------------|
| `--output` | Reads that aligned **better** to this reference |
| `--filtered-output` | Reads that aligned **worse** (won by another stream) |
| `--ambiguous-output` | Reads where no stream was clearly better |

All output files are bgzf-compressed BAM. The input BAM header (plus an
optional `@PG` line) is copied verbatim.

---

## Variant-aware rescue

Provide per-stream VCF/BCF files to enable rescue of reads that carry known
strain-specific variants. This is especially useful in single-alignment mode
to separate two strains from one alignment:

```bash
xenofilters strain_a.bam strain_b.bam \
  -s 0:strain_a_variants.vcf \
  -s 1:strain_b_variants.vcf \
  -o best_a.bam -o best_b.bam
```

Requirements for the VCF:
- Sample VCF: must contain `GT` and `GQ` FORMAT fields.
- Population VCF: must contain `AF` INFO field.
- BCF (binary VCF) is preferred for performance.
- Must be coordinate-sorted (variants within a stream are indexed for O(log n) lookup).

---

## Environment variables

| Variable | Effect |
|----------|--------|
| `RUST_LOG` | Override log level (e.g. `RUST_LOG=xenofilters=debug`) |

---

## Performance notes

- **Memory:** O(reads-per-fragment × streams). Typically < 100 MB for paired-end data.
- **Threads:** bgzf block (de)compression is parallelised with `-t`. The main
  record-processing loop is sequential by design (name-sorted streaming).
- **Throughput:** ~300–600 M read-pairs/hour on a 16-core server with fast storage.

---

## License

GPL-3.0 — see [LICENSE](LICENSE).
