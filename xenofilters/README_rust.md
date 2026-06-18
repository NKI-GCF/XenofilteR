# xenofilter — Rust implementation

Fast, streaming read classifier for xenograft / PDX data.

Compares multiple aligned BAM streams (host vs graft vs ...) per fragment and decides:

- best alignment → -o / stdout
- worse but still usable → --filtered-output
- ambiguous → --ambiguous-output

## Features (compared to original XenofilteR)

- Streaming / low-memory mode when read names are sorted identically
- Optional variant-aware rescue from VCF
- Rust performance & safety
- Custom BAM tags with decision score (planned)
- @PG header line (planned)

When --add-decision-tag is used, winning records receive an aux tag:

  XF:C:<phred>   phred-scaled confidence that this assignment is correct
                 (≈ 4.34 × log-likelihood ratio to the second-best alignment)
                 Range: 0–255 (capped). Higher = stronger evidence.
                 0 is never written (ambiguous cases get no tag unless threshold=0).

## Quick start

cargo build --release
./target/release/xenofilter host.bam graft.bam -o best.bam -f mouse.bam -a ambiguous.bam

## See also

Original XenofilteR (R): https://github.com/NKI-GCF/XenofilteR
