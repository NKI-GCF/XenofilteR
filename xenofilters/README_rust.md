# xenofilter — Rust implementation

Fast, streaming read classifier
 - to assign xenograft / PDX data to host (or graft).
 - to compare two alignment algorithms.
 - to assign viral integration data to host or pathogen.
 - to assign viral integration in xenograft data to host, graft or pathogen.

Compares multiple aligned BAM streams (host vs graft vs ...) per fragment and decides:

- alignment is better than in other species → -o / stdout
- ambiguous is as good or too close to other species→ --ambiguous-output
- alignment is worse than in other species → --filtered-output

## Features (compared to original XenofilteR)

- Streaming / low-memory mode when read names are sorted identically
- Optional variant-aware rescue from VCF
- Rust performance & safety
- Custom BAM tags with decision score
- @PG header line

When --add-decision-tag is used, winning records receive an aux tag:

  XF:C:<phred>   phred-scaled confidence that this assignment is correct
                 (≈ 4.34 × log-likelihood ratio to the second-best alignment)
                 Range: 0–255 (capped). Higher = stronger evidence.
                 0 is never written (ambiguous cases get no tag unless threshold=0).

## Quick start

cargo build --release
./target/release/xenofilter host.bam graft.bam -o best.bam -f mouse.bam -a ambiguous.bam
