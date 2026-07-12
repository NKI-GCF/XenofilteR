# xenofilters — Applications

## Currently supported

### PDX / xenograft read disambiguation

The primary use case. Human tumour tissue grafted into
immunocompromised mouse; WGS, WES, or RNA-seq reads assigned to
the correct reference. Extends to any host–graft pair:

| Graft   | Host          | Notes                          |
|---------|---------------|--------------------------------|
| Human   | Mouse         | Clinical PDX; most common      |
| Human   | Rat           | Cardiovascular, CNS models     |
| Human   | Zebrafish     | Larval xenograft; small genome |
| Mouse   | Rat           | Syngeneic cross-species models |
| Pig     | Human         | Xenotransplantation research   |

Reference genome availability determines feasibility. Zebrafish
and pig references are well-assembled; add `--stream-labels` for
clarity in multi-species output.

---

### Viral integration site mapping

When a retrovirus or DNA virus integrates into the host genome,
reads spanning the junction have one mate mapping to host and the
other to the viral reference. `--chimeric-pairs` routes both mates
to their respective outputs with `XC:Z:` tags.

**Currently working:**

- HPV (all high-risk strains): cervical, HNSC, anal carcinoma
- HBV: hepatocellular carcinoma integration sites
- HTLV-1: adult T-cell leukaemia proviral load and integration

**Setup:**

```bash
xenofilters human.bam hpv16.bam \
  --chimeric-pairs 0:1 \
  --stream-labels human hpv16 \
  --output human_best.bam --output hpv_best.bam \
  --ambiguous-output ambiguous.bam
```

Post-filter: `samtools view -d XC:hpv16 human_best.bam` extracts
junction-spanning reads for integration site inference.

---

### Three-stream: HPV + human tissue xenografted in mouse

The combination of viral integration and PDX modelling in the same
sample. Three alignment streams; chimeric pairs configured for
human↔HPV only; mouse competes normally.

```bash
xenofilters human.bam hpv.bam mouse.bam \
  --chimeric-pairs 0:1 \
  --stream-labels human hpv mouse
```

Assumption: viral integration events appear at low frequency;
most fragments assign to human or mouse via normal tournament.

---

### Cell line species authentication

Cell line misidentification (e.g. HeLa contamination) and inter-
species cross-contamination are widespread. Two streams — expected
species vs. contaminant panel — with low ambiguous threshold.
Contamination fraction estimated from routing counters.

Requires a contaminant reference. A `--stats-output` JSON provides
contamination percentage parseable by MultiQC or custom QC pipelines.

---

### Aligner comparison and benchmarking

Two or more streams from the same reads, aligned with different
aligners or parameter sets. xenofilters identifies reads where
aligners disagree (ambiguous or split routing) vs. where all agree.

Example: bwa-mem2 vs STAR on RNA-seq for splice-junction sensitivity.
Disagreed reads (routed to `--ambiguous-output`) represent the
informative subset for benchmarking.

Streams must be name-sorted; namesorted backend handles N ≤ 32.
No chimeric pairs; the `--ambiguous-threshold` controls sensitivity.

---

### Parental genome separation in F1 hybrids

Two streams: paternal and maternal reference (or two inbred strain
references). Reads assigned to whichever parent they align to more
accurately. With `--sample-variants` pointing to a phased VCF,
variant rescue improves assignment near heterozygous sites.

**Organisms with established inbred strain references:**

- Mouse (C57BL/6 × CAST/EiJ, B6 × DBA/2J; Sanger MGP VCFs available)
- Drosophila (DGRP lines)
- Arabidopsis (Col-0 × Ler)
- Yeast (S288C × SK1)

---

### Transposon / transgene vs. endogenous genome

A transgenic insert (Cre, GFP, human gene knock-in) aligns only to
the construct reference, not the host. One stream per genome. Reads
spanning the insertion boundary are chimeric.

Current limitation: the construct FASTA must be appended to one
reference or provided as a standalone FASTA for stream 2. No tool
change required.

---

### Microbiome host subtraction

Metagenomic or metatranscriptomic samples from host tissue contain
host reads that inflate sequencing cost and mask microbial signal.
Stream 0: host reference. Stream 1+: microbial pangenome or
pathogen reference. Host-winning reads discarded; microbial-winning
reads retained.

Current constraint: microbial reference must be a single FASTA.
A dereplicated pangenome (e.g. UHGG) works; query alignment must be
done first externally.

---

### RNA-seq contamination and ambient RNA removal

In single-nucleus RNA-seq, ambient RNA from lysed cells contaminates
droplets. If the experiment is a co-culture or xenograft, stream-level
ambient contamination can be estimated per barcode by treating each
species as a stream. Barcode-level routing counters are not yet
implemented (roadmap).

---

## Applications requiring minor code changes

### Long-read support (ONT / PacBio HiFi)

Current scoring assumes Illumina-style qualities and error profiles.
ONT R9/R10 and HiFi have distinct quality distributions; the penalty
arrays would need calibration. CIGAR strings are structurally
compatible. Chimeric-read detection is more important for long reads
(single-molecule chimeras from ligation artefacts).

Change required: quality-score penalty recalibration; optionally
a `--error-model` flag selecting Illumina, HiFi, or ONT presets.

---

### Bisulfite / methylation-aware scoring

WGBS reads contain C→T conversions at unmethylated CpGs. Mismatches
from bisulfite conversion should not be penalised equally to true
SNV mismatches. A strand-specific conversion-aware MD parser would
allow correct scoring.

Change required: new `BaseOp::BisulfiteConversion` variant in
`ScoreOpIter`; `--bisulfite` flag sets penalty to 0 for C→T on
the correct strand.

---

### Allele-specific expression quantification

Phased variant VCF already in `--phased-variants` (implemented but
experimental). Reads assigned to H1 or H2 at each heterozygous site
enable allele-specific count matrices. Requires per-read haplotype
tag in output (`HP:i:1` / `HP:i:2`), not yet emitted.

Change required: emit `HP` SAM tag from `scratch.read_haplotype`
in `write_record`.

---

### Single-cell doublet detection

In scRNA-seq, a doublet is a droplet capturing cells from two species.
Per-barcode routing counters (not per-fragment) would reveal barcodes
with reads routing to both species equally — doublets. Requires
barcode-aware counter accumulation keyed on `CB:Z:` tag.

Change required: barcode-level counter HashMap parallel to the
per-stream routing counters; `--barcode-tag CB` flag.

---

### CRISPR screen off-target scoring

Guide RNA off-targets align to both intended locus and off-target
loci. Two streams: on-target region FASTA vs. full genome with the
on-target masked. Reads that align better to the full genome are
candidate off-targets.

Change required: none to the disambiguation logic. Requires careful
reference construction externally. `--positive-regions` BED could
bias scoring toward expected cut sites.

---

### Structural variant breakpoint phasing

Long-read SVs produce chimeric alignments at breakpoints. With
`--chimeric-pairs`, reads spanning an SV junction are flagged
`XC:Z:`. If two streams represent the two alleles of a heterozygous
SV (constructed haplotype-resolved assemblies), xenofilters assigns
each read to its haplotype of origin.

Change required: none. Requires haplotype-resolved assemblies as
input (available for human from HPRC, others limited).

---

### Metagenome-assembled genome (MAG) dereplication

Multiple near-identical MAGs assembled from the same sample; reads
re-mapped to each. xenofilters assigns reads to the MAG they align
to most specifically. Reduces redundancy in MAG catalogues by
identifying reads that are genuinely specific to one assembly.

Change required: N > 2 streams in hashlookup (currently namesorted
only for N > 2). Memory scaling is the constraint.

---

## Speculative / exploratory

### Tumour clonal deconvolution

Two subclonal genomes (constructed from phased somatic variants) as
streams. Reads from one clone vs. the other. Variant rescue becomes
the primary disambiguation mechanism; NW scoring on somatic variants
distinguishes clones at heterozygous positions. Accuracy unknown
without empirical data. Requires high-confidence phased somatic VCFs.

### Horizontal gene transfer detection

Prokaryotic read aligning better to a recipient genome than to the
donor genus reference is a candidate HGT event. Multi-stream with
one stream per candidate donor clade. Computationally expensive
(many streams); CIGAR/MD subsumption helps prune early.

### Palaeogenomics: ancient DNA endogenous fraction

aDNA extracts contain modern environmental contaminants. Stream 0:
target ancient genome. Stream 1: modern human reference. Reads
routing to the ancient genome with low mismatch density (C→T
deamination artefacts correctly modelled) are endogenous. Requires
bisulfite-mode penalty adjustment (above) plus damage-pattern
scoring; not implemented.

### Reference genome quality assessment

Two assemblies of the same species (e.g. GRCh38 vs. T2T-CHM13).
Reads from a well-characterised sample assigned to whichever
assembly they align to with fewer mismatches. Regions consistently
ambiguous between assemblies flag collapsed repeats or assembly
errors. Aggregate ambiguous BED from `--ambiguous-output` and
`bedtools coverage`.
