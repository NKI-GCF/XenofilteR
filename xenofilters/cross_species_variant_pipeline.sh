#!/usr/bin/env bash

# ==============================================================================
# Monolithic Cross-Species Sequence Ambiguity & Variant Calling Pipeline
# Author: Bioinformatics & Systems Automation Specialist
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Logging Functions
# ------------------------------------------------------------------------------
log_info()  { echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [\033[0;34mINFO\033[0m] $1"; }
log_warn()  { echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [\033[0;33mWARN\033[0m] $1"; }
log_error() { echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [\033[0;31mERROR\033[0m] $1" >&2; }

# ------------------------------------------------------------------------------
# Usage & Help
# ------------------------------------------------------------------------------
usage() {
    cat << EOF
Usage: $0 -h <human_fasta> -m <mouse_fasta> -l <length> -e <mismatches> -t <threads> -o <out_dir> [-s <step_size>]

Required Arguments:
  -h STR    Path to Human reference FASTA
  -m STR    Path to Mouse reference FASTA
  -l INT    Sequence / K-mer length (L)
  -e INT    Allowed mismatches / edit distance (E)
  -t INT    Number of CPU threads to use
  -o STR    Output directory for final files

Optional Arguments:
  -s INT    Step size for window generation (default: 1 for exhaustive sliding window)
            [Pro-Tip: Set -s equal to -l for non-overlapping windows to speed up testing]
EOF
    exit 1
}

# ------------------------------------------------------------------------------
# Parameter Initialization & Parsing
# ------------------------------------------------------------------------------
HUMAN_FA=""
MOUSE_FA=""
LEN=""
MISMATCHES=""
THREADS=""
OUT_DIR=""
STEP=1

while getopts "h:m:l:e:t:o:s:" opt; do
    case "${opt}" in
        h) HUMAN_FA="$OPTARG" ;;
        m) MOUSE_FA="$OPTARG" ;;
        l) LEN="$OPTARG" ;;
        e) MISMATCHES="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        s) STEP="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$HUMAN_FA" || -z "$MOUSE_FA" || -z "$LEN" || -z "$MISMATCHES" || -z "$THREADS" || -z "$OUT_DIR" ]]; then
    log_error "Missing one or more required command-line arguments."
    usage
fi

# Standardize to absolute paths
HUMAN_FA=$(realpath "$HUMAN_FA")
MOUSE_FA=$(realpath "$MOUSE_FA")
mkdir -p "$OUT_DIR"
OUT_DIR=$(realpath "$OUT_DIR")

# ------------------------------------------------------------------------------
# Dependency Verification
# ------------------------------------------------------------------------------
log_info "Verifying execution dependencies in standard PATH..."
for cmd in samtools bedtools bcftools bowtie2 bowtie2-build awk; do
    if ! command -v "$cmd" &> /dev/null; then
        log_error "Critical pipeline dependency '$cmd' is missing."
        exit 1
    fi
done
log_info "All dependencies successfully verified."

# ------------------------------------------------------------------------------
# Resource Management & Emergency Cleanup Traps
# ------------------------------------------------------------------------------
TMP_DIR=$(mktemp -d -p "$OUT_DIR" .tmp_pipeline_XXXXXX)
cleanup() {
    log_info "Initiating cleanup sequence of intermediate workspace assets..."
    rm -rf "$TMP_DIR"
    log_info "Cleanup completed successfully."
}
trap cleanup EXIT INT TERM

# ------------------------------------------------------------------------------
# Core Execution Pipeline Engine
# ------------------------------------------------------------------------------
execute_pipeline_direction() {
    local SRC_FASTA="$1"
    local TGT_FASTA="$2"
    local SRC_LABEL="$3"
    local TGT_LABEL="$4"
    local OUT_BED_NAME="$5"
    local OUT_VCF_NAME="$6"

    log_info "============================================================"
    log_info "RUNNING DIRECTION: ${SRC_LABEL} Standard Sequences -> ${TGT_LABEL} Reference"
    log_info "============================================================"

    # --- Step 1: Reference Preparation & Indexing ---
    local INDEX_DIR="${OUT_DIR}/indexes/${TGT_LABEL}"
    local INDEX_PREFIX="${INDEX_DIR}/${TGT_LABEL}_idx"
    mkdir -p "$INDEX_DIR"

    if ls "${INDEX_PREFIX}"*.bt2 >/dev/null 2>&1 || ls "${INDEX_PREFIX}"*.bt2l >/dev/null 2>&1; then
        log_info "Valid Bowtie2 genome index detected for ${TGT_LABEL}."
    else
        log_info "Building new Bowtie2 index for ${TGT_LABEL} reference genome..."
        bowtie2-build --threads "$THREADS" "$TGT_FASTA" "$INDEX_PREFIX" > /dev/null
    fi

    if [[ ! -f "${SRC_FASTA}.fai" ]]; then
        log_info "Generating missing fasta index (.fai) for ${SRC_LABEL}..."
        samtools faidx "$SRC_FASTA"
    fi

    if [[ ! -f "${TGT_FASTA}.fai" ]]; then
        log_info "Generating missing fasta index (.fai) for ${TGT_LABEL}..."
        samtools faidx "$TGT_FASTA"
    fi

    # --- Step 2 & 4: K-mer Generation, Alignment & BAM Post-Processing ---
    # Setup flat penalty constraints to cleanly enforce explicit mismatch thresholds
    local MIN_SCORE=$(( -6 * MISMATCHES ))
    local SORTED_BAM="${TMP_DIR}/${SRC_LABEL}_on_${TGT_LABEL}.sorted.bam"

    log_info "Streaming ${SRC_LABEL} windows into ${TGT_LABEL} aligner with RG tags..."
    
    bedtools makewindows -g "${SRC_FASTA}.fai" -w "$LEN" -s "$STEP" | \
    bedtools getfasta -fi "$SRC_FASTA" -bed stdin -fo - | \
    bowtie2 -x "$INDEX_PREFIX" \
        -f -U - \
        --threads "$THREADS" \
        --end-to-end \
        --very-sensitive \
        --mp 6,6 --np 6 --rdg 10,6 --rfg 10,6 \
        --score-min "L,${MIN_SCORE},0" \
        --rg-id "RG_${SRC_LABEL}_to_${TGT_LABEL}" \
        --rg "SM:${SRC_LABEL}_sequences" \
        --rg "PL:ILLUMINA" 2> "${TMP_DIR}/${SRC_LABEL}_alignment.log" | \
    samtools view -h -F 4 -@ "$THREADS" - | \
    samtools sort -@ "$THREADS" -m 2G -T "${TMP_DIR}/sort_cache" -o "$SORTED_BAM" -

    log_info "Generating coordinate lookup index for the compiled BAM architecture..."
    samtools index -@ "$THREADS" "$SORTED_BAM"

    # --- Step 3: Ambiguity Filtering & Source BED Generation ---
    log_info "De-multiplexing read headers to build structural source BED files..."
    
    samtools view "$SORTED_BAM" | \
    awk -v max_e="$MISMATCHES" '
    BEGIN { FS="\t"; OFS="\t" }
    {
        nm = -1
        for (i=12; i<=NF; i++) {
            if ($i ~ /^NM:i:/) {
                split($i, a, ":")
                nm = a[3]
                break
            }
        }
        if (nm >= 0 && nm <= max_e) {
            qname = $1
            match(qname, /:[0-9]+-[0-9]+$/)
            if (RSTART > 0) {
                chrom = substr(qname, 1, RSTART-1)
                coord_part = substr(qname, RSTART+1)
                split(coord_part, coords, "-")
                print chrom, coords[1], coords[2]
            }
        }
    }' | sort -k1,1 -k2,2n -S 2G --parallel="$THREADS" | \
    bedtools merge -i stdin > "${OUT_DIR}/${OUT_BED_NAME}"

    # --- Step 5: Cross-Species Variant Calling ---
    log_info "Executing variant calling matrix via bcftools against ${TGT_LABEL} genome..."
    
    bcftools mpileup --threads "$THREADS" -f "$TGT_FASTA" -Ou "$SORTED_BAM" 2> /dev/null | \
    bcftools call --threads "$THREADS" -mv -Ov -o "${OUT_DIR}/${OUT_VCF_NAME}"

    log_info "Completed pipeline operations for direction: ${SRC_LABEL} -> ${TGT_LABEL}"
}

# ------------------------------------------------------------------------------
# Execution Management Block
# ------------------------------------------------------------------------------
log_info "Initializing Bi-Directional Cross-Species Mappability & Variant Genotyping Workflow"

# Direction 1: Human -> Mouse
execute_pipeline_direction \
    "$HUMAN_FA" \
    "$MOUSE_FA" \
    "human" \
    "mouse" \
    "human_ambiguous_zones.bed" \
    "human_on_mouse_variants.vcf"

# Direction 2: Mouse -> Human
execute_pipeline_direction \
    "$MOUSE_FA" \
    "$HUMAN_FA" \
    "mouse" \
    "human" \
    "mouse_ambiguous_zones.bed" \
    "mouse_on_human_variants.vcf"
