#!/usr/bin/env bash

# ==============================================================================
# Cross-Species Sequence Ambiguity Pipeline (Human <-> Mouse)
# Author: Bioinformatics & Systems Automation Specialist
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Logging Functions
# ------------------------------------------------------------------------------
log_info() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [\033[0;34mINFO\033[0m] $1"
}

log_warn() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [\033[0;33mWARN\033[0m] $1"
}

log_error() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [\033[0;31mERROR\033[0m] $1" >&2
}

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
  -o STR    Output directory for final configurations and BED files

Optional Arguments:
  -s INT    Step size for window generation (default: 1 for exhaustive sliding window)
            [Note: Setting -s equal to -l switches to non-overlapping tiles for speed]
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

# Validate required parameters
if [[ -z "$HUMAN_FA" || -z "$MOUSE_FA" || -z "$LEN" || -z "$MISMATCHES" || -z "$THREADS" || -z "$OUT_DIR" ]]; then
    log_error "Missing required arguments."
    usage
fi

# Convert paths to absolute paths
HUMAN_FA=$(realpath "$HUMAN_FA")
MOUSE_FA=$(realpath "$MOUSE_FA")
mkdir -p "$OUT_DIR"
OUT_DIR=$(realpath "$OUT_DIR")

# ------------------------------------------------------------------------------
# Environment & Dependency Checks
# ------------------------------------------------------------------------------
log_info "Verifying software dependencies..."
for cmd in samtools bedtools bowtie2 bowtie2-build awk; do
    if ! command -v "$cmd" &> /dev/null; then
        log_error "Required dependency '$cmd' is missing from PATH."
        exit 1
    fi
done
log_info "All dependencies verified successfully."

# ------------------------------------------------------------------------------
# Resource & Cleanup Management (Trap setup)
# ------------------------------------------------------------------------------
TMP_DIR=$(mktemp -d -p "$OUT_DIR" .tmp_ambiguity_XXXXXX)
cleanup() {
    log_info "Executing cleanup of intermediate temporary assets..."
    rm -rf "$TMP_DIR"
    log_info "Cleanup complete."
}
trap cleanup EXIT

# ------------------------------------------------------------------------------
# Pipeline Core Execution Function
# ------------------------------------------------------------------------------
run_ambiguity_mapping() {
    local SRC_FASTA="$1"
    local TGT_FASTA="$2"
    local SRC_LABEL="$3"
    local TGT_LABEL="$4"
    local FINAL_OUT_BED="$5"

    log_info "------------------------------------------------------------"
    log_info "Starting Processing Direction: ${SRC_LABEL} -> ${TGT_LABEL}"
    log_info "------------------------------------------------------------"

    # Step 1: Ensure Target Bowtie2 Index Exists
    local INDEX_DIR="${OUT_DIR}/indexes/${TGT_LABEL}"
    local INDEX_PREFIX="${INDEX_DIR}/${TGT_LABEL}_idx"
    mkdir -p "$INDEX_DIR"

    if ls "${INDEX_PREFIX}"*.bt2 >/dev/null 2>&1 || ls "${INDEX_PREFIX}"*.bt2l >/dev/null 2>&1; then
        log_info "Found existing Bowtie2 index for ${TGT_LABEL}."
    else
        log_info "Building Bowtie2 index for ${TGT_LABEL} (this may take some time)..."
        bowtie2-build --threads "$THREADS" "$TGT_FASTA" "$INDEX_PREFIX" > /dev/null
    fi

    # Step 2: Ensure Source FASTA is indexed
    if [[ ! -f "${SRC_FASTA}.fai" ]]; then
        log_info "Indexing source genome ${SRC_LABEL} with samtools..."
        samtools faidx "$SRC_FASTA"
    fi

    # Step 3: Stream Pipeline Execution
    # Calculate custom Bowtie2 penalty minimum score based on allowed mismatches
    # Mismatch penalty forced to constant 6. Perfect match is 0.
    local MIN_SCORE=$(( -6 * MISMATCHES ))

    log_info "Generating windows, matching, and extracting ambiguous regions via memory-safe stream..."
    
    # 1. Generate window coordinates on the fly
    # 2. Extract FASTA strings directly into standard output (attaches coordinates to header)
    # 3. Stream FASTA sequences directly into Bowtie2 standard input
    # 4. Drop unmapped reads instantly, keep only mapped candidates
    # 5. Filter precisely on Edit Distance (NM) and reverse engineer original coordinates
    bedtools makewindows -g "${SRC_FASTA}.fai" -w "$LEN" -s "$STEP" | \
    bedtools getfasta -fi "$SRC_FASTA" -bed stdin -fo - | \
    bowtie2 -x "$INDEX_PREFIX" \
        -f -U - \
        --threads "$THREADS" \
        --end-to-end \
        --very-sensitive \
        --mp 6,6 --np 6 --rdg 10,6 --rfg 10,6 \
        --score-min "L,${MIN_SCORE},0" \
        -k 1 2> "${TMP_DIR}/${SRC_LABEL}_bowtie2.log" | \
    samtools view -F 4 -@ "$THREADS" - | \
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
            # Handle complex chromosome names safely by anchoring to trailing coordinate format
            match(qname, /:[0-9]+-[0-9]+$/)
            if (RSTART > 0) {
                chrom = substr(qname, 1, RSTART-1)
                coord_part = substr(qname, RSTART+1)
                split(coord_part, coords, "-")
                print chrom, coords[1], coords[2]
            }
        }
    }' > "${TMP_DIR}/${SRC_LABEL}_raw.bed"

    # Step 4: Sort and Merge Overlapping Ambiguous Segments
    log_info "Sorting and consolidating ambiguous genomic intervals into final output..."
    sort -k1,1 -k2,2n -S 4G --parallel="$THREADS" "${TMP_DIR}/${SRC_LABEL}_raw.bed" | \
    bedtools merge -i stdin > "$FINAL_OUT_BED"

    log_info "Successfully completed direction: ${SRC_LABEL} -> ${TGT_LABEL}"
    log_info "Output saved to: ${FINAL_OUT_BED}"
}

# ------------------------------------------------------------------------------
# Orchestration Execution
# ------------------------------------------------------------------------------
log_info "Initializing Bi-directional Cross-Species Mapping Workflow"

# Direction A: Human -> Mouse
run_ambiguity_mapping "$HUMAN_FA" "$MOUSE_FA" "Human" "Mouse" "${OUT_DIR}/Human_ambiguous_to_Mouse.bed"

# Direction B: Mouse -> Human
run_ambiguity_mapping "$MOUSE_FA" "$HUMAN_FA" "Mouse" "Human" "${OUT_DIR}/Mouse_ambiguous_to_Human.bed"

log_info "=============================================================================="
log_info "Pipeline executed to completion. All final target artifacts successfully built."
log_info "=============================================================================="
