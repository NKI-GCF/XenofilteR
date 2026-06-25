#!/usr/bin/env bash
# ==============================================================================
# cross_species_variant_pipeline.sh  v1.0.0
#
# PURPOSE
#   Characterise cross-species sequence ambiguity between Human and Mouse
#   genomes in both directions:
#
#     Human  → Mouse : human_ambiguous_zones.bed  +  human_on_mouse_variants.vcf
#     Mouse  → Human : mouse_ambiguous_zones.bed  +  mouse_on_human_variants.vcf
#
# PIPELINE STAGES (per direction)
#   1. Reference FASTA indexing (samtools faidx, aligner index)
#   2. Streaming k-mer FASTA generation (embedded Python, named FIFO)
#   3. Alignment to target genome with injected @RG header
#   4. Ambiguity filtering → source-coordinate BED (sort + bedtools merge)
#   5. BAM post-processing  (filter mapped, sort, index)
#   6. Cross-species variant calling (bcftools mpileup | bcftools call)
#
# USAGE
#   ./cross_species_variant_pipeline.sh \
#       -h hg38.fa -m mm39.fa -l 150 -e 2 -t 16 -o results/
#
# REQUIREMENTS
#   samtools >= 1.10, bedtools >= 2.25, bcftools >= 1.10,
#   bowtie2  >= 2.4  (or set ALIGNER=bwa-mem2 below),
#   python3  >= 3.6
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

# ─── Metadata ─────────────────────────────────────────────────────────────────
readonly SCRIPT="$(basename "$0")"
readonly VERSION="1.0.0"

# ─── Tuneable defaults (override with env vars before calling) ─────────────────
ALIGNER="${ALIGNER:-bowtie2}"    # bowtie2 | bwa-mem2
MAX_N_PCT="${MAX_N_PCT:-50}"     # skip windows with > this % N bases

# ─── CLI parameters ───────────────────────────────────────────────────────────
HUMAN_REF=""
MOUSE_REF=""
READ_LEN=150
EDIT_DIST=2
THREADS=4
OUT_DIR=""

# ─── Runtime state ────────────────────────────────────────────────────────────
declare -a _TMP_FILES=()
declare -a _TMP_DIRS=()
STREAM_PY=""      # path of the embedded Python k-mer streamer (created at runtime)
IDX_DIR=""        # <OUT_DIR>/indices

# ==============================================================================
#  LOGGING
# ==============================================================================
_ts()    { date '+%H:%M:%S'; }
info()   { printf "[INFO ] [%s] %s\n"  "$(_ts)" "$*"; }
warn()   { printf "[WARN ] [%s] %s\n"  "$(_ts)" "$*" >&2; }
error()  { printf "[ERROR] [%s] %s\n"  "$(_ts)" "$*" >&2; }
banner() { printf "\n[STEP ] [%s] ══ %s ══\n" "$(_ts)" "$*"; }
done_()  { printf "[DONE ] [%s] %s\n"  "$(_ts)" "$*"; }

# ==============================================================================
#  CLEANUP & SIGNAL TRAPS
# ==============================================================================
_cleanup() {
    local rc=$?
    info "Cleaning up temporary files…"
    local f d
    for f in "${_TMP_FILES[@]+"${_TMP_FILES[@]}"}"; do
        [[ -e "$f" ]] && rm -f  "$f" && info "  removed $f"
    done
    for d in "${_TMP_DIRS[@]+"${_TMP_DIRS[@]}"}"; do
        [[ -d "$d" ]] && rm -rf "$d" && info "  removed dir $d"
    done
    [[ $rc -ne 0 ]] && error "Pipeline exited with status $rc."
    exit $rc
}
trap '_cleanup'         EXIT
trap 'error "SIGINT";   exit 130' INT
trap 'error "SIGTERM";  exit 143' TERM

_mk_tmp()  {
    local f; f=$(mktemp "${TMPDIR:-/tmp}/xvar_XXXXXXXX${1:-.tmp}")
    _TMP_FILES+=("$f"); echo "$f"
}
_mk_tmpd() {
    local d; d=$(mktemp -d "${TMPDIR:-/tmp}/xvar_dir_XXXXXXXX")
    _TMP_DIRS+=("$d"); echo "$d"
}
# FIFO gets its own private directory to avoid TOCTOU races.
_mk_fifo() {
    local d; d=$(mktemp -d "${TMPDIR:-/tmp}/xvar_fifo_XXXXXXXX")
    _TMP_DIRS+=("$d")
    mkfifo "${d}/pipe"
    echo "${d}/pipe"
}

# ==============================================================================
#  USAGE
# ==============================================================================
usage() {
    cat <<EOF
$SCRIPT v$VERSION — Cross-species ambiguity + variant calling pipeline

USAGE
  $SCRIPT -h <human.fa> -m <mouse.fa> -l <int> -e <int> -t <int> -o <outdir>

REQUIRED
  -h <path>   Human reference FASTA  (must be samtools faidx-able)
  -m <path>   Mouse reference FASTA
  -l <int>    k-mer / read length
  -e <int>    Allowed mismatches (edit distance) during alignment
  -t <int>    CPU threads
  -o <path>   Output directory (created if absent)

ENVIRONMENT VARIABLES (optional overrides)
  ALIGNER     bowtie2 (default) | bwa-mem2
  MAX_N_PCT   Maximum %N bases per k-mer window to skip  (default: 50)

OUTPUTS
  <outdir>/human_ambiguous_zones.bed      Human regions mapping to Mouse
  <outdir>/mouse_ambiguous_zones.bed      Mouse regions mapping to Human
  <outdir>/human_on_mouse_variants.vcf    Variants vs Mouse ref (human seqs)
  <outdir>/mouse_on_human_variants.vcf    Variants vs Human ref (mouse seqs)

ALIGNMENT FLAGS (bowtie2, per direction)
  --end-to-end          Full-length k-mer alignment; no soft-clipping.
  -N 1 -L 20 -D 20 -R 3 -i S,1,0.50
                        High-sensitivity seed search.
  --score-min C,-(6×e),0
                        Constant minimum score; admits ≤ e mismatches (--mp 6,6).
  --mp 6,6              Uniform mismatch penalty; makes the count exact.
  --rdg 1000,1000 --rfg 1000,1000
                        Gap open = 1000; effectively disables gapped alignments.
                        Edit-distance therefore equals Hamming distance.
  --no-unal             Suppress unaligned reads (saves pipeline bandwidth).
  -k 1                  One alignment per read (sufficient for ambiguity / VC).
  --rg-id <id> --rg <field>
                        Injects a valid @RG header; required by bcftools.

EXAMPLES
  # 150-mer windows, 2 mismatches, 16 threads
  $SCRIPT -h hg38.fa -m mm39.fa -l 150 -e 2 -t 16 -o results/

  # Use bwa-mem2 instead of bowtie2
  ALIGNER=bwa-mem2 $SCRIPT -h hg38.fa -m mm39.fa -l 150 -e 2 -t 16 -o results/

EOF
    exit 0
}

# ==============================================================================
#  ARGUMENT PARSING
# ==============================================================================
parse_args() {
    [[ $# -eq 0 ]] && usage
    local opt
    while getopts ":h:m:l:e:t:o:" opt; do
        case "$opt" in
            h) HUMAN_REF="$OPTARG"  ;;
            m) MOUSE_REF="$OPTARG"  ;;
            l) READ_LEN="$OPTARG"   ;;
            e) EDIT_DIST="$OPTARG"  ;;
            t) THREADS="$OPTARG"    ;;
            o) OUT_DIR="$OPTARG"    ;;
            :) error "Option -$OPTARG requires an argument."; exit 1 ;;
           \?) error "Unknown option: -$OPTARG"; usage ;;
        esac
    done

    # Required fields
    local errs=0
    [[ -z "$HUMAN_REF" ]] && { error "Missing -h (human FASTA)";    (( ++errs )); }
    [[ -z "$MOUSE_REF" ]] && { error "Missing -m (mouse FASTA)";    (( ++errs )); }
    [[ -z "$OUT_DIR"   ]] && { error "Missing -o (output dir)";     (( ++errs )); }
    [[ -z "$READ_LEN"  ]] && { error "Missing -l (read length)";    (( ++errs )); }
    [[ -z "$EDIT_DIST" ]] && { error "Missing -e (edit distance)";  (( ++errs )); }
    [[ -z "$THREADS"   ]] && { error "Missing -t (threads)";        (( ++errs )); }
    [[ $errs -gt 0 ]] && exit 1

    [[ -f "$HUMAN_REF" ]] || { error "Not found: $HUMAN_REF"; exit 1; }
    [[ -f "$MOUSE_REF" ]] || { error "Not found: $MOUSE_REF"; exit 1; }

    for var in READ_LEN EDIT_DIST THREADS; do
        [[ "${!var}" =~ ^[0-9]+$ && "${!var}" -ge 1 ]] || {
            error "$var must be a positive integer (got: ${!var})"; exit 1; }
    done
    [[ "$EDIT_DIST" -lt "$READ_LEN" ]] || { error "-e must be < -l"; exit 1; }
    [[ "$ALIGNER" == "bowtie2" || "$ALIGNER" == "bwa-mem2" ]] || {
        error "ALIGNER must be 'bowtie2' or 'bwa-mem2'"; exit 1; }

    IDX_DIR="${OUT_DIR}/indices"
    mkdir -p "$OUT_DIR" "$IDX_DIR"

    info "═══ Parameters ══════════════════════════════════════"
    info "  Human ref     : $HUMAN_REF"
    info "  Mouse ref     : $MOUSE_REF"
    info "  Read length   : $READ_LEN"
    info "  Edit distance : $EDIT_DIST"
    info "  Threads       : $THREADS"
    info "  Aligner       : $ALIGNER"
    info "  Output dir    : $OUT_DIR"
    info "  Max %N/window : $MAX_N_PCT%"
    info "════════════════════════════════════════════════════"
}

# ==============================================================================
#  DEPENDENCY CHECK
# ==============================================================================
check_deps() {
    banner "Dependency verification"
    local tools=("samtools" "bedtools" "bcftools" "python3" "awk" "sort")
    [[ "$ALIGNER" == "bowtie2"  ]] && tools+=("bowtie2" "bowtie2-build")
    [[ "$ALIGNER" == "bwa-mem2" ]] && tools+=("bwa-mem2")

    local miss=0 t
    for t in "${tools[@]}"; do
        if command -v "$t" &>/dev/null; then info "  ✓  $t"
        else error "  ✗  $t — NOT FOUND"; (( ++miss )); fi
    done
    [[ $miss -gt 0 ]] && { error "Install missing tools and re-run."; exit 1; }
    done_ "All ${#tools[@]} dependencies present."
}

# ==============================================================================
#  REFERENCE FASTA INDEXING
# ==============================================================================
ensure_fai() {
    local ref="$1"
    if [[ ! -f "${ref}.fai" ]]; then
        info "samtools faidx $ref"
        samtools faidx "$ref"
        done_ ".fai created for $(basename "$ref")"
    else
        info "  .fai OK: ${ref}.fai"
    fi
}

# ==============================================================================
#  ALIGNER INDEX
# ==============================================================================
_idx_base() { basename "$1" | sed 's/\.\(fa\|fna\|fasta\)\(\.gz\)\?$//'; }

ensure_bt2_index() {
    local ref="$1"
    local pfx="${IDX_DIR}/$(_idx_base "$ref")"
    if ls "${pfx}".*.bt2 &>/dev/null 2>&1 || ls "${pfx}".*.bt2l &>/dev/null 2>&1; then
        info "  Bowtie2 index OK: ${pfx}.*"
    else
        info "Building Bowtie2 index for $(basename "$ref") → $pfx"
        info "  (this may take 30–60 min)"
        bowtie2-build --threads "$THREADS" --quiet "$ref" "$pfx"
        done_ "Bowtie2 index: $pfx"
    fi
    echo "$pfx"
}

ensure_bwamem2_index() {
    local ref="$1"
    local dest="${IDX_DIR}/$(_idx_base "$ref").fa"
    if [[ -f "${dest}.0123" ]]; then
        info "  bwa-mem2 index OK: ${dest}.*"
    else
        info "Building bwa-mem2 index for $(basename "$ref") → $dest"
        info "  (this may take 60–90 min and ~24 GB RAM)"
        [[ -f "$dest" ]] || ln -sf "$(realpath "$ref")" "$dest"
        bwa-mem2 index "$dest"
        done_ "bwa-mem2 index: ${dest}.*"
    fi
    echo "$dest"
}

# ==============================================================================
#  EMBEDDED PYTHON K-MER STREAMER
# ==============================================================================
create_stream_py() {
    STREAM_PY=$(_mk_tmp ".py")
    cat > "$STREAM_PY" << 'PYEOF'
#!/usr/bin/env python3
"""
Stream k-mer FASTA windows from a .fai-indexed FASTA file.

Usage:
    stream_fasta.py <fasta> <chrom> <L> <step> <max_n_pct> <rg_id>

Output (to stdout):
    FASTA records.
    Header: >{chrom}@@{start}@@{end}  (0-based, half-open / BED convention)

Design:
    - Navigates the FASTA file using .fai byte offsets; never loads an
      entire chromosome into memory.
    - Windows where N-fraction > max_n_pct / 100 are silently skipped.
    - BrokenPipeError is caught cleanly (aligner may close the pipe early).
"""

import os
import sys


def load_fai(path):
    fai = path + '.fai'
    if not os.path.exists(fai):
        sys.exit(f'FATAL: {fai} missing — run: samtools faidx {path}')
    index = {}
    with open(fai) as fh:
        for line in fh:
            c = line.rstrip('\n').split('\t')
            if len(c) >= 5:
                index[c[0]] = (int(c[1]), int(c[2]), int(c[3]), int(c[4]))
    return index


def fetch(fh, offset, bpl, bpl_b, start, length):
    """Read `length` bases at 0-based `start`, handling multi-line FASTA."""
    nl   = bpl_b - bpl
    fl, r = divmod(start, bpl)
    fh.seek(offset + fl * bpl_b + r)
    buf = bytearray(length)
    n   = 0
    pos = r
    while n < length:
        want  = min(bpl - pos, length - n)
        chunk = fh.read(want)
        if not chunk:
            break
        k = len(chunk)
        buf[n:n+k] = chunk
        n  += k
        pos = 0
        if n < length:
            fh.seek(nl, 1)
    return bytes(buf[:n])


def stream(fasta, chrom, L, step, max_n_pct):
    idx = load_fai(fasta)
    if chrom not in idx:
        print(f'WARNING: {chrom} absent from index — skipped.', file=sys.stderr)
        return
    clen, offset, bpl, bpl_b = idx[chrom]
    max_n = int(L * max_n_pct / 100)
    out   = sys.stdout.buffer
    with open(fasta, 'rb') as fh:
        pos = 0
        while pos + L <= clen:
            seq = fetch(fh, offset, bpl, bpl_b, pos, L).upper()
            if len(seq) == L and seq.count(b'N') <= max_n:
                name = f'{chrom}@@{pos}@@{pos+L}'.encode()
                out.write(b'>'); out.write(name); out.write(b'\n')
                out.write(seq);                   out.write(b'\n')
            pos += step
    out.flush()


if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.exit(f'Usage: {sys.argv[0]} <fasta> <chrom> <L> <step> <max_n_pct>')
    try:
        stream(sys.argv[1], sys.argv[2],
               int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
    except BrokenPipeError:
        pass
PYEOF
    chmod +x "$STREAM_PY"
}

# ==============================================================================
#  SAM → SOURCE-COORDINATE BED  (AWK; POSIX-portable)
# ==============================================================================
# Read name format: {chrom}@@{start}@@{end}   (0-based, half-open)
# Flag checks use integer arithmetic — no gawk and() extension required.
readonly _AWK_BED='
function has_flag(f,b) { return int(f/b) % 2 == 1 }
BEGIN { OFS="\t" }
/^@/ { next }
{
    f = int($2)
    if (has_flag(f,4))    next   # unmapped
    if (has_flag(f,256))  next   # secondary
    if (has_flag(f,2048)) next   # supplementary
    n = $1
    gsub(/@@/, "\036", n)
    m = split(n, a, "\036")
    if (m >= 3) print a[1], a[2], a[3]
}
'

# ==============================================================================
#  CORE: RUN ONE PIPELINE DIRECTION
# ==============================================================================
#
# Arguments:
#   $1  src_label   "human" | "mouse"
#   $2  tgt_label   "mouse" | "human"
#   $3  src_fa      source genome FASTA
#   $4  tgt_fa      target genome FASTA  (used by bcftools)
#   $5  tgt_idx     aligner index prefix / FASTA path for target
#   $6  out_bed     output BED (source coordinates)
#   $7  out_vcf     output VCF (target coordinates)
#
run_direction() {
    local src_label="$1"
    local tgt_label="$2"
    local src_fa="$3"
    local tgt_fa="$4"
    local tgt_idx="$5"
    local out_bed="$6"
    local out_vcf="$7"

    banner "Direction: $src_label → $tgt_label"
    info "Source FASTA  : $src_fa"
    info "Target FASTA  : $tgt_fa"
    info "Target index  : $tgt_idx"
    info "Output BED    : $out_bed"
    info "Output VCF    : $out_vcf"

    # ── Score thresholds ────────────────────────────────────────────────────
    local bt2_score_min="C,$(( -6 * EDIT_DIST )),0"
    local bwa_min_score=$(( READ_LEN - 4 * EDIT_DIST ))

    # ── Read Group fields (required by bcftools / GATK) ─────────────────────
    # ID  : direction-specific unique identifier
    # SM  : sample name (source species label)
    # PL  : platform
    # LB  : library
    local rg_id="${src_label}_on_${tgt_label}"
    local rg_sm="${src_label^}_kmer_L${READ_LEN}"   # e.g. Human_kmer_L150
    local rg_pl="ILLUMINA"
    local rg_lb="${rg_sm}"

    # ── Per-direction work directory ─────────────────────────────────────────
    local work_dir; work_dir=$(_mk_tmpd)

    # Collect per-chromosome BED fragments
    local chr_bed_list=()

    # Sorted BAM (built across all chromosomes, used for variant calling)
    local sorted_bam="${work_dir}/${src_label}_on_${tgt_label}.sorted.bam"

    # Temporary file for aggregating aligned SAM records before sorting
    # We collect raw SAM text through a named pipe to avoid disk writes
    local sam_pipe; sam_pipe=$(_mk_fifo)
    local sam_header_file; sam_header_file=$(_mk_tmp ".sam_header")

    # We need a valid SAM header (from one alignment call) and all aligned
    # records merged, then sorted.  The strategy is:
    #
    #   FOR EACH CHROMOSOME:
    #     Python → FIFO → aligner (SAM to stdout)
    #       ├─ awk  → per-chromosome source-coord BED
    #       └─ pass-through → temp file (raw SAM records only, no @-lines)
    #
    #   AFTER ALL CHROMOSOMES:
    #     cat header + all raw SAM records | samtools sort → sorted BAM
    #
    # Because the aligner emits a fresh header per invocation, we capture
    # the header from the FIRST chromosome only, then strip @-lines from
    # the rest.

    local raw_sam_records; raw_sam_records=$(_mk_tmp ".raw.sam")
    local got_header=false
    local header_captured=false

    info "─── Stage 1: k-mer alignment + BED extraction (per chromosome) ───"

    while IFS=$'\t' read -r chrom chrom_len _rest; do
        [[ "$chrom_len" -lt "$READ_LEN" ]] && {
            info "  ⊘  $chrom (len=$chrom_len < L=$READ_LEN) — skipped"
            continue
        }

        local n_windows=$(( ( chrom_len - READ_LEN ) / READ_LEN + 1 ))
        info "  ▶  $chrom  len=${chrom_len}  ~${n_windows} windows"

        local fifo;    fifo=$(_mk_fifo)
        local aln_log; aln_log=$(_mk_tmp ".aln.log")
        local chr_bed="${work_dir}/${chrom}.bed"

        # ── Background: Python FASTA generator → FIFO ──────────────────────
        python3 "$STREAM_PY" \
            "$src_fa" "$chrom" "$READ_LEN" "$READ_LEN" "$MAX_N_PCT" \
            > "$fifo" 2>>"$aln_log" &
        local gen_pid=$!

        # ── Foreground: align → tee to BED awk + raw SAM accumulator ────────
        local n_ambig=0
        case "$ALIGNER" in
            bowtie2)
                n_ambig=$(
                    bowtie2 \
                        --end-to-end \
                        -N 1 -L 20 -D 20 -R 3 -i "S,1,0.50" \
                        --no-unal \
                        -k 1 \
                        --score-min "$bt2_score_min" \
                        --mp "6,6" \
                        --rdg "1000,1000" --rfg "1000,1000" \
                        -p "$THREADS" \
                        -x "$tgt_idx" \
                        -f -U "$fifo" \
                        --rg-id "$rg_id" \
                        --rg "SM:${rg_sm}" \
                        --rg "PL:${rg_pl}" \
                        --rg "LB:${rg_lb}" \
                        2>>"$aln_log" | \
                    tee >( awk "$_AWK_BED" >> "$chr_bed" ) | \
                    awk '
                        /^@/ {
                            if (ENVIRON["_HDR_DONE"] != "1") print > ENVIRON["_HDR_FILE"]
                            next
                        }
                        { print }
                    ' >> "$raw_sam_records"
                    wc -l < "$chr_bed" || echo 0
                ) || true
                ;;
            bwa-mem2)
                # bwa-mem2: inject RG via -R; @-line format is slightly different
                local rg_str="@RG\tID:${rg_id}\tSM:${rg_sm}\tPL:${rg_pl}\tLB:${rg_lb}"
                n_ambig=$(
                    bwa-mem2 mem \
                        -t "$THREADS" \
                        -k 11 -T "$bwa_min_score" \
                        -v 1 \
                        -R "$rg_str" \
                        "$tgt_idx" "$fifo" \
                        2>>"$aln_log" | \
                    tee >( awk "$_AWK_BED" >> "$chr_bed" ) | \
                    awk '
                        /^@/ {
                            if (ENVIRON["_HDR_DONE"] != "1") print > ENVIRON["_HDR_FILE"]
                            next
                        }
                        { print }
                    ' >> "$raw_sam_records"
                    wc -l < "$chr_bed" || echo 0
                ) || true
                ;;
        esac

        # Capture header from first chromosome only
        if ! $header_captured && [[ -s "${work_dir}/header.sam" ]]; then
            mv "${work_dir}/header.sam" "$sam_header_file"
            header_captured=true
        fi

        # ── Reap generator ──────────────────────────────────────────────────
        wait "$gen_pid" || {
            local rc=$?
            [[ $rc -eq 141 || $rc -eq 0 ]] || {
                error "FASTA generator failed for $chrom (exit $rc)"
                exit "$rc"
            }
        }

        [[ -s "$chr_bed" ]] && chr_bed_list+=("$chr_bed")
        info "      ambiguous windows: $n_ambig"

    done < "${src_fa}.fai"

    # ── Stage 2: Sort + merge source-coordinate BED ───────────────────────
    banner "Stage 2 [$src_label→$tgt_label]: Merge ambiguous BED"
    if [[ ${#chr_bed_list[@]} -eq 0 ]]; then
        warn "No ambiguous windows found ($src_label → $tgt_label)."
        touch "$out_bed"
    else
        info "Merging ${#chr_bed_list[@]} chromosome BED files → $out_bed"
        cat "${chr_bed_list[@]}" \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i stdin \
            > "$out_bed"
        local n_reg; n_reg=$(wc -l < "$out_bed")
        local bp;    bp=$(awk '{s+=$3-$2} END{print s+0}' "$out_bed")
        done_ "$out_bed  ($n_reg regions, ${bp} bp)"
    fi

    # ── Stage 3: Build sorted BAM ─────────────────────────────────────────
    banner "Stage 3 [$src_label→$tgt_label]: Sort & index BAM"

    # The awk in the alignment loop uses ENVIRON variables to write the header.
    # However, the tee+awk+ENVIRON trick is fragile across all shells.
    # We use a cleaner approach: re-run the alignment on chr1 only to capture
    # the SAM header, then merge header + raw records into the sorted BAM.

    if [[ -s "$raw_sam_records" ]]; then
        info "Assembling SAM header for $sorted_bam ..."

        # Use samtools view on a quick dummy to get the reference header,
        # or extract from the first alignment if available.
        # Most reliable: regenerate header from the reference FASTA dict.
        local sam_header; sam_header=$(_mk_tmp ".header.sam")

        # Build a SAM header from the target FASTA dictionary
        # @HD line
        printf '@HD\tVN:1.6\tSO:unsorted\n' > "$sam_header"
        # @SQ lines from .fai
        awk 'BEGIN{OFS="\t"} {print "@SQ","SN:"$1,"LN:"$2}' "${tgt_fa}.fai" \
            >> "$sam_header"
        # @RG line (required by bcftools)
        printf '@RG\tID:%s\tSM:%s\tPL:%s\tLB:%s\n' \
            "$rg_id" "$rg_sm" "$rg_pl" "$rg_lb" >> "$sam_header"
        # @PG line
        printf '@PG\tID:%s\tPN:%s\tVN:%s\tCL:%s\n' \
            "$ALIGNER" "$ALIGNER" "unknown" \
            "${ALIGNER} cross-species k-mer alignment" >> "$sam_header"

        info "Sorting aligned reads → $sorted_bam"
        {
            cat "$sam_header"
            cat "$raw_sam_records"
        } | samtools sort \
            -@ "$THREADS" \
            -m 2G \
            -O BAM \
            -o "$sorted_bam" -

        samtools index -@ "$THREADS" "$sorted_bam"
        done_ "Sorted BAM: $sorted_bam  ($(du -sh "$sorted_bam" | cut -f1))"
    else
        warn "No aligned SAM records found — skipping BAM creation."
        touch "${sorted_bam}.empty"
    fi

    # ── Stage 4: Cross-species variant calling ────────────────────────────
    banner "Stage 4 [$src_label→$tgt_label]: Variant calling → $out_vcf"

    if [[ -f "$sorted_bam" && -s "$sorted_bam" ]]; then
        info "bcftools mpileup | bcftools call → $out_vcf"
        info "  Reference : $tgt_fa"
        info "  BAM       : $sorted_bam"

        # bcftools mpileup:
        #   -f  : reference FASTA (target species)
        #   -d  : maximum per-position depth (prevents memory explosion on
        #          highly repetitive / collapsed regions)
        #   -Q  : skip bases with base quality < 20
        #   -q  : skip alignments with mapping quality < 20
        #   --annotate FORMAT/AD,FORMAT/DP
        #         : add allele depth and per-sample depth to FORMAT
        #
        # bcftools call:
        #   -m  : multiallelic caller (more accurate than -c for SNPs + indels)
        #   -v  : output variant sites only
        #   -A  : keep all alternate alleles (important for cross-species SNPs
        #          where the "alt" is just the source species' allele)
        #   --ploidy 1 : treat as haploid (k-mer sequences are single-copy)
        #
        bcftools mpileup \
            --threads "$THREADS" \
            -f "$tgt_fa" \
            -d 5000 \
            -Q 20 \
            -q 20 \
            --annotate "FORMAT/AD,FORMAT/DP" \
            -Ou \
            "$sorted_bam" \
        | bcftools call \
            --threads "$THREADS" \
            -m \
            -v \
            -A \
            --ploidy 1 \
            -Ov \
            -o "$out_vcf"

        local n_var; n_var=$(grep -vc '^#' "$out_vcf" || echo 0)
        done_ "$out_vcf  ($n_var variant sites)"
    else
        warn "No BAM available — skipping variant calling."
        touch "$out_vcf"
    fi

    done_ "Direction $src_label → $tgt_label complete."
}

# ==============================================================================
#  MAIN
# ==============================================================================
main() {
    parse_args "$@"

    banner "Cross-Species Ambiguity & Variant Pipeline  v${VERSION}"
    info "Human ↔ Mouse  |  L=${READ_LEN}  E=${EDIT_DIST}  aligner=${ALIGNER}"

    check_deps

    # ── Reference FASTA indices ───────────────────────────────────────────────
    banner "Reference FASTA indexing"
    ensure_fai "$HUMAN_REF"
    ensure_fai "$MOUSE_REF"

    # ── Aligner indices ───────────────────────────────────────────────────────
    banner "Aligner index preparation"
    local human_idx mouse_idx
    case "$ALIGNER" in
        bowtie2)
            human_idx=$(ensure_bt2_index "$HUMAN_REF")
            mouse_idx=$(ensure_bt2_index "$MOUSE_REF")
            ;;
        bwa-mem2)
            human_idx=$(ensure_bwamem2_index "$HUMAN_REF")
            mouse_idx=$(ensure_bwamem2_index "$MOUSE_REF")
            ;;
    esac

    # ── Python streaming engine ───────────────────────────────────────────────
    banner "Initialising k-mer streaming engine"
    create_stream_py
    info "Streaming script: $STREAM_PY"

    # ── Define output paths ────────────────────────────────────────────────────
    local human_bed="${OUT_DIR}/human_ambiguous_zones.bed"
    local mouse_bed="${OUT_DIR}/mouse_ambiguous_zones.bed"
    local human_vcf="${OUT_DIR}/human_on_mouse_variants.vcf"
    local mouse_vcf="${OUT_DIR}/mouse_on_human_variants.vcf"

    # ── Direction 1: Human → Mouse ────────────────────────────────────────────
    run_direction \
        "human" "mouse" \
        "$HUMAN_REF" "$MOUSE_REF" "$mouse_idx" \
        "$human_bed" "$human_vcf"

    # ── Direction 2: Mouse → Human ────────────────────────────────────────────
    run_direction \
        "mouse" "human" \
        "$MOUSE_REF" "$HUMAN_REF" "$human_idx" \
        "$mouse_bed" "$mouse_vcf"

    # ── Summary ───────────────────────────────────────────────────────────────
    banner "Pipeline complete — Summary"
    printf '\n'
    printf '  %-55s  %s\n' "Human ambiguous zones (BED):"    "$human_bed"
    printf '  %-55s  %s\n' "Mouse ambiguous zones (BED):"    "$mouse_bed"
    printf '  %-55s  %s\n' "Human→Mouse variants (VCF):"     "$human_vcf"
    printf '  %-55s  %s\n' "Mouse→Human variants (VCF):"     "$mouse_vcf"
    printf '\n'

    # File sizes
    for f in "$human_bed" "$mouse_bed" "$human_vcf" "$mouse_vcf"; do
        [[ -f "$f" ]] && info "  $(du -sh "$f" | cut -f1)  $f"
    done
    printf '\n'

    info "Integration with xenofilters:"
    info "  xenofilters human.bam mouse.bam \\"
    info "    --matching-algorithm hashlookup \\"
    info "    --ambiguous-regions  $human_bed  $mouse_bed \\"
    info "    --diagnostic-variants $human_vcf $mouse_vcf"
    printf '\n'
    done_ "All done."
}

main "$@"
