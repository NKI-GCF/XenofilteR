#!/usr/bin/env bash
# shellcheck disable=SC2034
# Pass 1: standard run. Pass 2: BQSR-corrected ambiguous re-scoring.
# Sorting before ApplyBQSR: NOT required (GATK processes sequentially).
set -euo pipefail

GRAFT_BAM=$1;
HOST_BAM=$2
GRAFT_REF=$3;
HOST_REF=$4
KNOWN_SITES=$5
OUTDIR=${6:-./results}; THREADS=${7:-16}

mkdir -p "$OUTDIR"

# -- Pass 1 --------------------------------------------------------------------
xenofilters "$GRAFT_BAM" "$HOST_BAM" \
  --output "$OUTDIR/pass1_graft.bam,$OUTDIR/pass1_host.bam" \
  --ambiguous-output "$OUTDIR/ambiguous_graft.bam,$OUTDIR/ambiguous_host.bam" \
  --threads "$THREADS" --score-threads "$THREADS"

# -- BSQR & Haplotype calling on pass-1 output --------------------------------
for STREAM in graft host; do
    REF_VAR="${STREAM^^}_REF"; REF="${!REF_VAR}"
    gatk BaseRecalibrator \
      -I "$OUTDIR/pass1_${STREAM}.bam" -R "$REF" \
      --known-sites "$KNOWN_SITES" \
      -O "$OUTDIR/${STREAM}.recal.table" --verbosity ERROR &
done
wait

# Merge tables (same run; NEEDS VERIFICATION: GatherBQSRReports cross-reference).
gatk GatherBQSRReports \
  -I "$OUTDIR/graft.recal.table" \
  -I "$OUTDIR/host.recal.table" \
  -O "$OUTDIR/combined.recal.table"

# Apply to ambiguous and pass1 BAM (no sort required before ApplyBQSR).
for OUTPUT in pass1 ambiguous; do
  for STREAM in graft host; do
    REF_VAR="${STREAM^^}_REF"; REF="${!REF_VAR}"
    gatk ApplyBQSR \
      -I "$OUTDIR/${OUTPUT}_${STREAM}.bam" -R "$REF" \
      --bqsr-recal-file "$OUTDIR/combined.recal.table" \
      -O "$OUTDIR/${OUTPUT}_recal_${STREAM}.bam"
  done
done
wait

# Winners: use as --sample-variants in pass 2 to improve variant rescue.
for STREAM in graft host; do
    REF_VAR="${STREAM^^}_REF"; REF="${!REF_VAR}"
    gatk HaplotypeCaller \
      -I "$OUTDIR/pass1_${STREAM}.bam" -R "$REF" \
      -O "$OUTDIR/pass1_${STREAM}.vcf.gz" \
      --verbosity ERROR &
done
wait

# -- Pass 2 --------------------------------------------------------------------
# Threshold auto → 0 (pass 2 detected from _xenoambig RG in input headers).
xenofilters \
  "$OUTDIR/ambiguous_recal_graft.bam" \
  "$OUTDIR/ambiguous_recal_host.bam" \
  --output "$OUTDIR/pass2_graft.bam,$OUTDIR/pass2_host.bam" \
  --sample-variants "0:$OUTDIR/pass1_graft.vcf.gz,1:$OUTDIR/pass1_host.vcf.gz" \
  --threads "$THREADS" --score-threads "$THREADS"

# -- Final merge ---------------------------------------------------------------
samtools merge -@ "$THREADS" -f \
  "$OUTDIR/final_graft.bam" \
  "$OUTDIR/pass1_graft.bam" "$OUTDIR/pass2_graft.bam"

samtools merge -@ "$THREADS" -f \
  "$OUTDIR/final_host.bam" \
  "$OUTDIR/pass1_host.bam" "$OUTDIR/pass2_host.bam"

samtools index "$OUTDIR/final_graft.bam"
samtools index "$OUTDIR/final_host.bam"

echo "Done. Final outputs:"
echo "  $OUTDIR/final_graft.bam"
echo "  $OUTDIR/final_host.bam"
echo "  $OUTDIR/pass2_discarded.bam  (ambiguous reads unresolved after BQSR)"
