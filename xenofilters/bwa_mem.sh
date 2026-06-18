#!/bin/bash
set -e

graft="$1"; shift;
host="$1"; shift;
bam="$1"; shift;
fq1="$1"; shift;

[ -n "$1" -a "${1:0:1}" != '-' ] && {
  fq2="$1"; shift
}

DEFAULT_BWA_ARGS="-M -t 4"

die() {
cat << EOF 1>&2
  Usage:

bwa=/path/to/bwa \\
samtools=/path/to/samtools \\
$0 /path/to/graft_reference.fa /path/to/host_reference.fa \$bam \$fq1 [\$fq2] [bwa arguments.. (default $DEFAULT_BWA_ARGS)]

Error: $1

EOF
  exit 1
}

BWA_ARGS="$@"
[ -z "$BWA_ARGS" ] && BWA_ARGS="$DEFAULT_BWA_ARGS"

[ -e "$graft" ] || die "'$graft' (graft reference):No such file"
[ -e "${graft}.pac" ] || die "$graft (graft reference):Not bwa indexed (hint: consider running 'bwa index $graft' first)"

[ -e "$host" ] || die "'$host' (host reference):No such file"
[ -e "${host}.pac" ] || die "$host (host reference):Not bwa indexed (hint: consider running 'bwa index $host' first)"

[ -n "$bam" ] || die "'$bam' (bam): requires filename"
[ ! -e "$bam" ] || die "$bam (bam):Already exists"

[ -e "$fq1" ] || die "'$fq1' (fastq1):No such file"
[ -z "$fq2" -o -e "$fq2" ] || die "'$fq2' (fastq2):No such file"

[ -n "$bwa" ] || bwa="$(which bwa)"
[ -x "$bwa" ] || die "$bwa (bwa): no executable found"

[ -n "$samtools" ] && samtools="$(which samtools)"
[ -x "$samtools" ] || die "$samtools (samtools): no executable found"

BN="${bam%.bam}"

rg="@RG\tID:$(mktemp -u | cut -d '.' -f 2)\tCN:$(whoami)\tPL:ILLUMINA\tSM:$BN\tLB:$BN"
target/release/xenofilter <($bwa mem $BWA_ARGS -R "$rg" $graft $fq1 $fq2) \
<($bwa mem $BWA_ARGS $host $fq1 $fq2) |
$samtools view -Sbu - | $samtools sort -m 1294967296 -@2 - -o $bam &&
$samtools index $bam

