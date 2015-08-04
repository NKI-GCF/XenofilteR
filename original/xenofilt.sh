#!/bin/bash

MM10_REF=$1
HG38_REF=$2

[ -z "$MM10_REF" ] && {
cat << EOF

Usage: ./xenofilt.sh xeno_reference.fasta host_reference.fasta

    xeno is foreign tissue, e.g. human, hg38 reference
    host is tissue, not of interest, e.g. mouse, mm10 reference

    all reference files need to be indexed with 'bwa index' if you want to use bwa.
    perl is required

EOF
exit 1
}


# install bwa and samtools (once), e.g.
[ -f bwa-0.7.12/bwa ] && alias bwa=bwa-0.7.12/bwa
[ -f samtools-0.1.19/samtools ] && alias samtools=samtools-0.1.19/samtools
donloaded=
if [ -z "$(which bwa)" ]; then
  downloaded=1
  wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2
fi
if [ -z "$(which samtools)" ]; then
  downloaded=1
  wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2
fi

if [ -n "$downloaded" ]; then
  ls -1 *.tar.bz2 | while read f; do
    tar xjf $f
    cd ${f%.tar.bz2}
    make
    cd -
  done
  [ -f bwa-0.7.12/bwa ] && alias bwa=bwa-0.7.12/bwa
  [ -f samtools-0.1.19/samtools ] && alias samtools=samtools-0.1.19/samtools
fi

# or if you want to use tophat/bowtie use that similarly.

# align and make sure the human (hg38) is sorted. make sure resulting file is called
# ${basename}-hg38.bam
# see http://bio-bwa.sourceforge.net/bwa.shtml

#your gzipped fastqs, if named differently change this and BN assignment below
for FQ1 in *_R1_001.fastq.gz; do

  BN=${FQ1%_R1_001.fastq.gz}

  # use 1 gig for sorting
  SAMMEM=1000000000
  # use 4 cpus
  NCPU=4

  #your gzipped fastqs
  FQ1=$(ls -1 ${BN}_R1_001.fastq.gz)

  #read 2 if you have it, for paired-end sequencing
  FQ2=$(ls -1 ${BN}_R2_001.fastq.gz 2>/dev/null)

  #select xenofilt script accordingly
  [ -z "$FQ2" ] && xenofilt_script=xenofilt_SE.pl || xenofilt_script=xenofilt_PE.pl

  HG38_BAM=${BN}-hg38
  HG38_SBAM=${BN}-hg38.bam

  # indexing is not strictly necessary
  bwa mem -M -t $NCPU -R "@RG\tID:$(< /dev/urandom tr -dc _A-Z-a-z-0-9 | head -c6)\tCN:NKI\tPL:ILLUMINA\tSM:${BN}\tLB:${BN}" $HG38_REF $FQ1 $FQ2 |
  samtools view -Sbu - |
  samtools sort -m $SAMMEM - $HG38_BAM
  samtools index $HG38_BAM&

  #indexing nor sort step is strictly necessary for mouse
  MM10_BAM=${BN}-mm10
  MM10_SBAM=${BN}-mm10.bam
  bwa mem -M -t $NCPU -R "@RG\tID:$(< /dev/urandom tr -dc _A-Z-a-z-0-9 | head -c6)\tCN:NKI\tPL:ILLUMINA\tSM:${BN}\tLB:${BN}" $MM10_REF $FQ1 $FQ2 |
  samtools view -Sbu - |
  samtools sort -m $SAMMEM - $MM10_BAM
  samtools index $MM10_SBAM&

  HG38_FILT_BAM=${BN}-hg38-filt
  HG38_FILT_SBAM=${BN}-hg38-filt.bam

  perl ./$xenofilt_script $MM10_SBAM $HG38_SBAM |samtools view -S -b - | samtools sort -m 2G - $HG38_FILT_BAM
  samtools index $HG38_FILT_SBAM&

done


