#!/usr/bin/env perl
# (c) Roel 2016

use strict;
use warnings;

# Run this to get slices of reads for mappings after xenofilter
#
# filename tags to_host means e.g. mapped to mm10, to_graft is e.g. mapped to hg38.
# filtered/nonfiltered is according to presence in the xf/ set.
# graftunmapped means it is present in the unmapped.bam of the hg38 alignment.
#


# usage:
# ./xf_slices.pl  $samp/accepted_htis.bam  $samp/unmapped.bam  xf/Filtered_bams/${samp}_Filtered.bam \
#    mm10/${samp}/accepted_htis.bam  slice/${samp}


my $graftBam = shift or die;
my $graftUnmapped = shift or die;
my $graftBamFilt =  shift or die;
my $hostBam = shift or die;
my $basename = shift or die;

my %h;
open(GFT, "samtools view -h $graftBam |") or die;
open(GFTUNM, "samtools view -h $graftUnmapped |") or die;
open(XF, "samtools view -h $graftBamFilt |") or die;
open(HST, "samtools view -h $hostBam |") or die;

open(FILTGFT, "| samtools view -S -b - > ${basename}_filt_from_graft.bam") or die;

open(NONFILTINGFT2HST, "| samtools view -S -b - > ${basename}_nonfiltered_to_host.bam") or die;
open(FILTINGFT2HST, "| samtools view -S -b - > ${basename}_filtered_to_host.bam") or die;
open(UNMINGFT, "| samtools view -S -b - > ${basename}_graftunmapped_to_host.bam") or die;

while (my $x = <GFTUNM>)
{
	$h{substr($x, 0, index($x, "\t"))} = -1 if substr($x, 0, 1) ne '@';
}
close GFTUNM;
warn "done with $graftUnmapped\n";

while (my $x = <XF>)
{
	if (substr($x, 0, 1) ne '@') {
		my $k = substr($x, 0, index($x, "\t"));
		while (my $g = <GFT>) {
			# store only the reads that were filtered. i.e. e.g. reads that were presumed mouse.
			my $r = substr($g, 0, index($g, "\t"));
			last if $r eq $k;
			if (substr($r, 0, 1) ne '@') {
				print FILTGFT $g;
				$h{$r} = 1;
			}
		}
	} else {
		print FILTGFT $x;
	}
}
close FILTGFT;
warn "done writing ${basename}_filt_from_graft.bam\n";

while (my $x = <HST>)
{
	if (substr($x, 0, 1) ne '@') {
		my $k = substr($x, 0, index($x, "\t"));
		if (exists $h{$k}) { # filtered reads in graft were presumed mouse. Here's the host mapping
			if ($h{$k} == -1) {
				print UNMINGFT $x;
			} else {
				print FILTINGFT2HST $x;
			}
		} else {
			# not presumed mouse (preserved in graft). Here's the host mapping
			print NONFILTINGFT2HST $x;
		}
	} else {
		print UNMINGFT $x;
		print FILTINGFT2HST $x;
		print NONFILTINGFT2HST $x;
	}
}
close FILTINGFT2HST;
close NONFILTINGFT2HST;
close UNMINGFT;
close HST;

