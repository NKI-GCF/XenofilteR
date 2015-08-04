#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: xenofilt.pl
#
#        USAGE: ./xenofilt.pl  
#
#  DESCRIPTION: filter xenograft reads
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: (c) Roel 2014 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 01/10/2014 12:25:31 PM
#     REVISION: ---
#===============================================================================
#
#Loop over primary host (mouse) alignments. Store separately for first and second reads
#in a pair the number of mismatches plus the number of clipped NTs multiplied by 1000 / read_quality.
#If the read is not mapped instead assign the exceeding balance of readlength * 1000.
#The readname + first/second [read in pair] are used as key.
#
#Subtract for the xenograft (human) alignments mismatches similarly for both reads.
#When the balance is complete for both primary reads, we keep these reads if:
#* both balances were positive, i.e. the (somewhat weighed) mismatches on the host
#(mouse) were greater than on the xenograft (human). This is counted as
#`matches_xeno_better'.
#* the sum for both reads is still positive. This is `read_comb_matches_xeno_better'
#
#If primary alignments are printed, non-primary - stored for xeno in memory, some
#pending - are printed as well.
#
#We loose the reads where this is not the case, or the result is `ambiguous'


use strict;
use warnings;
use utf8;

my $hostf = shift or die;
my $samplef = shift or die;
die unless -f $hostf;
die unless -f $samplef;

#die if -f $filtf;
open (HOST, "samtools view '$hostf'|") or die "cannot samtools view $hostf:$!";
open (SAMP, "samtools view -h '$samplef'|") or die "cannot samtools view $samplef:$!";

my %h;
sub print_stats {
	my $stats = shift;
	warn "count\tstat\n";
	foreach my $t (sort {$a cmp $b} keys %$stats) {
		warn join("\t", $stats->{$t}, $t)."\n";
	}
}
# create tag for host parsing
# store in hash entries which map perfectly to the host (e.g. mouse);

# Iterate over de host mappings. Skip secondary alignments en unmapped reads tenzij een van beide reads wel mapt.
# Voor beide reads sommeer het aantal mismatches en sla deze met de mapping quality op in een hash. De readname als key.
HOSTLINE: while (<HOST>) {
	chomp;
	my @e = split /\t/;
	next if $e[1] & 256; # ignore secondary alignments in host.
	next if ($e[1] & 4) and $e[2] eq '*'; # unmapped; associated if $e[2] ne '*';
	$e[0] =~ s/^HWI-[^:]++:[0-9]++:[^:]++://;
	my $rd = $e[1] & 128 ? 1 : 0;

	# sum mismatches (NM plus clips) and divide by quality
	if (($e[4] == 0) || ($e[5] eq '*')) { # entire fail is 1000 * readlength
		$h{$e[0]}->[$rd] = length($e[9]) * 1000;
	} else {
		my $i = $#e;
		$i-- while (($i > 10) && ($e[$i] !~ /^NM:i:([0-9]+)$/));
		if ($i != 10) {
			$i = $1;
			while ($e[5] =~ s/([0-9]++)S//o) { $i += $1; }
			$h{$e[0]}->[$rd] = $i * 1000 / $e[4]; # matches * 1000 / read_quality
		} else { # entire fail
			$h{$e[0]}->[$rd] = length($e[9]) * 1000;
		}
	}
}
close HOST;
warn "finished reading host...\n";

my (%pair, %stats, %alt, $L);
while ($L = <SAMP>) {
	last if $L !~ /^@/;
	print $L; #no newline needed
}
while ($L) {
	chomp $L;
	my $e = [split /\t/, $L];
	die unless $e->[0] =~ m/^HWI-[^:]++:[0-9]++:[^:]++:(.*)/;
	if (exists $h{$1}) {
		my ($id, $rd) = ($1, $e->[1] & 128 ? 1 : 0);
		if ($e->[1] & 256) { # non primary, unused but printed along if others evaluated matchin xeno
			if ((not exists $pair{$id}) || (ref($pair{$id}) eq 'ARRAY')) {
				push @{$alt{$id}}, $e;
			} elsif ($pair{$id} == 1) {
				print "$L\n";
			}
			next;
		}

		# subtract HG19 mismatches (NM plus clips) divided by quality
		if (($e->[4] == 0) || ($e->[5] eq '*')) {
			$h{$id}->[$rd] -= length($e->[9]) * 1000;
		} else {
			my $i = $#$e;
			$i-- while (($i > 10) && ($e->[$i] !~ /^NM:i:([0-9]+)$/));
			if ($i != 10) {
				$i = $1;
				my $cig = $e->[5];
				while ($cig =~ s/([0-9]++)S//o) { $i += $1; }
				$h{$id}->[$rd] -= $i * 1000 / $e->[4];
			} else {
				$h{$id}->[$rd] -= length($e->[9]) * 1000;
			}
		}
		if (exists $pair{$id}) {
			if (ref($pair{$id}) eq 'ARRAY') {
				my $balance = $h{$id};
				if (($balance->[0] > 0) && ($balance->[1] > 0)) { # both say it's more likely human
					$stats{matches_xeno_better}++;
					if (exists $alt{$id}) {
						print join("\t", @{$_})."\n" foreach @{$alt{$id}};
						delete $alt{$id};
					}
					print join("\t", @{$pair{$id}})."\n$L\n";
					$pair{$id} = 1;
				} elsif ($balance->[0] < 0 and $balance->[1] < 0) {
					$stats{matches_host_better}++;
					$pair{$id} = 0;
				} elsif ($balance->[0] + $balance->[1] > 0) {
					$stats{read_comb_matches_xeno_better}++;
					if (exists $alt{$id}) {
						print join("\t", @{$_})."\n" foreach @{$alt{$id}};
						delete $alt{$id};
					}
					print join("\t", @{$pair{$id}})."\n$L\n";
					$pair{$id} = 1;
				} elsif ($balance->[0] != $balance->[1]) {
					$pair{$id} = 0;
					$stats{read_comb_matches_host_better}++;
				} else {
					$pair{$id} = 0;
					$stats{ambiguous}++;
				}
			} elsif ($pair{$id} == 1) {
				print "$L\n";
			}
		} else {
			$pair{$id} = $e;
		}
	} else { # not mapped on host
		print "$L\n";
		$stats{xeno_mapct}++;
	}
} continue {
	$L = <SAMP>;
}
close SAMP;
print_stats(\%stats);

