#!/usr/bin/env perl

use strict;
use warnings;

use Fasta_reader;
use Nuc_translator;

my $usage = "usage: $0 TPSI.chains genome.fa\n\n";

my $TPSI_chains = $ARGV[0] or die $usage;
my $genome_fa = $ARGV[1] or die $usage;

my $MAX_ELE_LENGTH = 20000;

main: {

	my $fasta_reader = new Fasta_reader($genome_fa);

	my %genome = $fasta_reader->retrieve_all_seqs_hash();

	my $counter = 0;

	open (my $fh, $TPSI_chains) or die "Error, cannot open file $TPSI_chains";
	while (<$fh>) {
		chomp;
		if ( /^\#Chain/) {
			my $line = $_;
			
			my @x = split(/\t/);
			my $type = $x[1];
			my $TE_range = $x[2];
			my $scaffold = $x[3];
			my $genome_range = $x[4];
			my $strand = $x[5];
			
			my ($range_lend, $range_rend) = sort {$a<=>$b} split(/-/, $genome_range);
			
			my $ele_len = $range_rend - $range_lend + 1;
			if ($ele_len > $MAX_ELE_LENGTH) {
				print STDERR "warning, ele too large: $ele_len, skipping $line\n";
				next;
			}
			
			

			my $seq = substr($genome{$scaffold}, $range_lend - 1, $ele_len);
			
			if ($strand eq '-') {
				$seq = &reverse_complement($seq);
			}

			$seq =~ s/(\S{60})/$1\n/g;

			$counter++;	

			chomp $seq;
			
			$line =~ s/^\#Chain\s+//;

			print ">ele_$counter $line\n$seq\n";
			
		}
	}


	exit(0);
}

