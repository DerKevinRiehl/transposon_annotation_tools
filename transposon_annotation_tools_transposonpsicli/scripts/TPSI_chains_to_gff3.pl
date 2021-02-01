#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\nusage; $0 btab_file\n\n";

my $btab_file = $ARGV[0] or die $usage;

open (my $fh, $btab_file) or die "Error, cannot open file $btab_file";

my $chain_counter = 0;

while (<$fh>) {
	chomp;
	unless (/\w/) { next; }
	if (/^\#Chain/) { 
		$chain_counter++;
		next; 
	}
	

	my @x = split (/\t/);
	my $hit_acc = $x[0];
	my $contig_acc = $x[5];

	my ($hit_end5, $hit_end3) = ($x[6], $x[7]);
	my ($contig_end5, $contig_end3) = ($x[8], $x[9]);
	
	my $contig_orient = ($contig_end5 < $contig_end3) ? '+' : '-';

	my ($contig_lend, $contig_rend) = sort {$a<=>$b} ($contig_end5, $contig_end3);

	my $blast_score = $x[12];
	my $evalue = $x[19];

	my $match_count = sprintf("chain%05d", $chain_counter);
	print join ("\t", $contig_acc, "TransposonPSI", "translated_nucleotide_match", 
				$contig_lend, $contig_rend, $blast_score, $contig_orient, ".", 
				"ID=$match_count; Target=$hit_acc; E=$evalue $hit_end5 $hit_end3") . "\n";
}

close $fh;


exit(0);

