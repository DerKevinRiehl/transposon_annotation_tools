#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\nusage: $0 tblastn_chain_output\n\n";

my $tblastn_chain_output = $ARGV[0] or die $usage;


main: {
	
	my %genome_to_chains = &parse_chained_hits($tblastn_chain_output);

	foreach my $genome (sort keys %genome_to_chains) {
		
		my @chains = sort {$a->{genome_lend}<=>$b->{genome_lend}} @{$genome_to_chains{$genome}};
		
		for (my $i = 1; $i <= $#chains; $i++) {

			my $struct_i = $chains[$i];
			my ($genome_lend_i, $genome_rend_i) = ($struct_i->{genome_lend}, $struct_i->{genome_rend});

			for (my $j = $i - 1; $j >= 0; $j--) {

				my $struct_j = $chains[$j];
				my ($genome_lend_j, $genome_rend_j) = ($struct_j->{genome_lend}, $struct_j->{genome_rend});

				# make sure they don't overlap
				
				if ($genome_lend_i < $genome_rend_j && $genome_rend_i > $genome_lend_j) { next; } # overlap detected
				
				if ($struct_i->{score} + $struct_j->{sum_path_score} > $struct_i->{sum_path_score}) {
					$struct_i->{prev} = $struct_j;
					$struct_i->{sum_path_score} = $struct_i->{score} + $struct_j->{sum_path_score};
				}
			}
		}

		@chains = sort {$a->{sum_path_score}<=>$b->{sum_path_score}} @chains;

		my $best_chain = $chains[$#chains];
		
		my @report_chains;
		while (defined $best_chain) {
			push (@report_chains, $best_chain);
			$best_chain = $best_chain->{prev};
		}
		
		foreach my $chain (reverse @report_chains) { # order left to right
			print $chain->{chain_txt} . "\n";
		}	
	}
		
	exit(0);
}


####
sub parse_chained_hits {
	my ($chain_file) = @_;
	
	my %genome_to_chain_structs;
	
	my $curr_chain_struct = undef;

	open (my $fh, $chain_file) or die "Error, cannot open file $chain_file";
	while (<$fh>) {
		unless (/\w/) { next; }
		my $line = $_;
		chomp;
		if (/^\#/) {
			
			my ($chaintxt, $db_acc, $db_span, $genome_acc, $genome_span, $genome_orient, $score) = split (/\t/);
			
			my ($genome_lend, $genome_rend) = split (/-/, $genome_span);
			

			$curr_chain_struct = { genome_lend => $genome_lend,
								   genome_rend => $genome_rend,
								   score => $score,
								   sum_path_score => $score,
								   chain_txt => "",
							   };
			
			push (@{$genome_to_chain_structs{$genome_acc}}, $curr_chain_struct);
			
		}

		$curr_chain_struct->{chain_txt} .= $line;
	}
	close $fh;

	return (%genome_to_chain_structs);
}
