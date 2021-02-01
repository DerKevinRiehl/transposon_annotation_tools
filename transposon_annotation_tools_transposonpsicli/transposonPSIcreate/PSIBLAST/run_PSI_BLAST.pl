#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 refseq chkpFile database nuc|prot maxEvalue=1\n\n";

my $refseq = $ARGV[0] or die $usage;
my $chkp = $ARGV[1] or die $usage;
my $database = $ARGV[2] or die $usage;
my $type = $ARGV[3] or die $usage;
my $max_Evalue = $ARGV[4] || 1;

unless ($type eq 'nuc' || $type eq 'prot') { die $usage; }

main: {
	
	my $cmd;
	
	if ($type eq "nuc") {
		# do psitblastn w/ blastall
		$cmd = "blastall -i $refseq -d $database -p psitblastn -R $chkp -F F -M BLOSUM45 -t -1 -e $max_Evalue -v 10000 -b 10000 ";
	} else {
		# do psiblast w/ blastpgp
		$cmd = "blastpgp -i $refseq -d $database -R $chkp -j 1 -e $max_Evalue -v 10000 -b 10000 ";
	}
	
	my $ret = system($cmd);
	if ($ret) {
		die "Error, cmd $cmd died with ret $ret";
	}

	exit(0);
}
