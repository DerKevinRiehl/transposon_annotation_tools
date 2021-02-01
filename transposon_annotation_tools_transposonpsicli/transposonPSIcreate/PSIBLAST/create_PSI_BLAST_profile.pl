#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 refprot protDB\n\n";

my $refprot = $ARGV[0] or die $usage;
my $protDB = $ARGV[1] or die $usage;

main: {

	my $cmd = "blastpgp -i $refprot -d $protDB -j 2 -C $refprot.chkp";
	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	
	exit(0);
}




	
