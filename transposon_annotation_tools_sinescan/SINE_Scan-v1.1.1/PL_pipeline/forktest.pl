#!usr/bin/perl
use strict;
use Parallel::ForkManager;

my @a=('./PL_pipeline/tRNA-7SL-SINE-FINDER.py -T chunkwise -f fasta Megabat.sineKnowns.fa','./PL_pipeline/5S-SINE-FINDER.py -T chunkwise -f fasta -s 0 Megabat.sineKnowns.fa');
my $pm=Parallel::ForkManager->new(5);

foreach my $i (@a){
	$pm->start and next;
	system "python $i";
	$pm->finish;
}
$pm->wait_all_children;
print "11,11\n";
