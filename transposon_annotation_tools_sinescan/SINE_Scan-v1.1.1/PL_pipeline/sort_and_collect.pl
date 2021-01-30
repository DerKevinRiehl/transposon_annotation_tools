#!usr/bin/perl
use strict;

my %count=();
open in,$ARGV[0] or die "$!\n";
while(<in>){
	if($_=~/>/){
#		(my $r)=/\] (.*)/;
		my @a=split(/\|/,$_);
		my $r=$a[-2];
		#$r=lc($r);
		$count{$r}++;
	}
}
close in;

foreach my $i (sort{$count{$a}<=>$count{$b}} keys%count){
	print "$i\t$count{$i}\n";
}
