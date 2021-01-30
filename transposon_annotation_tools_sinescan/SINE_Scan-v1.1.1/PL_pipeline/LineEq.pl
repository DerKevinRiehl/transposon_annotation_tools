#!usr/bin/perl
use strict;

my @a=();
open in,$ARGV[0] or die "$!\n";
while(<in>){
	chomp $_;
	push @a,$_;
}
close in;

my @b=();
open in,$ARGV[1] or die "$!\n";
while(<in>){
	chomp $_;
	push @b,$_;
}
close in;

for(my $i=0;$i<@a;$i++){
	if($a[$i] ne $b[$i]){
		print "$i,$a[$i],$b[$i]\n";	
		exit 0;
	}
}
