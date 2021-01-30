#!usr/bin/perl
use strict;

open in,$ARGV[0] or die "$!\n";
open out,'>',$ARGV[0].".tmp" or die "$!\n";
my $u;
###7SL tRNA
my $a_box="[GA][CGA]TGG|TGGCTCACGCC|T[AG]G[CT]\\w{2}A\\w{3}G";
my $spacer_1="\\w{25,70}"; 
my $b_box="GTTC[AG]A|GTTCGAGAC|G[AT]TC[AG]A\\w{2}C";
##5S		
my $A_box="[ATC][AG]G[CT][CT]AAGC";
my $Spacer_1="\\w{20,50}";
my $B_box="[AG]TGG[AG][ATG]GAC";
while(<in>){
	if($_=~/>/){
		$u=$_;
	}else{
		my $s=substr($_,0,110);
		if($s=~/($a_box)($spacer_1)($b_box)/){
			print out "$u$_";
			next;
		}
		if($s=~/($A_box)($Spacer_1)($B_box)/){
			print out "$u$_";
#			print "5S RNA: $_";
		}
	}
}
close in;
close out;
