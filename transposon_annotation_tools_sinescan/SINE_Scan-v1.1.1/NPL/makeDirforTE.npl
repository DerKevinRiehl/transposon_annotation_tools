#!usr/bin/perl
use strict;

open IN,$ARGV[0] or die "$!\n";
my %name=();
my $workdir=$ARGV[1];
my $flag=0;
my $label=-1;
my $u="";
my $Flag=0;
while(<IN>){
	if($_=~/>/){
		if(!exists $name{$_}){
			$label++;
			$name{$_}=1;
			$u=$_;
			$flag=1;
			$Flag=1;
			if(-e "$workdir/$label"){
		        	unlink glob "$workdir/$label/* $workdir/$label/.*";
			}else{
        			mkdir "$workdir/$label", 0755 or warn "cannot create $workdir/$label directory:$!";
			}
		}else{
			$flag=0;
		}	
	}else{
		if($flag == 1){
			open fout,'>>',"$workdir/$label/$label.sine.fa" or die "$!\n";
			if($Flag == 1){
				print fout"$u";
				$Flag=0;
			}
			print fout "$_";
			close fout;
		}
	}
}
close IN;
