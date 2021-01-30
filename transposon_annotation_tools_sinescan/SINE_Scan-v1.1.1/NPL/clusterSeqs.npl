#!usr/bin/perl
use strict;
use File::Basename;

open in,$ARGV[0] or die "$!\n";
my $path=dirname($ARGV[0]);
my $sinefile=$ARGV[1];
my %sines=();
my %pos=();
if(-e $sinefile){
	system "rm $sinefile";
}
while(<in>){
	chomp $_;
	my $u=$_;
	if(-e "$path/$u/$u.sine.fa"){
		open IN,"$path/$u/$u.sine.fa" or die "$!\n";
		my $line;
		my $flag=0;
		my $Number=$u;
		while(defined($line=<IN>)){
			chomp $line;
			if($line=~/>/){
				my @a=split(/\|/,$line);
				$u=$line;
				my @b=split(/-/,$a[1]);
				my $w="$a[0],$a[1]";
				if(!exists $pos{$w}){
					$sines{$u}="";
					$flag=1;
					$pos{$w}=1;
				}else{
					$flag=0;
				}
			}else{
				if($flag == 1){
					$sines{$u}.=$line;
				}
			}
		}
		close IN;
	}
}
close in;
undef%pos;

open out,'>',$sinefile or die "$!\n";
foreach my $i (keys%sines){
	print out"$i\n$sines{$i}\n";
}
close out;
