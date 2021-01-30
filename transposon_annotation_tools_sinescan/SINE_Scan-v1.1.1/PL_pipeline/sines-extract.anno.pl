#!usr/bin/perl
use strict;

if(@ARGV < 2){
	warn "please enter : TE_list genome_file\n";
	exit 0;
} 

open IN,$ARGV[0] or die "$!\n";
my $line;
my %name=();
my @Name=();
while(defined($line=<IN>)){
	chomp $line;
	my @s=split(/\t/,$line);
	my $flag='+';
	if($s[8]-$s[9] < 0){
		$flag='-';
	}
	my $u="$s[0]:$flag";	
	my @a=split(/:/,$s[0]);
	$name{$u}=$a[0];
	push @Name,$u;
}
close IN;

open IN,$ARGV[1] or die "$!\n";
my %chr=();
my $r="";
while(defined($line=<IN>)){
	chomp $line;
	if($line=~/>/){
		$line=~s/>//;
		my @s=split(/\s+/,$line);
		if(exists $chr{$r}){
			my %happen=();
			foreach my $i (@Name){
				my @a=split(/:/,$i);
				my @b=split(/-/,$a[1]);
				my $w="$a[0]:$a[1]";
				if($name{$i} eq $r && !exists $happen{$w}){
					my $seq=substr($chr{$r},$b[0],$b[1]-$b[0]+1);
					$happen{$w}=1;
					print ">$i\n$seq\n";
				}
			}
			undef%happen;
		}
		$chr{$r}="";
		$chr{$s[0]}="";
		$r=$s[0];
	}else{
		$chr{$r}.=$line;
	}
}
close IN;

my %happen=();
foreach my $i (@Name){
	my @a=split(/:/,$i);
	my @b=split(/-/,$a[1]);
	my $w="$a[0]:$a[1]";
	if($name{$i} eq $r && !exists $happen{$w}){
		my $seq=substr($chr{$r},$b[0],$b[1]-$b[0]+1);
		$happen{$w}=1;
		print ">$i\n$seq\n";
	}
}
undef%happen;
undef%name;
undef%chr;
@Name=();

sub rev{
        my $x=reverse($_[0]);
        $x=~tr/ATCG/TAGC/;
        $x=~tr/atcg/tagc/;
        return($x);
}
