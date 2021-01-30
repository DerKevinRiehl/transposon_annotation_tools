####this script want to get a good sequence-seed which is used for genome scan in annotation step, here we get seeds from MSA###
#!usr/bin/perl
use strict;
use Bio::SimpleAlign;
use Bio::AlignIO;

if(@ARGV != 1){
	warn"please enter: sequence file \n";
	exit 0;
}

my $file=$ARGV[0];####MSA file
#my $cutoff=$ARGV[1];####consensus of each column cutoff

open IN,$file or die "$!\n";
my %seq=();
my $line;
my $u;
my $i=0;
while(defined($line=<IN>)){
	chomp $line;
	if($line=~/>/){
		$line=~s/>//;
		my @a=split(/\s+/,$line);
		$u=$a[0];
#		($u)=/(.*)_/;
		$seq{$u}="";		
	}else{
		$seq{$u}.=$line;
	}	
}
close IN;

my %position=();
open IN,$file.".order" or die "$!\n";
while(<IN>){
	chomp $_;
	if($_=~/>/){
		$_=~s/>//;
		$u=$_;
	}else{
		$position{$u}=$_;
	}
}
close IN;

my $outfile=$file.".better.fa";
my $matchNum=0;
###7SL tRNA
my $a_box="[GA][CGAT]TGG|TGGCTCACGCC|T[AG]G[CT]\\w{2}A\\w{3}G";
my $spacer_1="\\w{25,70}"; 
my $b_box="GTTC[AG]A|GTTCGAGAC|G[AT]TC[AG]A\\w{2}C";
##5S		
my $A_box="[ATC][AG]G[CT][CT]AAGC";
my $Spacer_1="\\w{20,50}";
my $B_box="[AG]TGG[AG][ATG]GAC";
open fout,">$outfile" or die "$!\n";
foreach my $i (keys%seq){
	my @a=split(/\t/,$position{$i});
	my $s=substr($seq{$i},$a[0],length($seq{$i})-$a[0]-$a[1]);
	my $S=reverse($s);
	$S=~tr/TACG/ATGC/;
	if($s=~/($a_box)($spacer_1)($b_box)/ or $S=~/($a_box)($spacer_1)($b_box)/){
		$matchNum++;	
	}elsif($s=~/($A_box)($Spacer_1)($B_box)/ or $S=~/($A_box)($Spacer_1)($B_box)/){
		$matchNum++;	
	}
	print fout ">$i\n$s\n";
}
close fout;
###whether exists known pattern####
if($matchNum == 0){
	exit 0;	
}

my $output=$file.".better.out";
system "/usr/bin/muscle -in $outfile -out $output -maxiters 1 -diags -quiet";
open IN,$output or die "$!\n";
while(defined($line=<IN>)){
	chomp $line;
	if($line=~/>/){
		$line=~s/>//;
		$_=$line;
		($u)=/(.*)_/;
		$seq{$u}="";		
	}else{
		$seq{$u}.=$line;
	}	
}
close IN;


my $cons=consensus($output);
#open fout,">$file.next.fa" or die "$!\n";
#print fout"$cons\n";
#close fout;
#$cons=`perl discard-poly-A.pl $file.extendseqbysolar.next.fa`;
my @s=split(/,/,$cons);
$file=~s/\.extendseq//;
open fout,">$file.for_annotation.fa" or die "$!\n";
for(my $i=0;$i<@s;$i++){
	my @a=split(/\//,$file);
	my $id=$a[-1]."_$i";
	print fout">$id\n$s[$i]\n";
}
close fout;
system "rm $output $outfile";

sub consensus{
        my $fl=$_[0];
        my $str = Bio::AlignIO->new(-file => $fl);
        my $aln = $str->next_aln();
	my $nice=0;
	my $nseq=$aln->no_sequences;
	if($nseq > 6){
		$nice=50;
	}else{
		$nice=30;
	}
	my $one=$aln->consensus_string(80)."\n";
	my $two=$aln->consensus_string($nice)."\n";############here one question left behind#######
	$one=~s/\?//g;
	$two=~s/\?//g;
	chomp $one;
	chomp $two;
	my @Seqs=();
	push @Seqs,$one;
	if(abs(length($one)-length($two)) > 0.1*length($one)){
		push @Seqs,$two;
	}
	my $CONS="";
	foreach my $cons (@Seqs){
		###roughly discard poly-A#####
		if($cons=~/^\w{0,2}[A|T]{5,}/){
			$cons=reverse($cons);
			$cons=~tr/ATCG/TAGC/;	
		}
		while($cons=~/[A|T]{5,}\w{0,2}$/ or $cons=~/^\w{0,2}[A|T]{5,}/){
			$cons=~s/[A|T]{5,}\w{0,2}$//;
			$cons=~s/^\w{0,2}[A|T]{5,}//;
		}
		$CONS.="$cons,";
	}
        return($CONS);
}
