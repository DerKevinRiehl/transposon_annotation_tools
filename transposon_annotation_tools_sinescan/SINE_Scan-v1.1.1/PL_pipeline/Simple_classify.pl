#!usr/bin/perl
use strict;
use Bio::AlignIO;
use Bio::Align::PairwiseStatistics;
use File::Basename;

if(@ARGV != 9){
	print "$0 <output-prefix> <genus_species_strain:care '-'> <cd-hit one> <cd-hit two> <blast-one> <rna-one> <rna-two> <outputDir> <family file>\n";
	print "Be sure TE database file and genome-mapped file have been formated for blast\n";
	exit(0);
}

my $script=dirname($0);
my $te="./SINEBase/SineDatabase.fasta";
my $prefix=$ARGV[0];
my $tRNA="./RNABase/RNAsbase.fasta";
my @A=split(/_/,$ARGV[1]);
#####cd-hit parameters#########
my $cd_one=$ARGV[2];##identity
my $cd_two=$ARGV[3];##overlap in length
#####blast filter parameter####
my $Blast_cut=$ARGV[4];
####rna align parameters########
my $rna_one=$ARGV[5];###identity
my $rna_two=$ARGV[6];##overlap in length
my $outputDir=$ARGV[7];###the directory of SINEs annotation files
my $Family=$ARGV[8];

my $n=@A;
my $speciesname=$ARGV[1];###to every users##########

system "/home/maohlzj/sine_te/cd-hist/cd-hit-est -i $Family  -o $Family.cluster -n 10 -c $cd_one -d 0 -r 1 -s $cd_two -aS $cd_two -aL $cd_two";

my $line;
my %cluster=();
my %represent=();
my %Represent=();
my %Num=();
open IN,$Family.".cluster.clstr" or die "$!\n";
my $u;
while(defined($line=<IN>)){
	chomp $line;
	if($line=~/^>/){
		$line=~s/>//;
		my @p=split(/\s+/,$line);
		$line=~s/\s/_/;
		$u=$line;
		$Num{$u}=$p[-1];
	}else{
		my @a=split(/\s/,$line);
		$a[2]=~s/>//;
		$a[2]=~s/\.{3}//;
		$cluster{$a[2]}=$u;
		if($line=~/\*/){
			$represent{$u}="$a[2],";
			$Represent{$u}=$a[2];
		}
	}
}
close IN;

my $out=$prefix.".tab";
system "/home/maohlzj/ncbi-blast-2.2.31+/bin/blastn  -task blastn -query $Family.cluster -db $te -max_target_seqs 100000 -evalue 1e-10 -dust no -outfmt 6 -out $out";

open IN,$te or die "$!\n";
my %length=();
my $u="";
while(defined($line=<IN>)){
	chomp $line;
	if($line=~/>/){
		$line=~s/>//;
		$u=$line;
		$length{$u}=0;	
	}else{
		$length{$u}+=length($line);
	}
}
close IN;

open IN,"$Family.cluster" or die "$!\n";
my %seq=();
while(defined($line=<IN>)){
	chomp $line;
	if($line=~/>/){
		$u=$line;
		$u=~s/>//;
	}else{
		$seq{$u}.=uc($line);
	}
}
close IN;


open IN,$out or die "$!\n";
my %mapped=();
my %E_value=();
while(defined($line=<IN>)){
	chomp $line;
	my @a=split(/\t/,$line);
	my $identity=$a[2]*$a[3];
	my $A=length($seq{$a[0]});
	my $B=$length{$a[1]};
	my $c=$identity/$A;
	my $d=$identity/$B;
	if($c > $Blast_cut and $d > $Blast_cut){
		if(exists $mapped{$a[0]}){
			if($a[-2] < $E_value{$a[0]}){
				$E_value{$a[0]}=$a[-2];	
				$mapped{$a[0]}=$a[1];
			}
		}else{
			$mapped{$a[0]}=$a[1];
		}
	}
}
close IN;
undef%E_value;
undef%length;


#open OUT_1,'>',$prefix.".mappedToknown.gff" or die "$!\n";
my %unmapped=();
my %ID=();
foreach my $j (keys%represent){
	if(exists $mapped{$Represent{$j}}){
		my @P=split(/,/,$represent{$j});
		my $Order=0;
		foreach my $i (@P){	
			my @s=split(/:/,$i);
			my @a=split(/-/,$s[1]);	
			my @b=split(/\|/,$mapped{$Represent{$j}});
			###classification name####
			$b[1]=~s/5S/RSS/;
			$b[1]=~s/7SL/RSL/;
			$b[1]=~s/tRNA/RST/;
			$b[1]=~s/Unknown/RSX/;
#			my $id="$s[0]:$s[1]|$speciesname|$cluster{$i}|$b[1]_KnownBase";
			my $id="$b[1]-$speciesname-$j-$Order#SINE  match=$mapped{$Represent{$j}}";
			$Order++;
			$ID{$i}=$id;
#			print OUT_1 "$s[0]\tSINEscan\tSINE\t$a[0]\t$a[1]\t$s[-1]\tgene_id=",$id,"\n";
		}
	}else{
		$unmapped{$Represent{$j}}=$represent{$j};
	}
}
#close OUT_1;
my $C=keys%cluster;
my $D=keys%mapped;
undef%mapped;

my $flag=0;
open IN,"$Family.cluster" or die "$!\n";
my %seqs_unmapped=();
while(defined($line=<IN>)){
	chomp $line;
	if($line=~/>/){
		$line=~s/>//;
		if(exists $unmapped{$line}){
			$flag=1;
#			while($line=~/\|/){	
#				$line=~s/\|/_/;
#			}
			$u=$line;
			$seqs_unmapped{$u}="";
		}else{
			$flag=0;
		}
	}elsif($flag == 1){
		$seqs_unmapped{$u}.=uc($line);
	}
}
close IN;
my $RRR=keys%Represent;


###tRNA blast###
open IN,$tRNA or die "$!\n";
my %trna=();
while(defined($line=<IN>)){
	chomp $line;
	if($line=~/>/){
		$line=~s/>//;
		$trna{$line}="";
		$u=$line;	
	}else{
		$trna{$u}.=$line;
	}
}
close IN;

my %R=();
my %r_score=();
foreach my $j (keys%seqs_unmapped){
	open OUT_1,'>',"$prefix.one.fa" or die "$!\n";
	my $J=$j;
	my $seq=substr($seqs_unmapped{$j},0,120);
	print OUT_1">$J\n$seq\n";
	my $A=length($seq);
#	my $A=length($seqs_unmapped{$j});
	close OUT_1;
	$R{$j}="RSX";
	print "$j\n";
#	system "mkdir $dir";
	foreach my $k (keys%trna){
		open OUT_2,'>',"$prefix.two.fa" or die "$!\n";
		my $B=length($trna{$k});
		print OUT_2">$k\n$trna{$k}\n";	
		close OUT_2;
		system "/home/maohlzj/EMBOSS-6.6.0/emboss/stretcher -asequence $prefix.one.fa -bsequence $prefix.two.fa -outfile $prefix.align -aformat fasta";
#		system "perl Sine_tRNA.pl align $A $B >>$trna_con";
		my $AlignFile=$prefix.".align";
		my $align=Bio::AlignIO->new(-file => $AlignFile,-format=>'fasta');
		my $stat=Bio::Align::PairwiseStatistics->new();
		my $aln=$align->next_aln();
		my $gap=$stat->number_of_gaps($aln);
		my $difference=$stat->number_of_differences($aln);
		my $L=$stat->number_of_comparable_bases($aln);
			
		my $P1=($L-$difference)/$A;
		my $P2=($L-$difference)/$B;
		if($P2 > $rna_one and $L > $rna_two){ 
			$_=$k;
			if(exists $R{$j}){
				if($r_score{$j} < $P2){
					$R{$j}=$k;
					$r_score{$j}=$P2;
				}
			}else{
				$R{$j}=$k;
				$r_score{$j}=$P2;
			}
		#	last;
		}
	}
	system "rm $prefix.one.fa $prefix.two.fa $prefix.align";
}
undef%seqs_unmapped;
undef%trna;

###unmapped sine annotation###

foreach my $j (keys%unmapped){
	my @s=split(/,/,$unmapped{$j});
	###classification name####
	my $Order=0;
	foreach my $i (@s){
		my $id="";
		if(exists $R{$j}){
			my $rn=$R{$j};
			if($rn=~/tRNA/){
				$rn="RST";
			}elsif($rn=~/5S/){
				$rn="RSS";
			}elsif($rn=~/7SL/){
				$rn="7SL";
			}else{
				$rn="RSX";
			}
			$id="$rn-$speciesname-$cluster{$j}-$Order#SINE  match=NEW";
			$Order++;
		}else{
			my $rn="RSX";
			$id="$rn-$speciesname-$cluster{$j}-$Order#SINE  match=NEW";
			$Order++;
		}
		$ID{$i}=$id;
	}
}
undef%R;
undef%mapped;

####output three files gff,all fasta, represent fasta####
print "$outputDir\n";
open OUT_3,'>',"$outputDir/".$speciesname.".represent.sine.fasta" or die "$!\n";
foreach my $i (sort{$Num{$a}<=>$Num{$b}} keys%Num){
	my @p=split(/,/,$represent{$i});	
	foreach my $line (@p){
		if($Represent{$cluster{$line}} eq $line){
			$flag=1;
			my @A=split(/-/,$ID{$line});
			$_=$ID{$line};
			(my $r)=/\s+(.*)/;
			print OUT_3 ">$A[0]-$A[1]-$A[2]#SINE  $r\n$seq{$line}\n";
		}
	}
}
undef%ID;
undef%Represent;
undef%cluster;
undef%represent;
close OUT_3;
system "rm $prefix.tab $Family.*";
