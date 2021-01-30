####here we want to find TSDs and check boundary######
#!usr/bin/perl
use strict;
use lib '/mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/modules';
use Statistics::Basic qw(:all);
use File::Basename;
use Bio::SimpleAlign;
use Bio::AlignIO;

if(@ARGV<3){
	print "$0 <mult fasta file> <size Flank> <size TE end>\n";
	exit(0);
}

my $file=$ARGV[0];
my @a=split(/\//,$file);
my $script=dirname($file);
$script=~s/$a[-2]$//;
my $PassEndList="$script/passList";
my $ListNum=$a[-2];
my $sizeFlank=$ARGV[1];
my $sizeEnd=$ARGV[2];
my $outfile=$file.".SixfyFiftySixty";

open fout,">$outfile.fasta" or die "$!\n";
open fin,"<$file" or die $!;
$/=">";
my $label=<fin>;
$/="\n";
while($label=<fin>){
	$label =~ s/^>//;
        $label =~ s/\s*$//;
        $/=">";
        my $seq=<fin>;
        $/="\n";
        $seq =~ s/>$//;
        $seq =~ s/\s+//g;
        my $head=substr($seq,0,$sizeFlank);
        my $tail=substr($seq,-1*$sizeFlank,$sizeFlank);
        my $homo=substr($seq,$sizeFlank,$sizeEnd).substr($seq,-1*$sizeFlank-$sizeEnd,$sizeEnd);
	my $total=$head.$homo.$tail;
#	my $total=$head.$tail;
	my $L=length($total);
	my @A=split(/\s+/,$label);
	print fout">$label L$sizeFlank M$sizeEnd R$sizeFlank\n$total\n";
}
close fout;
close fin;
system "/usr/bin/muscle -in $outfile.fasta -out $outfile.msa.fasta -maxiters 1 -diags -quiet";
####find 60 50 60 positions in MSA #####
my @positionA=();
my @positionB=();
my $str = Bio::AlignIO->new(-file => "$outfile.msa.fasta");
my $aln = $str->next_aln();
my $nseq= $aln->num_sequences;
foreach my $seq ($aln->each_seq) {
	my $count=0;
	for(my $pos=1;$pos<=$aln->length;++$pos){
       		my $res = $seq->subseq($pos, $pos);
		if($res ne '-'){
			$count++;
		}	
		if($count == 60){
			push @positionA,$pos;
		}elsif($count == 110){
			push @positionB,$pos;
		}
	}
}

my $meanA=int(0.5+mean(@positionA));####left truncated point
my $meanB=int(0.5+mean(@positionB));####right truncated point

#my ($consensus,$score)=split(/\t/,Similarity("$outfile.msa.fasta"));
my @consensus=Similarity("$outfile.msa.fasta");
my $n=@consensus;
my $L=0;##left region
my $M=0;##mediate region
my $R=0;###right region
my @score=();
for(my $i=1;$i<$n;$i++){
	push @score,$consensus[$i];
	if($i == $meanA){
		$L=mean(@score);
		@score=();
	}
	if($i == $meanB){
		$M=mean(@score);
		@score=();
	}
}
$R=mean(@score);
######judgement condition########
######high identical region and low identical region,TSD finder####
open OUT,">>$PassEndList" or die "$!\n";
print "L=$L,H=$M,R=$R\n";
if(($L<=0.6 && $R<=0.6 && $M>=0.75) || ($M-$L>=0.3 && $M-$R>=0.3)){
#	print "L=$L,H=$M,R=$R\n";
	print OUT"$ListNum\n";
}
close OUT;


sub Similarity{
        my $fl=$_[0];
        my $str = Bio::AlignIO->new(-file => $fl);
        my $aln = $str->next_aln();
#	print $aln->consensus_string(0)."\n";
	my $nseq= $aln->num_sequences;
        my @consensus;
        my @conperc;
        for(my $pos=1;$pos<=$aln->length;++$pos){
                my %count;
                foreach my $seq ($aln->each_seq) {
                        my $res = $seq->subseq($pos, $pos);
                        $count{$res}+=1;
                }
                my $cons=(reverse sort {$count{$a}<=>$count{$b}} keys %count)[0];
                my $perc=$count{$cons}/$nseq;
		my $gap=$count{'-'}/$nseq;
		delete $count{'-'};
		my $n=keys%count;
		if($nseq >= 8 && $n <= 2 && $perc < 0.8 && $cons ne '-' && $gap < 0.2){
			$perc=0.8;
		}
		#######here question????#######
                if($cons!~/[ACGTagct]/){
                        $cons="N";
			$perc=0;
                }
	       	push(@consensus,$cons);
        	push(@conperc,$perc);
        }
	return(@conperc);
}

