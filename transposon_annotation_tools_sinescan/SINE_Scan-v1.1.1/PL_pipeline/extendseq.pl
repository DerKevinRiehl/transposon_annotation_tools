use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO;


my $flank=1000;

if(@ARGV<2){
	print "$0 <bls> <Genome> [extendsize]\n";
	exit(0);
}

if($ARGV[2]=~/^\d+$/){
	$flank=$ARGV[2];
}

my $database=$ARGV[1];
my $seqin = Bio::SeqIO->new(-format=>'fasta', -file=>$database);
my %genome;
while(my $seq = $seqin->next_seq){
	$genome{$seq->id}=$seq->seq;
}

my $n=0;
open fin,"<$ARGV[0]" or die $!;
while(<fin>){
	my @line=split(/\s+/,$_);
	my $extendseqs;
	my $strand='+';
	if($line[9]-$line[8]<0){
		my $t=$line[8];
		$line[8]=$line[9];
		$line[9]=$t;
		$strand='-';
	}
	if($line[8]-$flank>=1){
		$extendseqs=substr($genome{$line[1]},$line[8]-1-$flank,$line[9]-$line[8]+1+2*$flank);
	}else{
		$extendseqs=substr($genome{$line[1]},0,$line[9]+$flank);
	}
	$extendseqs=uc($extendseqs);
	if($strand eq "-"){
		$extendseqs=reverse($extendseqs);
		$extendseqs=~tr/ATCG/TAGC/;
	}
	$n+=1;
	$line[12]=~s/,/-/g;
	$line[12]=~s/;$//;
	$line[12]=~s/;/../g;
	my $id="$line[5]|$line[7]--$line[8]|$line[12]|$line[4]";
	print ">$line[0]\_$n  $line[1]:$line[8]:$line[9]  $flank  $strand\n";
#	print ">$id  flank=$flank $line[0]\n";
	print FormatSequence(\$extendseqs,50);
}
close fin;

sub FormatSequence
{
    my($seq,$colume)=@_;

    my $len=length($$seq);
    my $out;
    my $i;
    for($i=0;$i<($len);$i+=$colume)
    {
         $out.=substr($$seq,$i,$colume)."\n";
    }
    return $out;
}
