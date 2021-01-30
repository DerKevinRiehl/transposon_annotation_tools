####here we want to find TSDs and check boundary######
#!usr/bin/perl
use strict;
use lib './modules';
use Statistics::Basic qw(:all);
use Bio::SimpleAlign;
use Bio::AlignIO;

if(@ARGV<3){
	print "$0 <mult fasta file> <size Flank> <size TE end>\n";
	exit(0);
}

my $file=$ARGV[0];
my @a=split(/\//,$file);
my $PassEndList="$a[0]/pass.manual.checkList";
my $ListNum=$a[1];
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
system "/home/maohlzj/sine_te/softwares/muscle -in $outfile.fasta -out $outfile.msa.fasta -maxiters 1 -diags -quiet";
####find 60 50 60 positions in MSA #####
my @positionA=();
my @positionB=();
my $str = Bio::AlignIO->new(-file => "$outfile.msa.fasta");
my $aln = $str->next_aln();
my $nseq= $aln->no_sequences;
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
my $SINEs=0;
open OUT,">>$PassEndList" or die "$!\n";
if(($L<=0.6 && $R<=0.6 && $M>=0.75) || ($M-$L>=0.3 && $M-$R>=0.3)){
	$SINEs=1;
#	print "L=$L,H=$M,R=$R\n";
	print OUT"$ListNum\n";
}
close OUT;

####TSD-finder ####
my $state=0;
if($SINEs == 1){
	####it pass the boundary condition######
	my $String=$aln->consensus_string(30);
	my @A=();
	my @B=();
	my $flag=-1;
	for(my $i=0;$i<$n;$i++){
		if($consensus[$i] >= 0.8){
			push @A,$i;
		}
	}
	my $r="$A[0],";
	for(my $j=1;$j<@A;$j++){
		###here I should consider high region and low region's boundary####
		if($A[$j]-$A[$j-1] > 3){
			push @B,$r;
			$r="";
		}
		$r.="$A[$j],";
	}
	push @B,$r;
	my $W="";
	for(my $i=0;$i<@B;$i++){
		my @c=split(/,/,$B[$i]);
		if(@c > 10){
			$W.="$B[$i],";
		}
	}
	my @coordinate=split(/,/,$W);
	my @sign=();
	push @sign,$coordinate[0];
	push @sign,$coordinate[-1];	
	####boundary condition: close enough and meidate all gaps##########
	for(my $i=0;$i<@B;$i++){
		my @c=split(/,/,$B[$i]);
		if(@c <= 10 && $c[-1] < $coordinate[0]){
			my $string=substr($String,$c[-1]+1,$coordinate[0]-$c[-1]-1);
			if($string!~/[ACGTagct]/){
				push @sign,$c[0];
				push @sign,$c[-1];
			}	
		}
		if(@c <= 10 && $c[0] > $coordinate[-1] && $c[0]-$coordinate[-1] < 5){
			my $string=substr($String,$coordinate[-1]+1,$c[0]-$coordinate[-1]-1);
			if($string!~/[ACGTagct]/){
				push @sign,$c[0];
				push @sign,$c[-1];
			}	
		}
	}
	@sign=sort{$a<=>$b} @sign;
	###TSD boundary condition####?
	my $start=$sign[0];
	my $end=$sign[-1];
	my $R=substr($String,$start,$end-$start+1);
#	print "$start,$end\n$String\n@B\n$R\n";
	####TSD up and down###
	my $mite=MITE_test($R);
	if($mite == 1){
		$state=2;
		print "$state\n";
	}
	my $counter=0;
	my $nseq= $aln->no_sequences;
	my %order=();########each sequence 's name - boundary position#########
	
	foreach my $seq ($aln->each_seq) {
       		if($state == 2){
			last;
		}
		my $res = $seq->subseq(1,$start);
		my $name=$seq->display_name();
		while($res=~/-/){
			$res=~s/-+//;
		}
		my $S=length($res);
		my $Res = $seq->subseq($end,$n);####3 side
		while($Res=~/-/){
			$Res=~s/-+//;
		}
		my $E=length($Res);
		$order{$name}="$S\t$E";
		####TSD length condition #############
		if($S < 40 or  $E < 40){
			$state=-1;
			print "$state\n";
			last;
		}
	}
	foreach my $seq ($aln->each_seq) {
		if($state == -1 || $state == 2){
			last;
		}
       		my $res = $seq->subseq(1,$start);
		while($res=~/-/){
			$res=~s/-+//;
		}
		$res=substr($res,-15,15);####5 side
		my $Res = $seq->subseq($end,$n);####3 side
		while($Res=~/-/){
			$Res=~s/-+//;
		}
#######moving frame to detect TSD#################moving the 'TSD' of 5' side over 3' side sequence###########
		for(my $i=6;$i<25;$i++){
			my $s=0;
			if($i-15 > 0){
				$s=$i-15;
			}else{
				$s=0;
			}
			my $ser=substr($res,-$i,$i);
			my $Ser=substr($Res,$s,$i-$s);
			open fout, ">$file.tsd.fa" or die "$!\n";
			print fout">one\n$ser\n>two\n$Ser\n";
			close fout;
			my $str = Bio::AlignIO->new(-file => "$file.tsd.fa");
			my $ALN=$str->next_aln();
			my $cons=$ALN->consensus_string(80);	
			my $CONS=$cons;
#			print "$ser\n$Ser\n$CONS\n\n";	
			while($cons=~/\?/){
				$cons=~s/\?//;
			}
			my $OKbase=0;
			my $percent=length($cons)/length($ser);
			if($i < 8 and $CONS=~/\w{4,}/ and $percent > 0.7){
				$OKbase=1;
			}elsif($CONS=~/\w{6,}/ and $percent > 0.5){
				$OKbase=1;
			}
			if($OKbase == 1){
#				print "$ser\n$Ser\n$CONS\n";	
				$counter++;
				last;
			}
		}
		if($counter == 2 or ($nseq < 7 and $counter == 1)){
			$state=1;
			print "$state\n";
			open fout,">$file.order" or die "$!\n";
			foreach my $i (keys%order){
				print fout">$i\n$order{$i}\n";
			}
			close fout;
			last;
		}
	}
	undef%order;
	if(-e "$file.tsd.fa"){
		system "rm $file.tsd.fa";
	}
}else{
	print "$L,$M,$R\n";
}
#system "rm $outfile.fasta $outfile.msa.fasta";

sub Similarity{
        my $fl=$_[0];
        my $str = Bio::AlignIO->new(-file => $fl);
        my $aln = $str->next_aln();
#	print $aln->consensus_string(0)."\n";
	my $nseq= $aln->no_sequences;
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

sub MITE_test{
	my $str=$_[0];
	while($str=~/\?/){
		$str=~s/\?//;
	}
	my $head=substr($str,0,10);
        $head=reverse($head);
        $head=~tr/ATCG/TAGC/;
        $head=~tr/atcg/tagc/;
	my $tail=substr($str,-10,10);
	my $a=0;###head start
	my $b=0;###tail start
	for(my $i=5;$i<15;$i++){
		if($i >= 10){
			$a=0;
			$b=$i-10;
		}else{
			$a=10-$i;
			$b=0;
		}	
		my $len=10-abs(10-$i);
		my $one=substr($head,$a,$len);
		$one=~s/N/M/;
		my $two=substr($tail,$b,$len);
		###consensus OK
		my $con="";
		for(my $j=0;$j<$len;$j++){
			my $p=substr($one,$j,1);
			my $q=substr($two,$j,1);
			if($p eq $q){
				$con.=$p;
			}else{
				$con.='?';
			}
		}
		if($con=~/\w{5,}/){
			return 1;
#			last;
		}
	}
	return 0;
}
