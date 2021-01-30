########renormalization group to approximate obvious boundaries of specific DNA elements########
#!usr/bin/perl
use strict;
use lib '/mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/modules';
use Statistics::Basic qw(:all);
use Bio::SimpleAlign;
use Bio::AlignIO;
use File::Basename;

if(@ARGV < 2){
	print "$0 <mult fasta file> <name>\n";
	exit(0);
}

my $file=$ARGV[0];
my @a=split(/\//,$file);
my $ListNum=$a[-2];
my ($one,$two)=(dirname($file),basename($file));
$one=~s/\/$a[-2]//;
####global sequences, but I suppose local would be good###
system "/usr/bin/muscle -in $file -out $file.msa.fasta -maxiters 1 -diags -quiet";
my @consensus=Similarity("$file.msa.fasta");

#######boundary detection##########	
my $str = Bio::AlignIO->new(-file => "$file.msa.fasta");
my $aln = $str->next_aln();
my $n=$aln->length;
my $String=$aln->consensus_string(30);
my @A=();
my @B=();
my $flag=-1;
for(my $i=0;$i<$n;$i++){
	if($consensus[$i] >= 0.7){
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
my @island=();####high quality island boundaries#####
for(my $i=0;$i<@B;$i++){
	my @c=split(/,/,$B[$i]);
	if(@c > 10){
		push @island,$i;
#		$W.="$B[$i],";
	}
}

my $stop_point=$island[-1];
for(my $j=0;$j<@island-1;$j++){
	my $k=$island[$j];
	my $K=$island[$j+1];
	my @c=split(/,/,$B[$k]);
	my @C=split(/,/,$B[$K]);	
	if(abs($C[0]-$c[-1]) > 20){
		my $str=substr($String,$c[-1],$C[0]-$c[-1]+1);
		my $gapR=gaps_compute($str);
		if($gapR == 0){
			$stop_point=$k;
			last;
		}
	} 	
}

@island=();
for(my $i=0;$i<$stop_point+1;$i++){
	my @c=split(/,/,$B[$i]);
	if(@c > 10){
		push @island,$i;
		$W.="$B[$i],";
	}
}

my @coordinate=split(/,/,$W);
my @sign=();
push @sign,$coordinate[0];
push @sign,$coordinate[-1];	
####boundary condition: close enough and meidate all gaps##########
###OK, I do loop (iteration) here######
#print "@B\n@sign\n";
for(my $i=$island[0]-1;$i>-1;$i--){
	my @c=split(/,/,$B[$i]);
	my $string=substr($String,$c[-1]+1,$coordinate[0]-$c[-1]-1);
#	print "$consensus[134]\n$i\n@c\n$string\n";
#	exit 0;
#	print "1,$string\n";
	my $R=good_block($B[$i]);
	my $r=gaps_compute($string);	
	if(($r == 1 && $R == 0) or ($R == 2 && length($string) < 20)){
	#if(($r == 1 && $R == 0) or ($R == 2)){
#	if($string!~/[ACGTagct]/){
		$coordinate[0]=$c[0];
		push @sign,$c[0];
		push @sign,$c[-1];
	}else{
		last;
	}	
}
for(my $i=$island[-1]+1;$i<@B;$i++){
	my @c=split(/,/,$B[$i]);
#	if($c[0]-$coordinate[-1] < 10){
	my $string=substr($String,$coordinate[-1]+1,$c[0]-$coordinate[-1]-1);
	my $R=good_block($B[$i]);
	my $r=gaps_compute($string);	
	if(($r == 1 && $R == 0) or ($R == 2 && length($string) < 20)){
	#if(($r == 1 && $R == 0) or ($R == 2)){
#	if($string!~/[ACGTagct]/){
		$coordinate[-1]=$c[-1];
		push @sign,$c[0];
		push @sign,$c[-1];
	}else{
		last;
	}	
#	}
}
@sign=sort{$a<=>$b} @sign;
if(scalar(@sign) < 2){
	print "0\n";
}
my $start=$sign[0];
my $end=$sign[-1];
#print "@consensus\n";
#print "$String\n";
my $R=substr($String,$start,$end-$start+1);
while($R=~/\?/){
	$R=~s/\?//;
}
$R=uc($R);
my $H=substr($R,0,10);
my $E=substr($R,-10,10);
if($H=~/^\w{0,2}[A|T]{5,}/){
	if($E=~/[A|T]{5,}\w{0,2}$/){
		$H=~s/A|T//g;
		$E=~s/A|T//g;
		if(length($H) < length($E)){
			$R=reverse($R);
			$R=~tr/ATCG/TAGC/;	
		}
	}else{
		$R=reverse($R);
		$R=~tr/ATCG/TAGC/;	
	}
}
while($R=~/[A|T]{5,}\w{0,2}$/ or $R=~/^\w{0,2}[A|T]{5,}/){
	$R=~s/[A|T]{5,}\w{0,2}$//;
	$R=~s/^\w{0,2}[A|T]{5,}//;
}
###7SL tRNA
my $a_box="[GA][CGAT]TGG|TGGCTCACGCC|T[AG]G[CT]\\w{2}A\\w{3,4}G";
my $spacer_1="\\w{25,70}"; 
my $b_box="GTTC[AG]A|GTTCGAGAC|G[AT]TC[AG]A\\w{2}C";
##5S		
my $A_box="[ATC][AG]G[CT][CT]AAGC";
my $Spacer_1="\\w{20,50}";
my $B_box="[AG]TGG[AG][ATG]GAC";
my $s=$R;
my $flag=0;
foreach my $i ($aln->each_seq){
	my $s = $i->subseq($start+1,$end+1);###5 side
	while($s=~/-{1,}/){
		$s=~s/-{1,}//;
	}
	my $rs= reverse($s);
	$rs=~tr/TACG/ATGC/;
	$s=substr($s,0,110);
	$rs=substr($rs,0,110);
	if($s=~/($a_box)($spacer_1)($b_box)/ or $s=~/($A_box)($Spacer_1)($B_box)/){
		$flag=1;	
		last;	
	}elsif($rs=~/($a_box)($spacer_1)($b_box)/ or $rs=~/($A_box)($Spacer_1)($B_box)/){
		$flag=1;
#		$R=reverse($R);
#		$R=~tr/ATCG/TAGC/;	
		last;	
	}
}

my $len_consensus=length($R);
if($flag == 0){
	$len_consensus=0;
}
#$end=$start+$len_consensus-1;
#print "$start $end\n";
#print "$R\n";
print "$len_consensus\n";
####here we shoulg give conditions to judge convergence based on ends of this cluster#####
my $repfile=$file;
$repfile=~s/extendseq/fa/;
my $seqname=$ARGV[1];
$seqname=~s/,/\|/g;
system "echo '>$seqname' >$repfile";
system "echo $R >>$repfile";

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

sub gaps_compute{
	my $str=$_[0];
	my $b=length($str);	
	while($str=~/\?/){
		$str=~s/\?//;
	}
	my $a=length($str);
	my $percent=int($a*100/$b)/100;
	if($percent <= 0.2){
		return 1;
	}else{
		return 0;
	}
}

sub good_block{
	my $str=$_[0];
	my @a=split(/,/,$str);
	my $block_size=$a[-1]-$a[0]+1;
	my $high_points=@a;
	my $percent=int($high_points*100/$block_size)/100;
	if($percent < 0.5 or $high_points < 2){
		return 1;
	}elsif($percent >= 0.8){
		return 2;
	}else{
		return 0;
	}
}
