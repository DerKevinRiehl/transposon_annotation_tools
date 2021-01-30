use strict;
if(@ARGV<1){
	print "$0 <SINE-FINDER output>\n";
	exit(0);
}

my %strand;
$strand{"F"}="+";
$strand{"R"}="-";

open fin, "<$ARGV[0]" or die $!;
my $suffix="";
if($ARGV[0] =~/5smatch/){
	$suffix="5S";
}
my $label;my $seq;
$/=">";
$label=<fin>;
$/="\n"; 
while($label=<fin>){
	$label =~ s/^>//;
	$label =~ s/\s*$//;               
	$/=">";
	$seq=<fin>;  
	$/="\n";
	$seq =~ s/>$//;
	$seq =~ s/\s+//g;
	my @lab=split(/\s+/,$label);
	my $scaf=$lab[0];
	my $chain=$strand{$lab[1]};
	my $pos=$lab[2];
	$pos=~s/:/-/;
	my $tsd=0;
	if($lab[3]=~/TSD-len=(\d+);/){
		$tsd=$1;
	}
	my $len=length($seq);
	my $seq=substr($seq,$tsd,$len-2*$tsd);
	my $id="$scaf|$pos|$chain";
	if($suffix ne ""){
		$id=$id."|".$suffix;
	}
	print ">$id\n$seq\n";
}
close fin;
