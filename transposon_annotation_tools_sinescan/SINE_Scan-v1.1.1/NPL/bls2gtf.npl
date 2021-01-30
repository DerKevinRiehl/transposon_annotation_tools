use strict;
if(@ARGV<1){
	print "$0 <bls>\n";
	exit(0);
}

my %entry;
open fin,"<$ARGV[0]" or die $!;
while(<fin>){
	chomp;
	if(/^\s*$/){
		next;
	}
	my @x=split(/\t/,$_);
	my $id=$x[0];
	my $flag='+';
	if($x[9]-$x[8]<0){
		$flag='-';
		my $t=$x[8];
		$x[8]=$x[9];
		$x[9]=$t;
	}
	$entry{$id}{"$x[1]\tTE\tclass/family\t$x[8]\t$x[9]\t.\t$flag\t."}=1;
}
close fin;

foreach my $id(sort keys %entry){
	my $n=0;
	foreach my $info(sort keys %{$entry{$id}}){
		++$n;
		print "$info\tgene_id \"$id.$n\";\n";
	}
}
