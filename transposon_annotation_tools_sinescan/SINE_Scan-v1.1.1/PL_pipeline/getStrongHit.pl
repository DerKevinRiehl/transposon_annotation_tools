## this program merge overklapping regions that match TEs
#then use the best hit in TE base for name of that region
use strict;
use File::Basename;
if(@ARGV<4){
	print "$0 <TE fasta file> <genomic sequence> <output prefixs> <cpu>\n";
	print "Be sure TE fasta file and genome file have been formated to blast db\n";
	exit(0);
}

my $script=dirname($0);
my $te=$ARGV[0];
my $genome=$ARGV[1];
my $prefix=$ARGV[2];
my $cpu=$ARGV[3];

print "find all homologous region in genomes\n";
system "/usr/bin/blastn -task blastn -db $genome -query $te -max_target_seqs 100000 -evalue 1e-5 -dust no -out $prefix.te2genome.bls -num_threads $cpu -outfmt 6";
print "filter blast to get region match >9";
#system "perl ./PL_pipeline/solar/solar.pl -cCd -1 $prefix.te2genome.bls >$prefix.te2genome.solar";
open in,$te or die "$!\n";
my %len=();
my $u;
while(<in>){
	chomp $_;
	if($_=~/>/){
		my @a=split(/\s+/,$_);
		$u=$a[0];
		$u=~s/>//;
		$len{$u}=0;
	}else{
		$len{$u}+=length($_);
	}
}
close in;

open fin,"<$prefix.te2genome.bls" or die $!;
open fout,">$prefix.te2genome.filter" or die $!;
while(<fin>){
	my @x=split(/\t/,$_);
	###mao correct it for solar question 20150329
	my $query_length=$len{$x[0]};
	if($x[7]-$x[6]+1>=0.9*$query_length && abs($x[9]-$x[8])+1>=0.90*$query_length){
		print fout $_;
	}
}
close fin;
close fout;

print "get GFF file\n";
system "perl /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/bls2gtf.pl $prefix.te2genome.filter >$prefix.te2genome.gtf";


my $gff=basename("$prefix.te2genome.gtf");
if(!-s "$prefix.te2genome.gtf"){
	print "these SINE candidates don't have any strong copies.\n";
	exit 0;
}
print "Sort GFF file\n"; 
system "sort -k1,1 -k4,4n $prefix.te2genome.gtf >$prefix.sortGFF";
print "Merge overlapping regions using bedtools\n";
system "/mnt/home1/miska/riehl/anaconda3/envs/python27/bin/bedtools merge -i $prefix.sortGFF >$prefix.mergebed";
print "Get sequences of these merge regions\n";
system "/mnt/home1/miska/riehl/anaconda3/envs/python27/bin/bedtools getfasta -fi $genome -bed $prefix.mergebed -fo $prefix.mergebed.fasta";
print "Blast these merged region to TE base. Only the first best hit is kept\n";
system "/usr/bin/blastn -task blastn -db $te -query $prefix.mergebed.fasta -max_target_seqs 1 -evalue 1e-5 -out $prefix.mergebed.TEBase.bls -num_threads $cpu -outfmt 6";
print "Get info of best hits (blast output format)\n";
#system "perl ./PL_pipeline/solar/solar.pl $prefix.mergebed.TEBase.bls >$prefix.assignRegionToTE";
###read length###
if(!-s "$prefix.mergebed.TEBase.bls"){
	print "these SINE candidates don't have any strong copies.\n";
	exit 0;
}
my %len=();
my $u;
open in,$te or die "$!\n";
while(<in>){
	chomp $_;
	if($_=/>/){
		$u=$_;
		$u=~s/>//;
		$len{$u}=0;
	}else{
		$len{$u}+=length($_);
	}
}
close in;

open in,"$prefix.mergebed.fasta" or die "$!\n";
while(<in>){
	chomp $_;
	if($_=/>/){
		$u=$_;
		$u=~s/>//;
		$len{$u}=0;
	}else{
		$len{$u}+=length($_);
	}
}
close in;
####
print "check complex region. These regions are much longer than its best matched TE (cutoff = 2*max(TE size))\n";
#They are extracted to a additional fasta file for further inspection

#open fin,"<$prefix.assignRegionToTE" or die $!;
open fin,"<$prefix.mergebed.TEBase.bls" or die "$!\n";
open fout1,">$prefix.assignRegionToTE.normal" or die $!;###useful now
open fout2,">$prefix.assignRegionToTE.complex" or die $!;
while(<fin>){
#	my @x=split(/\s+/,$_);
	my @x=split(/\t/,$_);
	if($len{$x[0]}>2*$len{$x[1]} || $len{$x[1]}<0.5*$len{$x[0]}){
		print fout2 $_;
	}else{
		print fout1 $_;
	}
}
close fout1;
close fout2;
close fin;

print "Done.\n";
