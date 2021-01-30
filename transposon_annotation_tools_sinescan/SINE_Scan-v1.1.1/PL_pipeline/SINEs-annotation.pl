#!usr/bin/perl
use strict;
use Getopt::Long;

###total parameters in pipeline of SINEs-annotaion####
my($help,$TEfile,$genomefile,$out_prefix,$species_name,$identity_coverage,$cdhit_identity,$cdhit_length,$trna_identity,$trna_overlap,$outputDir,$cpu,$log);
Getopt::Long::GetOptions(
	"h" => \$help,
	"i=s" => \$TEfile,
	"g=s" => \$genomefile,
	"o=s" => \$out_prefix,
	"n=s" => \$species_name,
	"D=f" => \$identity_coverage,
	"c=f" => \$cdhit_identity,
	"a=f" => \$cdhit_length,
	"e=f" => \$trna_identity,
	"p=i" => \$trna_overlap,
	"z=s" => \$outputDir,
	"k=i" => \$cpu,
	"l=s" => \$log
);

help() if $help;
if($TEfile eq "" or $genomefile eq "" or $out_prefix eq "" or $species_name eq "" or $log eq ""){
	warn "please use -h to know some message\n";
	exit 0;
}

if($identity_coverage eq ""){
	$identity_coverage=0.8;
}
$identity_coverage*=100;

if($cdhit_identity eq ""){
	$cdhit_identity=0.8;
}

if($cdhit_length eq ""){
	$cdhit_length=0.8;
}

if($trna_identity eq ""){
	$trna_identity=0.6;
}

if($trna_overlap eq ""){
	$trna_overlap=60;
}

if($cpu eq ""){
	$cpu=2;
}


###Recheck a-box and b-box###
=pod
system "perl /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/ABbox_check.pl $TEfile";
system "cat $TEfile >$TEfile.tmp2";
system "cat $TEfile.tmp >$TEfile";
system "cat $TEfile.tmp2 >$TEfile.tmp";
system "wc $TEfile >$TEfile.note";
open in,$TEfile.".note" or die "$!\n";
my $Empty=<in>;
if($Empty=~/^0/){
	print "no SINE candidtates passes RNA-verification!\n";
	system "rm $TEfile.tmp2 $TEfile.note";
	exit 0;
}
close in;
system "rm $TEfile.tmp2 $TEfile.note";
=cut
#######

system "/usr/bin/makeblastdb -in $TEfile -dbtype nucl";
system "/usr/bin/makeblastdb -in /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/SINEBase/SineDatabase.fasta -dbtype nucl";

#######main pipeline####################
my $SINEs=$out_prefix.".sines";
my $seq=$out_prefix.".seq.sine.fa";
my $List=$out_prefix.".assignRegionToTE.normal";
system "perl /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/getStrongHit.pl $TEfile $genomefile $out_prefix $cpu >$log";
if(-s $List){
	system "perl /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/sines-extract.anno.pl $List $genomefile >$seq";
	system "perl /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/classification-stronghit.anno.pl $seq $out_prefix $species_name $cdhit_identity $cdhit_length $identity_coverage $trna_identity $trna_overlap $outputDir $TEfile >>$log";
} 
system "rm $out_prefix.*";
#########################################
sub help(){
	warn "please read these information\n";
print STDERR <<EOF;
	Usage: perl SINEs-annotation.pl -i inputfile(list) -o out-prefix;
	Additional options:
	-h	:help information
	-i	:input file SINEs fasta 
	-g	:genome file
	-o	:out-prefix
	-n	:species_name,genus_species_version/strain,care '_' to sperate name
	-D	:sines sequence identity percentage,default is 0.9
	-c	:cdhit identical cutoff,default is 0.8
	-a	:cdhit length overlap, default is 0.8
	-e	:trna identity percentage,default is 0.6
	-p	:number trna comparable bases(overlap,only match and mismatch without gaps) ,default is 60
	-z	:the directory of SINEs annotation output files
	-k	:the cpu number, default is 2
	-l	:the path of Module Three's logfile
EOF
}
