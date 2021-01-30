###a draft pipeline for search SINEs###
#!usr/bin/perl
use strict;
use lib './modules';
use Parallel::ForkManager;
use Benchmark qw(:all);
use File::Basename;
use Getopt::Long qw(:config no_ignore_case bundling);

###total parameters in this pipeline####
my $one_start=Benchmark->new;
my($help,$cpu_num,$step,$genome,$workdir,$name,$infile,$annotation,$min_copy_num,$extend_size,$end_size,$MITE_base,$cutoff,$cdlength,$B_identity,$RNAi,$RNAc,$outputDir,$min_D,$max_flank,$min_M,$min_con,$siteLen,$H_block,$H_sites,$continue);
GetOptions(
	'h'=>\$help,
	'k=i'=>\$cpu_num,
	's=i'=>\$step,
	'g=s'=>\$genome,###genomic file
	'd=s'=>\$workdir,
	'o=s'=>\$name,
	'i=s'=>\$infile,###user provided unvalidated SINE candidate file
	'a=s'=>\$annotation,###user provided varified SINEs file
	'n=i'=>\$min_copy_num,
	'F=i'=>\$extend_size,
	'E=i'=>\$end_size,
	'c=f'=>\$cutoff,###for cd-hit identical percentage
	'p=f'=>\$cdlength,###length coverage for cd-hit
	'b=f'=>\$B_identity,####for blast filter
	't=f'=>\$RNAi,####for RNA identical percentage
	'r=i'=>\$RNAc,####for RNA length overlap
	'z=s'=>\$outputDir,
	'D=f'=>\$min_D,###minimum difference
	'S=f'=>\$max_flank,###maximum score of flanking sequences
	'I=f'=>\$min_M,###minimum score of intermediate sequences
	'C=f'=>\$min_con,###minimum identical percent of a column
	'l=i'=>\$siteLen,###maximum distance between two conserved sites
	'L=i'=>\$H_block,###minimum length of a high conserved block
	'H=f'=>\$H_sites,###minimum percent of conserved sites in high conserved block
	'w=i'=>\$continue
);

my $boundary=1;###beta paramter, revising
my @para=();
if($continue == 2 || $continue == 3){
	if(-e "para.log"){
		print "SINE_Scan will be started from Step $continue.\n";
		print "Reading previous parameter settings. New parameters in command line are ignored. If you really need to change parameters, do this in the file para.log\n";
		open IN,"para.log" or die "$!\n";
		while(<IN>){
			chomp $_;
			my @a=split(/=/,$_);
			push @para,$a[1];
		}
		close IN;
###parameter setting########
		if($continue == 3){
			$step=3;
		}else{
			$para[0]=~s/1//;
			$step=$para[0];
		}
		$genome=$para[1];
		$workdir=$para[2];
		$name=$para[3];
		$infile=$para[4];
		if($step=~/2/ and !-e $infile){
			print"Error. Run SINE_Scan from Step $continue requires -i parameter is properly set.\n";
			exit(1);
		}
		if($step!~/2/){
			$annotation=$para[5];
			$infile="";
			if(!-e $annotation){
				print"Error. Run SINE_Scan from Step $continue requires -a parameter is properly set.\n";
				exit(1);
			}
		}else{
			if(-e $para[5]){
				system "rm $para[5]";
			}
			$annotation="";
		}
		$min_copy_num=$para[6];
		$extend_size=$para[7];
		$end_size=$para[8];
		$cutoff=$para[9];
		$cdlength=$para[10];
		$B_identity=$para[11];
		$RNAi=$para[12];
		$RNAc=$para[13];
		$outputDir=$para[14];
		$boundary=$para[15];
		$cpu_num=$para[16];
		$min_D=$para[17];
		$max_flank=$para[18];
		$min_M=$para[19];
		$min_con=$para[20];
		$siteLen=$para[21];
		$H_block=$para[22];
		$H_sites=$para[23];
		if($genome eq "" or $workdir eq "" or $name eq "" or $step eq ""){
			print"Error. The file para.log contains invalid values for some parameters.\n";
			exit(1);
		}
	}else{
		print"Error. The file para.log cannot be found.\n";
		exit(1);
	}
}


if($help or $genome eq "" or $workdir eq "" or $name eq "" or $step eq ""){
	help();
}

if($step eq ""){
	$step=123;
}

if($cpu_num eq ""){
	$cpu_num=2;
}

if($boundary eq ""){
	$boundary=0;
}

if($min_copy_num eq ""){
	$min_copy_num=5;
}

if($extend_size eq ""){
	$extend_size=60;
}

if($end_size eq ""){
	$end_size=25;
}

if($cutoff eq ""){
	$cutoff=0.8;
}

if($cdlength eq ""){
	$cdlength=0.8;
}

if($B_identity eq ""){
	$B_identity=0.8;
}

if($RNAi eq ""){
	$RNAi=0.6;
}

if($RNAc eq ""){
	$RNAc=60;
}

if($outputDir eq ""){
	$outputDir='./';
}

if($min_D eq ""){
	$min_D=0.3;
}

if($max_flank eq ""){
	$max_flank=0.6;
}

if($min_M eq ""){
	$min_M=0.75;
}

if($min_con eq ""){
	$min_con=0.8;
}

if($siteLen eq ""){
	$siteLen=4;
}

if($H_block eq ""){
	$H_block=10;
}

if($H_sites eq ""){
	$H_sites=0.6;
}

if($continue eq ""){
	$continue=1;
}
##########main pipeline#############
#####some trivial judgments#####
if($step!~/1/ and $step=~/2/ and $infile eq ""){
	print "Error. -i parameter should be properly set if begin from Step Two.\n";
	exit(1);
}

if($step!~/2/ and $step=~/3/ and $annotation eq ""){
	print "Error. -a parameter should be properly set if begin from Step Three.\n";
	exit(1);
}

if(-d "$workdir"){
	if($continue == 2 || $continue == 3){
		print "$workdir exists. Rerunning SINE_Scan from Step $continue.\n";
	}else{
		print "$workdir exists. clean it.\n";
		system "rm -rf $workdir/*";
	}
}else{
        print "$workdir does not exists. build it.\n";
	mkdir "$workdir", 0755 or print "cannot create $workdir directory:$!";
}

my $G_name=basename($genome);
my @S=split(/\./,$G_name);
if($S[-1] eq 'fasta'){

}elsif($S[-1] eq 'fas'){
	$G_name=~s/\.fas$/\.fasta/;
}elsif($S[-1] eq 'fa'){
	$G_name=~s/\.fa$/\.fasta/;
}elsif($S[-1] eq 'fna'){
	$G_name=~s/\.fna$/\.fasta/;
}elsif($S[-1] eq 'fsa'){
	$G_name=~s/\.fsa$/\.fasta/;
}else{
	$G_name.=".fasta";
}

if(-f $genome){
	if($continue != 2 and $continue != 3){
		system "cp $genome $workdir/$G_name";
	}
}else{	
	print "Error. $genome doesn't exist or the format isn't FASTA.\n";
	exit(1);
}
#$genome=basename($genome);
my $Genome=$G_name;
$genome="$workdir/$G_name";
my $line=$genome;
if($line=~/fasta$/){
	$line=~s/\.fasta$//;
}
my $sinefile=$line.".sine.fa";
my $SFfile=$line."-matches.fasta";
my $S5Ffile=$line."-5smatches.fasta";
my $logfile2=$line.".moduleTwo.logfile";
my $logfile3=$line.".moduleThree.logfile";
####process user input SINE file which needs to be identified########
if($infile ne ""){	
	if(-f $infile and -e $infile){
		if($continue != 2 and $continue != 3){
			system "cp $infile $sinefile";
		}
	}else{
		print "Error. $infile doesn't exist or the format isn't FASTA.\n";
		exit(1);
	}
}

###########process user input SINE annotation file which needs to be annotated##########
my $Workdir=basename($workdir);
my $anno_file=$Workdir.".for_annotation.fa";########annotation file in workdir
$anno_file="$workdir/$anno_file";
if($annotation ne ""){
	if(-e $annotation and -f $annotation){
		if($continue != 2 and $continue != 3){
			system "cp $annotation $anno_file";
		}
	}else{
		print "Error. The file $annotation cannot be found or its format isn't FASTA.\n";
		exit(1);
	}
}

open fout,">para.log" or die "$!\n";
print fout"step=$step\ngenome=$genome\nworkdir=$workdir\nname=$name\ninfile=$sinefile\nannotation=$anno_file\nmin_copy_num=$min_copy_num\nextend_size=$extend_size\nend_size=$end_size\ncutoff=$cutoff\ncdlength=$cdlength\nB_identity=$B_identity\nRNAi=$RNAi\nRNAc=$RNAc\noutputDir=$outputDir\nboundary=$boundary\ncpu_num=$cpu_num\nmin_D=$min_D\n$max_flank=$max_flank\nmin_M=$min_M\nmin_con=$min_con\nsiteLen=$siteLen\nH_block=$H_block\nH_sites=$H_sites\n";
close fout;

####Step One: ab initial identification of SINE candidates. Currently only SINE-Finder is used here. Other tools (RepeatModeler, RepeatScout etc.) will be added in future. User can run them separately and use -i parameter to perform verification  #######
if($step=~/1/){
	if($infile eq "" and $annotation eq ""){	
		print "Step One: Run SINE-Finder.\n";
		my $pm=Parallel::ForkManager->new(5);
		my @exec=('/mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/7SLandtRNA-sine_finder.py','/mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/5S-sine_finder.py -s 0');
		my $m=0;
		Finder_Loop:
		foreach my $process (@exec){
			my $pid=$pm->start and next Finder_Loop;
			system "python $process -T chunkwise -f fasta $genome";
			$pm->finish;
		}
		$pm->wait_all_children;	
		if(!-f $SFfile and !-f $S5Ffile){
			print "Error. SINE-Finder output file $SFfile doesn't exist.\n";
			exit(1);
		}
		if(-e $SFfile){
			system "perl /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/getSINE-noTSD.pl $SFfile >$sinefile";
		}
		if(-e $S5Ffile){
			system "perl /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/getSINE-noTSD.pl $S5Ffile >>$sinefile";
		}
		print "Finish Step One.\n";
	}elsif($infile ne "" and $annotation ne ""){
		print "Incorrect paramter setting: since validated SINE sequences are provided,  parameter combination -s 3 -a $annotation shoud be used. If you want to validate sequences in file $infile. use -s 23 -i $infile combination\n";
		exit(1);
	}else{
		print "Incorrect paramter setting: since you provide unvalidated SINE sequences, parameter combination -s 23 -i $infile shoud be used.\n";
		exit(1);
	}
}
my $one_stop=Benchmark->new;
my $time_cost=timediff($one_stop,$one_start);
print "Time cost for SINE-Finder Module One:",timestr($time_cost),"\n";

		
		
######################
#Build genome database####
my $second_start=Benchmark->new;
if($step =~/2/ or $step=~/3/){
	my @files=();
	if(-e $genome.".nal"){
		open in,$genome.".nal" or die "$!\n";
		while(<in>){
			if($_=/DBLIST/){
				chomp $_;
				my @s=split(/\s+/,$_);
				shift @s;
				foreach my $i (@s){
					my @a=split(/\./,$i);
					my $w=$genome.".$a[-1]";
					push @files,$w;
				}
			}
		}
		close in;	
	}else{
		push @files,$genome;
	}
	my $flag=0;
	foreach my $j (@files){
		if(-e $j.".nhr" and -e $j.".nin" and -e $j.".nsq"){
			$flag=1;
		}else{
			$flag=0;
		}
	}
	
	if($flag == 1){
		print "Genomic database exists.\n";
	}else{
		print "Build BlAST database for the genomic sequences.\n";
		system "/usr/bin/makeblastdb -in $genome -dbtype nucl";
		if(-e $genome.".nal"){
			open in,$genome.".nal" or die "$!\n";
			open out,'>',$genome.".nal2" or die "$!\n";
			while(<in>){
				while($_=~/$workdir\//){
					$_=~s/$workdir\//\.\//;
				}
				print out"$_";
			}
			close in;
			close out;
			system "mv $genome.nal2 $genome.nal";
		}
	}
}

if($boundary == 1 && $step=~/1/){
	my $RGdir="$workdir/RG";
	mkdir "$RGdir", 0755 or print "cannot create $workdir directory:$!";
	system "perl /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/rg_mainscript.pl $sinefile $genome 5 $RGdir 60 25 $cpu_num >$RGdir/logfile";
	#system "cp -r $workdir/RG .";
	if(-e "$workdir/RG.error.log"){
		print "Your genomic dataset has no reasonable SINE candidates\n";
		exit (0);
	}else{
		system "rm -rf $workdir/RG";
	}
#	system "mv $RGdir RG_total/";
}
my $second_stop=Benchmark->new;
$time_cost=timediff($second_stop,$second_start);
print "Time cost for Boundary fix:",timestr($time_cost),"\n";

###################
###Step Two: checking SINE candidates by identifying sequnence signals of TE amplification###
my $third_start=Benchmark->new;
if($step=~/2/){
	if($annotation eq ""){
		if($sinefile ne "" and -e $sinefile and -s $sinefile and -f $sinefile){
			print "Step Two: SINE Verification.\n";
			system "perl /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/mainpipeline.pl $sinefile $genome $min_copy_num $workdir $extend_size $end_size $min_D $max_flank $min_M $min_con $siteLen $H_block $H_sites $cpu_num >$logfile2";
			print "Finish Step Two.\n";
		}else{
			print "Step Two: Could not read SINE candidate file $sinefile. Please check this file. Stop.\n";
			if(-e "formatdb.log"){
				system "mv formatdb.log $workdir";
			}
			if(-e "error.log"){
				system "mv error.log $workdir";
			}
			exit 0;
		}
	}else{
		print "Since -a parameter is used, -s can only be set as 3. You can use -a $annotation -s 3 combination\n";
		if(-e "formatdb.log"){
			system "mv formatdb.log $workdir";
		}
		if(-e "error.log"){
			system "mv error.log $workdir";
		}
		exit 0;
	}
}

my $third_stop=Benchmark->new;
$time_cost=timediff($third_stop,$third_start);
print "Time cost for Module Two :",timestr($time_cost),"\n";

###################################
####Optional step: manual verification#########
##################################

###Step Three SINEs annotation########
my $fourth_start=Benchmark->new;
if($step=~/3/){
	my $out_anno=$line.".anno";
	if($anno_file ne "" and -e $anno_file and -s $anno_file and -f $anno_file){
		print "Step Three: SINE classification and genome-wide annotation!\n";
		system "perl /mnt/beegfs/home1/miska_copy/riehl/AnnotationTask/SoftwarePackages/BenchmarkAnnotator/sineScan/SINE_Scan-v1.1.1/PL_pipeline/SINEs-annotation.pl -i $anno_file -g $genome -o $out_anno -n $name -z $outputDir -D $B_identity -c $cutoff -a $cdlength -e $RNAi -p $RNAc -k $cpu_num -l $logfile3";
		print "Finish Step Three.\n";
	}else{
		print "Step Three: Cannot read validated SINE candidate file. Stop.\n";
		print "If you ran Step 2, please check the file with the suffix of for_annotation.fa. Empty file means no good candidate was identified.\n";
		print "If you use your own validated SINEs (-a parameter), check this file.\n";
		if(-e "formatdb.log"){
			system "mv formatdb.log $workdir";
		}
		if(-e "error.log"){
			system "mv error.log $workdir";
		}
		exit 0;
	}
}
my $fourth_stop=Benchmark->new;
$time_cost=timediff($fourth_stop,$fourth_start);
print "Time cost for Module Third :",timestr($time_cost),"\n";
print "SINE_Scan pipeline finished. Running mode=$step.\n";
if(-e "formatdb.log"){
	system "mv formatdb.log $workdir";
}
if(-e "error.log"){
	system "mv error.log $workdir";
}
print "SINE_Scan pipeline finished. Running mode=$step.\n";


sub help{

print STDERR <<EOF;

SYNOPSIS: SINE_Scan_process.pl <options>
 
Options:
        -g	string	        Genomic file, the file's suffix must be 'fasta', 'fa' or 'fas'. [Required] 
	-d	string          Directory which the pipeline uses. Intermidiate and final results files will be put in this directory. [Required]
	-o	string          Genome ID used in final output. An example(the authors' advice):genus_species_strain(or version). [Required]
	-z	string		The directory where final results are written. [Default: ./] 
	-k	integer		The number of cpu numbers. [Default: 2]
	-s	integer         Mode of running SINE_Scan. This pipeline has three modules, you can run them all or seperately. To run the entire pipeline, type 123. To run the first two modules, type 12. [default: 123]
	-i	string	        FASTA format file containing SINE candidates for verification. This file is needed if you skip Module One to run Module Two. [Default: NULL]
	-a	string          FASTA format file containing SINE candidates for classification and homology search. This file is needed if you skip Module One and Two to run Module Three directly. [Default: NULL]
	-n	integer	        Minimum copy number of SINEs in genomic data. [Default: 5]
	-F	integer	        Length of sequences flanking to SINE ends. This region will be anlyzed to determine SINE insertion events. [Default: 60]
	-E	integer	        Length of SINE ends used to determine SINE insertion events. [Default: 25]
	-c	float	        Minimum identity of SINEs within a family. This value should be >= 0.8. [Default: 0.8]
	-p	float           Minimum proportion of matched region of two SINEs within a family. [Default: 0.8]
	-b	float           Identity percentage of a shorter sequence of two sequences for blast filter. [default: 0.8]
	-t	float           Minimum identity to define a SINE matching a RNA well. [Default: 0.6]
	-r	integer         Minimum percentage of a RNA matching a SINE, which defines a SINE matching a RNA well. [default: 60]
	-D	float		Minimum differences between SAQ of TE ends and flanking regions. [default: 0.3]
	-S	float		Maximum SAQ of flanking regions. [default: 0.6]
	-I	float		Minimum SAQ of TE ends. [default: 0.75]
	-C	float           Minimum identity in one column to define a conserved site. [Default: 0.8]
	-l	integer		Maximum distance between two conserved sites. [Default: 4]	
	-L	integer		Minimum length of a highly conserved block. [Default: 10]
	-f	float		Minimum percent of conserved sites in highly conserved block. [Default: 0.6]
	-w	integer		Continue to run the pipeline from the point it stops. The value of -w only could be assigned 2 or 3, which means the pipeline begins from Module Two or Three. Use this parameter only based on an aborted previous run. If this paramter is used, no other parameter is needed. 
	-h	string          Show this help

An example : perl SINE_Scan_process.pl -g genome_file(fasta) -d workdir -s 123 -o species_name;

EOF
	exit 0; 
}
