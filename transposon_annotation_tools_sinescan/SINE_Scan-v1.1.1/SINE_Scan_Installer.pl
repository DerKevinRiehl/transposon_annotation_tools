#!usr/bin/perl
use Getopt::Std;
#-----------------------------------------------------

getopts("d:a:f:b:c:M:l:e:S:R:n:h:");

$SINE_Scan_Dir   = defined $opt_d ? $opt_d : ".";
$SINEFINDER      = defined $opt_a ? $opt_a : "";
$Formatdb        = defined $opt_f ? $opt_f : "";
$Blastall        = defined $opt_b ? $opt_b : "";
$cd_hit_est      = defined $opt_c ? $opt_c : "";
$Muscle          = defined $opt_M ? $opt_M : "";
$Bedtools	 = defined $opt_l ? $opt_l : "";
$stretcher       = defined $opt_e ? $opt_e : "";
$SINEs_dataFile  = defined $opt_S ? $opt_S : "./SINEBase/SineDatabase.fasta";
$RNA_dataFile    = defined $opt_R ? $opt_R : "./RNABase/RNAsbase.fasta";
$Help            = defined $opt_h ? $opt_h : "";

if($Help or $SINEFINDER eq "" or $SINE_Scan_Dir eq "" or $cd_hit_est eq "" or $Muscle eq "" or $stretcher eq ""
   or $Formatdb eq "" or $Blastall eq "" or $SINEs_dataFile eq "" or $RNA_dataFile eq "" or $Bedtools eq ""){
	usage();
}

my $sinefinder1="";##7sl and tRNA SINE
my $sinefinder2="";##5S RNA

###check parameters####
if(! -d $SINE_Scan_Dir){
	print "Error. Make sure that SINE_Scan directory is input correctly.\n";
	exit(1);
}

if(! -f $SINEs_dataFile){
	print "Error. Make sure that Known_SINEs data file with FASTA format is input correctly.\n";
	exit(1);
}

if(! -f $RNA_dataFile){
	print "Error. Make sure that RNA data file with FASTA format is input correctly.\n";
	exit(1);
}

if(! -e $SINEFINDER){
	print "Error. Make sure that the path of SINE-FINDER.py is correctly input and you have the authority to run it.\n";
	exit(1);
}else{
	system "dos2unix $SINEFINDER";
	$sinefinder1=$SINE_Scan_Dir."/PL_pipeline/7SLandtRNA-sine_finder.py";
	my $sine_finder_patch=$SINE_Scan_Dir."/SINE_Finder-v1.1-7SLandtRNA.patch";
	system "cp $SINEFINDER $sinefinder1";
	system "patch $sinefinder1 <$sine_finder_patch";
	$sinefinder2=$SINE_Scan_Dir."/PL_pipeline/5S-sine_finder.py";
	my $sine_finder_patch=$SINE_Scan_Dir."/SINE_Finder-v1.1-5sRNA.patch";
	system "cp $SINEFINDER $sinefinder2";
	system "patch $sinefinder2 <$sine_finder_patch";
}

if(! -x $Formatdb){
	print "Error. Make sure that the path of makeblastdb is correctly input and you have the authority to run it.\n";
	exit(1);
}

if(! -x $Blastall){
	print "Error. Make sure that the path of blastn is correctly input and you have the authority to run it.\n";
	exit(1);
}

if(! -x $Bedtools){
	print "Error. Make sure that the path of bedtools2 is correctly input and you have the authority to run it.\n";
	exit(1);
}

if(! -x $cd_hit_est){
	print "Error. Make sure that the path of cd-hit-est is correctly input and you have the authority to run it.\n";
	exit(1);
}

if(! -x $stretcher){
	print "Error. Make sure that the path of EMBOSS stretcher is correctly input and you have the authority to run it.\n";
	exit(1);
}

if(! -x $Muscle){
	print "Error. Make sure that the path of muscle is correctly input and you have the authority to run it.\n";
	exit(1);
}
#system "cp $SINE_Scan_Dir/NPL/*.npl $SINE_Scan_Dir";
@NPL_Files=('ABbox_check.npl','rg_mainscript.npl','RG_boundary.npl','CheckEnd.npl','clusterSeqs.npl','betterSeq-seeds.npl','New_checkboundaryAndTSD.npl','classification-stronghit.anno.npl','extendseq.npl','getSINE-noTSD.npl','getStrongHit.npl','mainpipeline.npl','makeDirforTE.npl','SINEs-annotation.npl','SINE_Scan_process.npl','sines-extract.anno.npl','bls2gtf.npl');
foreach(@NPL_Files) {
	if(-e "$SINE_Scan_Dir/NPL/$_"){	
	}else{
		print "$_ can not be found! Please check the NPL directory.\n";
		exit 0;
	}
}

#-----------------------------------------------------
#@Raw_Files = glob "*.npl";
foreach(@NPL_Files) {
	$NPL = $_;
	open(RF, "<$SINE_Scan_Dir/NPL/$NPL")||die"$!\n";
	$PL = $NPL;
	$PL =~ s/\.npl/\.pl/;
	open(PL, ">$SINE_Scan_Dir/$PL")||die"$!\n";
	while(<RF>) {
		chomp;
		$Line = $_;
		$Line =~ s/_scan_/$SINE_Scan_Dir/;
		$Line =~ s/_sinefinder1_/$sinefinder1/;
		$Line =~ s/_sinefinder2_/$sinefinder2/;
		$Line =~ s/_formatdb_/$Formatdb/;
		$Line =~ s/_blastall_/$Blastall/;
		$Line =~ s/_cd_hit_est_/$cd_hit_est/;
		$Line =~ s/_muscle_/$Muscle/;
		$Line =~ s/_bedtools_/$Bedtools/;
		$Line =~ s/_stretcher_/$stretcher/;
		$Line =~ s/_SINEsbase_/$SINEs_dataFile/;
		$Line =~ s/_RNAbase_/$RNA_dataFile/;		
		print(PL "$Line\n");
	}
	close(RF);
	close(PL);

	system "chmod 755 $PL\n";
}
#system "mv $SINE_Scan_Dir/*.npl $SINE_Scan_Dir/NPL";
system "mv $SINE_Scan_Dir/*.pl  $SINE_Scan_Dir/PL_pipeline";
system "mv $SINE_Scan_Dir/PL_pipeline/SINE_Scan_process.pl $SINE_Scan_Dir";
system "mv $SINE_Scan_Dir/PL_pipeline/SINE_Scan_Installer.pl $SINE_Scan_Dir";
print "Install SINE_Scan finished.\n\n";

sub usage {
    print "\nUsage of SINE_Scan_Installer.pl\n";
    print STDERR <<"    _EOT_";

SYNOPSIS: SINE_Scan_process.pl <options>

All parameters, except for -n and -h, are required, we advise users that you should use the absoulte path.

	-d	string		The directory of SINE_Scan pipeline
	-a	string		The path of SINE-FINDER.py 
	-f	string		The path of "makeblastdb" command, which belongs to blast+, not legacy blast
	-b	string		The path of "blastn" command, which belongs to blast+, not legacy blast
	-c	string		The path of "cd-hit-est" command
	-M	string		The path of "muscle" command
	-l	string		The path of  "bedtools" command
	-e	string		The path of EMBOSS "stretcher" command
	-S	string		The Known SINEs database file, you could give your own databasefile,like "/dir/SINEsdir/SINEs.fa". To use the database provided by SINE_Scan,input 'SINE_Scan/directory/SineDatabase.fasta'
	-R	string		FASTA format file containing tRNA,5SrRNA and 7SLRNA. You can use your own database like "/dir/RNAdir/RNA.fa", To use the database provided by SINE_Scan,input 'SINE_Scan/directory/RNAbase.fasta'
	-h	string          show usage information

An example: SINE_Scan_Installer.pl -d /home/SINE/SINE_Scan/ -a /home/software/SINE-FINDER.py -f /usr/local/ncbi-blast+/makeblastdb -b /usr/local/ncbi-blast+/blastn -c /home/software/cd-hit/cd-hit-est -M /home/software/muscle -l /home/software/bedtools2/bin/bedtools -e /home/EMBOSS-6.6.0/emboss/stretcher -S /home/SINE/SINE_Scan/SINEdatabase/Sines.fasta -R /home/SINE/SINE_Scan/RNAdatabase/RNAs.fasta

****NOTE****
Input the FULL PATH (i.e. that begins with root /) for each paremeter, NOT relative path (see above example).

Please make sure that correct path of each prereqired program is provided. Otherwise, errors will be reported when the program is called. If that happens, you should reinstall the pipeline again using this script.


    _EOT_
    exit(1)
}
