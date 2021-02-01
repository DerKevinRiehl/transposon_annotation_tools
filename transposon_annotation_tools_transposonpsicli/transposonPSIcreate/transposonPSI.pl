#!/usr/local/bin/perl

use lib ($ENV{EUK_MODULES});
use strict;
use Fasta_reader;
use File::Basename;
use Cwd;

my $usage = "usage: $0 fastaFile prot|nuc\n\n";

my $dbDir = "/home/bhaas/CVS/ANNOTATION/EUK_GENOME_DEVEL/RepAnnotator/transposonPSIcreate/TPSI_library";

my $max_Evalue = "1e-5";


main: {

    my $fastaFile = $ARGV[0] or die $usage;
    my $dbType = $ARGV[1] or die $usage;
    unless ($dbType eq "prot" || $dbType eq "nuc") {
	die $usage;
    }
    
    mkdir ("transposonPSI") or die "Error, cannot mkdir transposonPSI";
    
    my $runDir = cwd();
    open (TOPHIT, ">$fastaFile.topHits") or die $!;
    
    my $fasta_reader = new Fasta_reader($fastaFile);
    
    while (my $seqObj = $fasta_reader->next()) {
	my $accession = $seqObj->get_accession();
	my $fasta_seq = $seqObj->get_FASTA_format();
	
	print STDERR "processing $accession.\n";
	
	&run_transposonPSI($accession, $fasta_seq, $dbType);
	
    }

    close TOPHIT;

    exit(0);

}
					
####
sub run_transposonPSI {
    my ($accession, $fasta_seq, $dbType) = @_;
    
    
    my $seqDir = $accession;
    $seqDir =~ s/\W/_/g;
    my $clean_acc = $seqDir;
    
    $seqDir = "transposonPSI/$clean_acc";
    mkdir $seqDir or die "Error, cannot mkdir $seqDir\n";
    
    my $seqFilename = "$seqDir/$clean_acc.seq";
    open (SEQFILE, ">$seqFilename") or die $!;
    print SEQFILE $fasta_seq;
    close SEQFILE;
    
    my $isProt = ($dbType eq "nuc") ? "F" : "T";
    my $cmd = "formatdb -i $seqFilename -p $isProt";
    system $cmd;
    die "Error, $cmd (ret $?)" if $?;
    

    
    my @hit_summary;
    
    foreach my $refseq (<$dbDir/*.refSeq>) {
	
	print STDERR "\tblast against $refseq\n";
	
	my $chkpFile = $refseq;
	$chkpFile =~ s/\.refSeq/\.chkp/;
	my $coreFilename = basename($refseq);
	
	my $outputFile = "$seqDir/$clean_acc.$coreFilename.psitblastn";
	
	my $cmd;
	if ($dbType eq "nuc") {
	    # do psitblastn w/ blastall
	    $cmd = "blastall -i $refseq -d $seqFilename -p psitblastn -R $chkpFile -F F -M BLOSUM45 -t -1 -e $max_Evalue > $outputFile";
	} else {
	    # do psiblast w/ blastpgp
	    $cmd = "blastpgp -i $refseq -d $seqFilename -R $chkpFile -j 1 > $outputFile";
	}

	print "CMD: $cmd\n";
	system $cmd;
	die "Error $cmd $?" if $?;
	
	## covert to btab output
	my $cmd = "BPbtab < $outputFile > $outputFile.btab";
	system $cmd;
	die "Error $cmd $?" if $?;
	
	open (BTAB, "$outputFile.btab") or die $!;
	while (<BTAB>) {
	    unless (/\w/) { next;}
	    unless (/error/i) {
		chomp;
		my @x = split (/\t/);
		push (@hit_summary, [@x]);
	    }
	}
	close BTAB;
	
    }
    
    my $summary_file = "$seqFilename.btab";
    @hit_summary = reverse sort {$a->[12]<=>$b->[12]} @hit_summary;

    my $top_hit;
    open (SUMMARY, ">$summary_file") or die $!;
    foreach my $hit (@hit_summary) {
	my $btab_entry = join ("\t", @$hit);
	unless ($top_hit) {
	    $top_hit = $btab_entry;
	}
	print SUMMARY "$btab_entry\n";
    }
    close SUMMARY;
    
    if ($top_hit) {
	print TOPHIT $top_hit . "\n";
    } else {
	print TOPHIT "// $accession lacks T-PSI match.\n";
    } 
    
}


	


