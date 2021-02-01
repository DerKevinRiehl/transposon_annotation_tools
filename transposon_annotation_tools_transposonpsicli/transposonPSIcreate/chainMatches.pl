#!/usr/local/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case bundling);

my ($btabFile, $btabDir, $maxOverlap, $maxSegDist, $swap_db_query, $follow_status);


my $usage = <<__EOUSAGE__;



################# Chain Matches  ##############################
##                                                           ##
####            For use with Transposon-PSI btab files     ####
#
#  parameters:
#  
#  --btabFile      : a single btab file
#    or
#  --btabDir       :directory containing the btab files (all search results from the same asmbl_id)
#
#  *Advanced options:
# 
#  --maxOverlap   :max bp's of overlap between 
#                    adjacent hsps (default 50)
#  --maxSegDist   :max bp's between chained segments (default 500)
#
# Flags:
#  --swap_db_query  :with other (non-TPSI searches), where the genome sequence is the query, set this option.l
#  --follow_status  :shows percentage done.
#################################################################

The algorithms proceeds in two steps:

    1.  the maximum scoring chain of non-overlapping matches is found for all
          hits on a genome sequence.

    2.  neighboring hits in the same orientation (query and match) are chained together.

## All btabs should include search results for only a single assembly!!!!!!  


__EOUSAGE__


    ;

&GetOptions ("btabFile=s" => \$btabFile,
             "btabDir=s" => \$btabDir,
             "maxOverlap=i" => \$maxOverlap,
             "maxSegDist=i" => \$maxSegDist,
             "swap_db_query" => \$swap_db_query,
             "follow_status" => \$follow_status,
             
             );


my $MAX_PREV_COMPARISONS = 1000;


unless ($btabDir || $btabFile) {
    die $usage;
}

unless (defined($maxOverlap)) {
    $maxOverlap = 50;
}
unless (defined($maxSegDist)) {
    $maxSegDist = 500;
}

if ($maxOverlap < 0 || $maxSegDist < 0) {
    die $usage;
}


my @btabs;
if ($btabFile) {
    @btabs = ($btabFile);
} else {
    @btabs = glob ("$btabDir/*btab");
}

my @structs;

## struct definition:
##   btab => btab_line
##   db_acc => transposon acc
##   contig_acc => genome sequence acc
##   orient => [+-]
##   genome_lend => lend
##   genome_rend => rend
##   db_lend => lend
##   db_rend => rend
##   score => blast score (bits)
##   path_score => Dynamic programming path score 
##   link => pointer to the previous struct in the path to its highest score


## Read each btab file and create structs for each hit:
foreach my $btab (@btabs) {
    open (BTAB, $btab) or die $!;
    while (<BTAB>) {
	unless (/\w/) { next;}
	if (/^error/i) { next;}
	my $btab_line = $_;
	chomp;
	my @x = split (/\t/);
	my ($db_acc, $genome_acc, $db_lend, $db_rend, $genome_end5, $genome_end3, $score) = ($x[0], $x[5], $x[6], $x[7], $x[8], $x[9], $x[12]);
	
    if ($swap_db_query) {
        ($db_acc, $genome_acc) = ($genome_acc, $db_acc);
        ($db_lend, $db_rend, $genome_end5, $genome_end3) = ($genome_end5, $genome_end3, $db_lend, $db_rend);
    }
    

	if ($db_lend > $db_rend) {
	    die "Error, db coords reversed: $db_lend, $db_rend !! ";
	}

	my $orient = ($genome_end5 < $genome_end3) ? '+' : '-';

	my ($genome_lend, $genome_rend) = sort {$a<=>$b} ($genome_end5, $genome_end3);
	my $struct = {
	    btab => $btab_line,
	    db_acc => $db_acc,
	    contig_acc => $genome_acc,
	    orient => $orient,
	    genome_lend => $genome_lend,
	    genome_rend => $genome_rend,
	    db_lend => $db_lend,
	    db_rend => $db_rend,
	    score => $score,
	    path_score => $score,
	    link => 0
	    };
	push (@structs, $struct);
    }
    close BTAB;
}


@structs = sort {$a->{genome_lend}<=>$b->{genome_rend}} @structs;

## Chain hits using DP, find the highest scoring path
for (my $i = 1; $i <= $#structs; $i++) {

    my $struct_i = $structs[$i];
    my ($genome_lend_i, $genome_rend_i, $score_i, $path_score_i) = ($struct_i->{genome_lend}, 
								    $struct_i->{genome_rend}, 
								    $struct_i->{score},
								    $struct_i->{path_score});


    if ($follow_status) {
        my $percent_done = sprintf ("%.2f", $i / $#structs * 100);
        print STDERR "\r   $percent_done   ";
    }
    
    for (my $j = $i-1; $j >= 0; $j--) {
	

        if ($i - $j > $MAX_PREV_COMPARISONS) { 
            last;
        }
        
        ## J always comes before I (except in alphabet)
        
	my $struct_j = $structs[$j];
	
	my ($genome_lend_j, $genome_rend_j, $path_score_j) = ($struct_j->{genome_lend}, 
							      $struct_j->{genome_rend}, 
							      $struct_j->{path_score});
	
	$genome_rend_j -= $maxOverlap;
	
	if ($genome_rend_j < $genome_lend_i) {
	    # possible link
	    
	    my $score = $score_i + $path_score_j;
	    if ($score > $path_score_i) {
		$path_score_i = $struct_i->{path_score} = $score;
		$struct_i->{link} = $struct_j;
	    }
	}
    }
}


## find the maximum path score:
my $max_score = 0;
my $max_scoring_struct = 0;

foreach my $struct (@structs) {
    my $path_score = $struct->{path_score};
    if ($path_score >= $max_score) {
	$max_score = $path_score;
	$max_scoring_struct = $struct;
    }
}

## get all nodes in path
my @pathStructs;
my $link = $max_scoring_struct;
while ($link != 0) {
    push (@pathStructs, $link);
    $link = $link->{link};
}

# re-sort for genome end5
@pathStructs = sort {$a->{genome_lend}<=>$b->{genome_lend}} @pathStructs;

# chain together neighboring matches in the same orientation and order:
foreach my $struct (@pathStructs) {
    $struct->{link} = 0;
}
for (my $i = 1; $i <= $#pathStructs; $i++) {

    my $first_struct = $pathStructs[$i-1];
    my $first_db_acc = $first_struct->{db_acc};
    my $first_orient = $first_struct->{orient};
    my $first_db_lend = $first_struct->{db_lend};
    my $first_db_rend = $first_struct->{db_rend};
    my $first_genome_rend = $first_struct->{genome_rend};
    
    my $second_struct = $pathStructs[$i];
    my $second_db_acc = $second_struct->{db_acc};
    my $second_orient = $second_struct->{orient};
    my $second_db_lend = $second_struct->{db_lend};
    my $second_db_rend = $second_struct->{db_rend};
    my $second_genome_lend = $second_struct->{genome_lend};

    if ($first_db_acc eq $second_db_acc && $first_orient eq $second_orient && 
	
	(## foward orientation  (allow for small $maxOverlap region of overlap)
	 ($first_orient eq '+' && $first_db_rend - $maxOverlap < $second_db_lend) 
	 ||
	 ## reverse orientation
	 $first_orient eq '-' && $second_db_rend - $maxOverlap < $first_db_lend)
 	
	## within segment distance limit on genome
	&& $second_genome_lend - $first_genome_rend <= $maxSegDist
	
	) { 
	$second_struct->{link} = $first_struct;
    }
}


## print final entries, separating the chains:
for (my $i = 0; $i <= $#pathStructs; $i++) {
    if ($pathStructs[$i]->{link} == 0 && $i != 0) {
	print "\n"; #separate chains with newlines
    }
    print $pathStructs[$i]->{btab};
}

exit(0);


    
    



