#!/usr/local/bin/perl

use lib ("PerlLib", ## common functions 
         "TransposonPerlLib"); ## custom functions for here.

use DBI;
use Mysql_connect;
use CGI;
use CGI::Carp qw(fatalsToBrowser); 
use strict;
use Data::Dumper;
use TransposonDB_API;


$|++;

my $cgi = new CGI();
print $cgi->header();

my %params = $cgi->Vars();
my $genbank_gi = $params{genbank_gi} or die "Error, require genbank_gi\n";
my $delta = $params{delta};

my ($dbproc) = &connect_to_db("bhaas-lx","TransposonDB","access","access");

print $cgi->start_html();
print "<a href='TransposonDB.cgi'>TransposonDB</a>";
print "<center>\n";


my ($genbank_acc, $organism_name, $seqTitle);
my $query = "select genbank_acc, genbank_gi, organism_name, seqTitle from GenomeSeqRecord where genbank_gi = ?";
my $result = &first_result_sql($dbproc, $query, $genbank_gi);
if ($result) {
    ($genbank_acc, $genbank_gi, $organism_name, $seqTitle) = @$result;
    print "<h1>$seqTitle</h1><h2>($genbank_acc, gi:$genbank_gi) of <i>$organism_name</i></h2>\n";
} else {
    die "Error, no record for gi:$genbank_gi\n";
}


## Illustration of region:




## get list of TEs for this sequence:
my $query = "select tf.feat_id, tf.end5, tf.end3 from TransposonFeature tf, GenomeSeqRecord gs where gs.transposon_seq_ID = tf.gseq_id and tf.feat_type = 'TE' and gs.genbank_gi = ? order by tf.end5";
my @results = &do_sql_2D($dbproc, $query, $genbank_gi);
foreach my $result (@results) {
    my ($TE_id) = @$result;
    my $proteins_fasta_text;


    print "<img src=\"illustrate_region.cgi?genbank_gi=$genbank_gi&TE=$TE_id&delta=$delta\" alt=\"TE-image\"/>\n";

    print "<table border=1>" 
	. "<tr><th colspan=3 bgcolor='#006666'><font color='#ffffff'>Transposon($TE_id)</font></th></tr>\n"
	. "<tr><th>Feature Type</th><th>Coordinates</th><th>Annotations</th></tr>\n";
    ## get TE components:
    my $query = "select tf.feat_id, tf.feat_type, tf.end5, tf.end3 from TransposonFeature tf, TransposonLink tl where tl.component_ID = tf.feat_id and tl.element_ID = ? order by tf.end5";
    my @results = &do_sql_2D($dbproc, $query, $TE_id);
    foreach my $result (@results) {
	my ($feat_id, $feat_type, $end5, $end3) = @$result;
	my $coordstring = "$end5-$end3";
	
	## get annot features:
	my $query = "select qualifier, annotText from FeatureAnnots where feat_id = ?";
	my @results = &do_sql_2D($dbproc, $query, $feat_id);
	my %feature_annots;
	foreach my $result (@results) {
	    my ($qualifier, $annotText) = @$result;
	    $feature_annots{$qualifier} = $annotText;
	}
	

	if ($feat_type eq "CDS") {
	    $coordstring = "";
	    my $query = "select end5, end3 from CDS_coords where orf_ID = ?";
	    my @results = &do_sql_2D($dbproc, $query, $feat_id);
	    foreach my $result (@results) {
		my ($end5, $end3) = @$result;
		$coordstring .= "$end5-$end3,";
	    }
	    chop $coordstring;
	    $coordstring = "$coordstring";
	    
	    ## get the orf data:
	    my $query = "select sequence from FeatureSequence where feat_id = ? and seqType = ?";
	    my $translation = &very_first_result_sql($dbproc, $query, $feat_id, 'P');
	    	    
	    $translation =~ s/(\w{60})/$1\n/g;
	    my $gi = $feature_annots{db_xref};
	    my $acc = $feature_annots{protein_id};
	    my $note = $feature_annots{note};
	    my $header = "$acc $gi $note {$organism_name}";
	    $header =~ s/\s+/ /g;
	    $header =~ s/^\s//;
	    $proteins_fasta_text .= ">$header\n$translation\n";
	    
	}
	my $annot_text = "<pre>";
	foreach my $key (sort keys %feature_annots) {
	    my $val = $feature_annots{$key};
	    $annot_text .= "$key => $val\n";
	}
	$annot_text .= "</pre>";
	
	print "<tr><td>$feat_type</td><td>$coordstring</td><td>$annot_text</td></tr>\n";
    }
    print "</table>\n";
    print "</center>\n";
    print "<pre>$proteins_fasta_text</pre>\n";
}



print $cgi->end_html();

$dbproc->disconnect;


