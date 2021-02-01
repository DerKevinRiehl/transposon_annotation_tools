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
print $cgi->header;
print $cgi->start_html();


my ($dbproc) = &connect_to_db("bhaas-lx","TransposonDB","access","access");

my $query = "select gsr.genbank_acc, gsr.genbank_gi, gsr.organism_name, gsr.seqTitle, tf.feat_id, tf.end5, tf.end3 from GenomeSeqRecord gsr, TransposonFeature tf where gsr.transposon_seq_ID = tf.gseq_id and tf.feat_type = 'TE' order by gsr.transposon_seq_ID, tf.end5";
my @results = &do_sql_2D($dbproc, $query);

print "<table border=1><tr><th>gb_acc</th><th>gi_num</th><th>org_name</th><th>description</th><th>feat_id</th><th>end5</th><th>end3</th></tr>\n";
foreach my $result (@results) {
    my ($genbank_acc, $genbank_gi, $organism_name, $description, $feat_id, $end5, $end3) = @$result;
    print "<tr><td>$genbank_acc</td><td><a href=\"TE_report.cgi?genbank_gi=$genbank_gi\">$genbank_gi</a></td><td>$organism_name</td><td>$description</td><td>$feat_id</td><td>$end5</td><td>$end3</td></tr>\n";
}
print "</table>\n";

print $cgi->end_html();


$dbproc->disconnect;

exit;
