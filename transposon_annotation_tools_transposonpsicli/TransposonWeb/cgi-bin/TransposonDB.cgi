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



my $cgi = new CGI();

print $cgi->header;
print $cgi->start_html();

print <<__EOHTML;

<center>
<h1>TransposonDB</h1>
</center>
<ul>
    <li><a href="list_database_entries.cgi">Database Contents</a></li>
    <li><a href="transposon_blast.cgi">Blast TransposonDB</a></li>
    <li><a href="genbank_to_TransposonDB_loader.cgi">Import a Genbank Entry</a></li>
</ul>

__EOHTML

    ;

print $cgi->end_html();

exit(0);




