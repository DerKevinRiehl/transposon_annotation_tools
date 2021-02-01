#!/usr/local/bin/perl

use lib ("PerlLib", ## common functions 
	 "TransposonPerlLib"); ## custom functions for here.

use DBI;
use Mysql_connect;
use CGI;
use strict;
use Bio::DB::GenBank;
use Data::Dumper;
use TransposonDB_API;

$| = 1;

our $DB_SEE = 0;

open (STDERR, ">&STDOUT");

my $cgi = new CGI;
my %params = $cgi->Vars;

print $cgi->header;

## arg processing
my $LOAD_ENTRY = 0;
my $GB_ACC = $params{GB_ACC};
$GB_ACC =~ s/\s//g;

## See if we need to create a new TE (ie. load entry)
foreach my $key (keys %params) {
    if ($key =~ /checkbox/ && $params{$key}) {
	$LOAD_ENTRY = 1;
	last;
    }
}



main: {

    print $cgi->start_html();
    print "<center>\n";
    
    &print_accession_form();
    
    if ($LOAD_ENTRY) {
	&process_Transposon_Entry();
    } elsif ($GB_ACC) {
	&check_exists_already($GB_ACC);
	&print_GB_ACC_report($GB_ACC);
    }
    
    print $cgi->end_html();
    exit(0);
    
}


####
sub print_accession_form {
    print <<__EOACCFORM;
    
    <form name=gb_acc_form action="genbank_to_TransposonDB_loader.cgi" method=get >
    Genbank Accession: <input type=text value="" name="GB_ACC" />
    <input type=submit value="Get Record" />
    <input type=reset value="clear" />
    </form>
    <hr>


__EOACCFORM

}



####
sub print_GB_ACC_report {
    print "<form name=genbankReport action=\"genbank_to_TransposonDB_loader.cgi\" method=post />\n";
    print "<table border=1 width=100% >\n";
    print "<tr><td align=center >Report for $GB_ACC </td></tr>\n";
    my $gb = new Bio::DB::GenBank();
    # this returns a Seq object :
    if ($gb) {
	my $seq1 = $gb->get_Seq_by_acc($GB_ACC);
	
	my $alphabet = $seq1->alphabet();
	if ($alphabet ne "dna") {
	    die "<b><font color='#ff0000'>Error, you've requested an $alphabet record.  Only dna records are to be loaded here.</font></b>\n";
	}
	
	my $accession = $seq1->accession();
	my $gi_num = $seq1->primary_id();
	print "<tr><td>Accession: $accession, GI:$gi_num</td></tr>\n";
	print &hidden_ele("accession", $accession);
	print &hidden_ele("gi_number", $gi_num);
	
	my $species = $seq1->species();
	my $organism_name = $species->common_name();
	print "<tr><td>Organism: $organism_name</td></tr>\n";
	print &hidden_ele("organism", $organism_name);
	
	my $description = $seq1->description();
	print "<tr><td>Description: $description</td></tr>\n";
	print &hidden_ele("description", $description);
	print "<tr><td bgcolor='#006666'>&nbsp;</td></tr>\n";
	print "<tr><td><table border=1 width=100% >\n";
	print "<tr><th>Type</th><th>Attributes</th><th>Coordinates</th></tr>\n";
	

	my $feat_count = 0;
	foreach my $feat ($seq1->all_SeqFeatures()) {
	    $feat_count++;
	    dump_feature ($feat, $feat_count);
	}
	print "</table></td></tr>\n";
	my $sequence = $seq1->seq();
	print "<input type=hidden name=sequence value=\"$sequence\" />\n";

    } else {
	print "Sorry, no report found for $GB_ACC\n";
    }

    print "<tr><td align=center bgcolor='#000000' ><input type=submit value=\"Load new TE\" /> <input type=reset value=clearForm /> </td></tr>\n";
    print "</table>\n</form>\n";
    

}


##
sub dump_feature {
    my ($feat, $feat_num) = @_;
    
    my $feat_primary_tag = $feat->primary_tag();
    print "<tr><td><input type=checkbox name=\"$feat_num:checkbox\" unchecked />&nbsp;" . $feat_primary_tag . "&nbsp;&nbsp;</td>";
    print &hidden_ele ("$feat_num:primary_tag", $feat_primary_tag);
    print "<td>";
    my $hidden_tags_text = "";

    foreach my $tag ( $feat->get_all_tags() ) {
	my @values = $feat->get_tag_values($tag);
	if ($tag eq "translation") {
	    print "$tag => length: " . length($values[0]) . "\n";
	} else {
	    
	    print "$tag => " . join(' ',@values) . "<br>\n";
	}
	$hidden_tags_text .= &hidden_ele("$feat_num:$tag", $values[0]);
    }
    print "</td>";

    $hidden_tags_text .= &hidden_ele("$feat_num:simple_coords", $feat->start() . "-" . $feat->end());
    
    my @locations = $feat->location->each_Location();
    my @coords;
    foreach my $location (@locations) {
        push (@coords, $location->start() . "-" . $location->end());
    }
    my $coordString = join (", ", @coords);
    print "<td>$coordString</td>";
    $hidden_tags_text .= &hidden_ele("$feat_num:coords", $coordString);

    print "</tr>\n";
    print $hidden_tags_text;
    

}


####
sub hidden_ele {
    my ($token_name, $token_val) = @_;
    return ("<input type=hidden name=\"$token_name\" value=\"$token_val\" />\n");
}

####
sub process_Transposon_Entry {
    print "Processing Entry!!\n";
    #print "<pre>\n" . &strip_html_chars(Dumper (\%params));

    my %assemblyFeatures;
    my %feature_indices_to_load;
    my %feature_index_to_feature_data;
    foreach my $key (keys %params) {
	my $value = $params{$key};
	if ($key =~ /:/) {
	    # feature component
	    my ($feature_index, $feature_attribute) = split (/:/, $key, 2);
	    if ($feature_attribute eq "checkbox") { 
		if ($value) {
		    $feature_indices_to_load{$feature_index} = 1;
		}
	    } else {
		$feature_index_to_feature_data{$feature_index}->{$feature_attribute} = $value;
	    }
	} else {
	    #core attribute
	    $assemblyFeatures{$key} = $value;
	}
    }

    ## Load the Features
    &Load_TE (\%assemblyFeatures, \%feature_indices_to_load, \%feature_index_to_feature_data);
}


####
sub strip_html_chars {
    my $text = shift;
    $text =~ s/[<>&]/?/g;
    return ($text);
}


####
sub Load_TE {
    my ($assemblyFeatures_href, $feature_indices_to_load_href, $feature_index_to_feature_data_href) = @_;
    my ($dbproc) = &connect_to_db("bhaas-lx", "TransposonDB","access","access");
    
    ## check accession and GI-number to see if it already exists:
    my $accession = $assemblyFeatures_href->{accession};
    my $gi_number = $assemblyFeatures_href->{gi_number};
    die "Error, no accession plus gi_number ($accession, $gi_number)\n" unless ($accession && $gi_number);
    
    my $description = $assemblyFeatures_href->{description};
    
    my $organism_name = $assemblyFeatures_href->{organism};


    my $sequence = $assemblyFeatures_href->{sequence} or die "Error, need genomic sequence.\n";
    unless (-d "genomeSequences") {
	mkdir ("genomeSequences") or die "Error, cannot create directory genomeSequences.\n";
	chmod (0755, "genomeSequences");
    }
    my $fasta_file = "genomeSequences/$gi_number.seq";
    unless (-s $fasta_file) {
	open (FILE, ">$fasta_file") or die "Error, cannot write $fasta_file\n";
	my $seqLen = length($sequence);
	$sequence =~ s/(\w{60})/$1\n/g;
	print FILE ">$gi_number $accession $organism_name len=$seqLen\n$sequence";
	close FILE;
	chmod (0444, $fasta_file);
	## run repfind on it:
	my $cmd = "repfind -f -p -l 10 $fasta_file > $fasta_file.repfind";
	my $ret = system ($cmd);
	if ($ret) {
	    print "Sorry, refind may have failed on sequence $fasta_file (ret: $ret)\n";
	}
	chmod (0444, "$fasta_file.repfind");
    }
    
    
    my $transposon_seq_ID = &TransposonDB_API::load_GenomeSeqRecord($dbproc, 
								    'genbank_acc' => $accession,
								    'genbank_gi' => $gi_number,
								    'organism_name' => $organism_name,
								    'seqTitle' => $description);
    
    ## Create the TE
    my $TE_id = &TransposonDB_API::load_TransposonFeature($dbproc, 
							  'gseq_id' => $transposon_seq_ID, 
							  'feat_type' => 'TE',
							  'end5' => 0, #tmp coords, update later.
							  'end3' => 0);
    


    
    
    my @all_coords;
    my %orient_counter;
    ## Add the individual features of the TE:

    
    foreach my $feature_id (keys %$feature_indices_to_load_href) {
	my $feature_ref = $feature_index_to_feature_data_href->{$feature_id};
	
	my $feat_type = $feature_ref->{"primary_tag"} or die "Error, no primary tag for feature $feature_id\n" . Dumper ($feature_ref);
	
	my $simple_coords = $feature_ref->{"simple_coords"};
	my ($simple_end5, $simple_end3) = split (/-/, $simple_coords);
	if ($simple_end5 < $simple_end3) {
	    $orient_counter{'+'}++;
	} else {
	    $orient_counter{'-'}++;
	}
	
	my $feat_id = &TransposonDB_API::load_TransposonFeature ($dbproc, 
								 'gseq_id' => $transposon_seq_ID,
								 'feat_type' => $feat_type,
								 end5 => $simple_end5,
								 end3 => $simple_end3);


	&TransposonDB_API::load_FeatureAnnots ($dbproc, $feat_id, %$feature_ref);
	
	push (@all_coords, $simple_end5, $simple_end3);

	&TransposonDB_API::load_TransposonLink($dbproc, 'element_ID' => $TE_id, 'component_ID' => $feat_id);
	
	if ($feat_type eq "CDS") {
	    
	    ## must populate the FeatureSequence and CDS_coords
	    
	    my $gi_number = $feature_ref->{"db_xref"};
	    my $translation = $feature_ref->{"translation"};
	    my $prot_acc = $feature_ref->{"protein_id"};
	    my $gene = $feature_ref->{"gene"};
	    my $product = $feature_ref->{"product"};
	    
	    &TransposonDB_API::load_FeatureSequence($dbproc, 
					     'feat_id' => $feat_id,
					     'seqType' => 'P',
					     'sequence' => $translation);
	    
	    my $coordListing = $feature_ref->{"coords"};
	    my @coordSets = split (/,/, $coordListing);
	    foreach my $coordSet (@coordSets) {
		my ($end5, $end3) = split (/-/, $coordSet);
		&TransposonDB_API::load_CDS_coords($dbproc,
						   'orf_ID' => $feat_id,
						   'end5' => $end5,
						   'end3' => $end3);
	    }
	}
	
    }
    @all_coords = sort {$a<=>$b} @all_coords;
    my $lend = shift @all_coords;
    my $rend = pop @all_coords;

    my @orients = sort {$orient_counter{$a}<=>$orient_counter{$b}} keys %orient_counter;
    my $orient = pop @orients;

    my ($end5, $end3) = ($orient eq "+") ? ($lend, $rend) : ($rend, $lend);

    ## update TE element:
    &TransposonDB_API::update_TransposonFeature($dbproc, 'feat_id' => $TE_id, end5 => $end5, end3 => $end3);

    print "Entry Loaded.\n";
    
}


####
sub check_exists_already {
    my ($gb_acc) = @_;
    my ($dbproc) = &connect_to_db("bhaas-lx","TransposonDB","access","access");
    my $query = "select count(*) from GenomeSeqRecord where genbank_acc = ?";
    my $count = &very_first_result_sql($dbproc, $query, $gb_acc);
    if ($count) {
	print "<font color='#ff0000'>Warning! Entry $gb_acc Exists in TransposonDB Already.</font><p>\n";
    } else {
	print "<p><font color='#0000ff'>Record for $gb_acc does not currently exist in TransposonDB.\n</font></p>";
    }
}

