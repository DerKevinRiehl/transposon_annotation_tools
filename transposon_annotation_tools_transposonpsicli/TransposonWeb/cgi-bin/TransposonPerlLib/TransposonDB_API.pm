package TransposonDB_API;

use Mysql_connect;
use strict;
use Data::Dumper;

our %supported_qualifiers;

BEGIN {
    my @qualifiers = qw (allele db_xref function gene label locus_tag note product protein_id 
			 pseudo standard_name rpt_family rpt_type rpt_unit transposon);
    foreach my $qual (@qualifiers) {
	$supported_qualifiers{$qual} = 1;
    }

}


####
sub load_GenomeSeqRecord {
    my ($dbproc, %keyvals) = @_;
    
    ## if entry already exists, return the identifier for it.
    my $query = "select transposon_seq_ID from GenomeSeqRecord where genbank_acc = ? and genbank_gi = ?";
    my $transposon_seq_ID = &very_first_result_sql($dbproc, $query, $keyvals{genbank_acc}, $keyvals{genbank_gi});
    unless ($transposon_seq_ID) {
	## insert new entry
	my $query = "insert GenomeSeqRecord (genbank_acc, genbank_gi, organism_name, seqTitle) values (?,?,?,?)";
	&RunMod($dbproc, $query, $keyvals{genbank_acc}, $keyvals{genbank_gi}, 
		$keyvals{organism_name}, $keyvals{seqTitle});
	
	$transposon_seq_ID = &get_last_insert_id($dbproc);
    }

    return ($transposon_seq_ID);
}


####
sub load_TransposonFeature {
    my ($dbproc, %keyvals) = @_;
    
    my $query = "insert TransposonFeature (gseq_id, feat_type, end5, end3) values (?,?,?,?)";
    &RunMod($dbproc, $query, $keyvals{gseq_id}, $keyvals{feat_type}, $keyvals{end5}, $keyvals{end3});
    my $id = &get_last_insert_id($dbproc);
    
    return ($id);
}


####
sub update_TransposonFeature {
    my ($dbproc, %keyvals) = @_;
    my $feat_id = $keyvals{feat_id};
    delete ($keyvals{feat_id});
    my @values;
    my $query = "update TransposonFeature set ";
    foreach my $key (keys %keyvals) {
	push (@values, $keyvals{$key});
	$query .= " $key = ?,";
    }
    chop $query;
    
    $query .= " where feat_id = $feat_id";
    
    &RunMod($dbproc, $query, @values);
}

sub load_FeatureSequence {
    my ($dbproc, %keyvals) = @_;

    ## accepted values for seqType:  
    #   P :protein
    #   G :genomic sequence
    #   C :cds

    my $query = "insert FeatureSequence (seqType, sequence, feat_id) values (?,?,?)";
    &RunMod($dbproc, $query, $keyvals{seqType}, $keyvals{sequence}, $keyvals{feat_id});
    
    my $fseq_id = &get_last_insert_id($dbproc);
    return ($fseq_id);
}



####
sub load_TransposonLink {
    my ($dbproc, %keyvals) = @_;
    
    my $query = "insert TransposonLink (element_ID, component_ID) values (?,?)";
    &RunMod($dbproc, $query, $keyvals{element_ID}, $keyvals{component_ID});
    my $link_id = &get_last_insert_id($dbproc);
    return ($link_id);
}

####
sub load_CDS_coords {
    my ($dbproc, %keyvals) = @_;
    
    my $query = "insert CDS_coords (orf_ID, end5, end3) values (?,?,?)";
    &RunMod($dbproc, $query, $keyvals{orf_ID}, $keyvals{end5}, $keyvals{end3});
    
    my $cds_id = &get_last_insert_id($dbproc);
    return ($cds_id);
}


#### 
sub load_FeatureAnnots {
    my ($dbproc, $feat_id, %keyvals) = @_;
    
    foreach my $qualifier (keys %keyvals) {
	if ($supported_qualifiers{$qualifier}) {
	    my $annotText = $keyvals{$qualifier};
	    my $query = "insert FeatureAnnots (feat_id, qualifier, annotText) values (?,?,?)";
	    &RunMod($dbproc, $query, $feat_id, $qualifier, $annotText);
	}
    }
}




1; #EOM
