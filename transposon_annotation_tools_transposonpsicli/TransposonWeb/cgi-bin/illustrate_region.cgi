#!/local/perl/bin/perl

use lib ("PerlLib", ## common functions 
         "TransposonPerlLib"); ## custom functions for here.

use DBI;
use Mysql_connect;
use CGI;
use CGI::Carp qw(fatalsToBrowser); 
use strict;
use Data::Dumper;
use TransposonDB_API;
use GD;

$|++;

our $DEBUG = 0;

my $cgi = new CGI();
print $cgi->header(-type=> 'image/gif');

my %params = $cgi->Vars;

my $gi = $params{genbank_gi} or die "Error, need genbank_gi as parameter\n";
my $minRepLen = $params{MIN_REPEAT_LEN} || 12;
my $TE_id = $params{TE};

my ($dbproc) = &connect_to_db("bhaas-lx","TransposonDB","access","access");

my $repeat_file = "genomeSequences/$gi.seq.repfind";

my $start = $params{start} || 0;
my $end   = $params{end} || 0;
my $orig_end = $end;
my $delta = ($params{delta} =~ /\d+/) ? $params{delta} : 1000;


if ($TE_id) {
    my $query = "select end5, end3 from TransposonFeature where feat_id = ?";
    my $result = &first_result_sql($dbproc, $query, $TE_id);
    ($start, $end) = @$result;
}

if ($start && $end) {
    ($start, $end) = sort {$a<=>$b} ($start, $end);
    $start -= $delta; 
    $end += $delta;
} 
else {
    ($start, $end) = (1,0);
}

if ($start < 1) { 
    $start = 1;
    $end = $orig_end;
}

## get features of sequence:
my $query = "select tf.feat_id, tf.feat_type, tf.end5, tf.end3 from TransposonFeature tf, GenomeSeqRecord gs where tf.gseq_id = gs.transposon_seq_ID and gs.genbank_gi = ? order by tf.end5 desc";
my @results = &do_sql_2D($dbproc, $query, $gi);
my @features_to_draw;
foreach my $result (@results) {
    my ($feat_id, $feat_type, $end5, $end3) = @$result;
    push (@features_to_draw, [$feat_id, $feat_type, $end5, $end3]);
}

my $num_features_to_draw = scalar (@features_to_draw);

my $error;

unless(-s $repeat_file) {
    die "Error: Can't find the repeat listings at: $repeat_file.<br>";
}


####################  
# Parse repeat file.
####################
open (REP, "$repeat_file") or croak "Error: Can't open $repeat_file\n";
my @repeats;
my $sequence_length;
my %seen;

while (<REP>) {
    if (/^\#( \d+)?/) { 
	if ($1) {
	    $sequence_length = $1;
	    print STDERR "seqLen: $sequence_length\n" if $DEBUG;
	}
	next;
    }
    unless (/\w/) {next;}
    chomp;
    s/^\s+//; #trim leading ws
    
    my @x = split (/\s+/);
    
    my ($repLen, $ptA, $type, $repLenB, $ptB, $something, $pvalue) = split (/\s+/);
    
    if ($repLen < $minRepLen) { next;}

    my ($rep1_c1, $rep1_c2, $rep2_c1, $rep2_c2) = 
	($ptA, $ptA+$repLen-1, $ptB, $ptB+$repLenB-1);
    
    
    if ($rep1_c1 < $rep2_c2 && $rep1_c2 > $rep2_c1) {
	#overlap, avoid these (simple, or tandem)
	next;
    }
    
    ## avoid exact tandems:
    if ($rep2_c1 == $rep1_c2 + 1) { next;}
    
    #avoid dups.
    my @entries = sort (@x);
    my $key = join ("_", @entries);
    if ($seen{$key}) {
	next;
    } else {
	$seen{$key} = 1;
    }
    
    ## only examine repeats within specified limits.
    if ($end != 0) { ## make sure repeats are within the TE boundaries.
	unless (($rep1_c1 >= $start && $rep1_c1 <= $end) &&
		($rep2_c1 >= $start && $rep2_c1 <= $end) 
		) {
	    next;
	}
    }
    
    my $rep = { rep1_c1 => $rep1_c1,
		rep1_c2 => $rep1_c2,
		rep2_c1 => $rep2_c1,
		rep2_c2 => $rep2_c2,
		type => $type };
    push (@repeats, $rep);
}
close REP;

#####################
# Set sequence range.
#####################
unless ($end) {
    ### default setting.
    $end = $sequence_length;
}

#################
# Image drawing #
#################

## Set up the Scalable Image settings.
my $settings = { IMAGE_X_SIZE => 700,  #default image length
                 DRAW_PANEL_SCALER => 1.0, #percentage of image to draw matches, rest for text-annotation of match
                 ELEMENT_VERTICAL_SPACING => 10, #amount of vertical space consumed by each match
                 TICKER_TOGGLE => 1, #toggle for displaying a ticker for protein length, default is on.
	         ELEMENTS => [], #array ref for holding each element (as an array reference, also);
		 SEQ_START => $start, #crop viewing area by setting start-stop
	         SEQ_STOP => $end, #if not specified use max coordinate
		 DRAW_PANEL_SPACING => 50 #eliminates vertical white_space around image
		 }; 




my $image_x = $settings->{IMAGE_X_SIZE};
my $draw_panel_scaler = $settings->{DRAW_PANEL_SCALER};
my $element_vspacing = $settings->{ELEMENT_VERTICAL_SPACING};
my $draw_panel_size = $draw_panel_scaler * $image_x;
$settings->{DRAW_PANEL_SIZE} = $draw_panel_size;
my $text_panel_size = $image_x - $draw_panel_size;



## determine number of sections needed.
## Want 30 sections for repeat viewer
## want 4 for ticker

#my $number_of_sections = int ($image_height/$settings->{ELEMENT_VERTICAL_SPACING} +0.5);

my $total_number_sections = 40 + $num_features_to_draw;

my $repeat_arc_section = $total_number_sections - 30;
my $image_height = $total_number_sections * $settings->{ELEMENT_VERTICAL_SPACING};
my $ticker_section_number = $total_number_sections - 32;


## Calculate image_x coords for repeat features
foreach my $repeat (@repeats) {
    foreach my $type ("rep1_c1", "rep1_c2", "rep2_c1", "rep2_c2") {
	$repeat->{"${type}_tr"} = &coord_transform($settings, $repeat->{$type});
    }
}

## assign repeats to rows in the image:
my $max_row = &assign_repeats_to_rows(@repeats);
my ($im, $white, $blue, $black, $green, $red, $purple, $yellow); #image and colors declared.
my $found_repeats = 0;
my $imagefile;
my $repeat_list_ref;
my $orient_list_ref;


 main: {
     ## create GD image.
     $im = new GD::Image($image_x, $image_height);
     $white = $im->colorAllocate (255,255,255);
     $blue= $im->colorAllocate (0,0,255);
     $black = $im->colorAllocate (0,0,0);
     $green = $im->colorAllocate (0,255,0);
     $red = $im->colorAllocate (255,0,0);
     $purple = $im->colorAllocate(154,52,188);
     $yellow = $im->colorAllocate(255,255,0);

     &create_ticker ($settings, $im, $ticker_section_number);
     
     if (@repeats) {
	 &create_repeat_image($im, $repeat_arc_section);
     }
     

     if (@features_to_draw) {
	 my $section_start = 2; #skip one
	 my $section_end = $total_number_sections - $ticker_section_number;
	 &draw_seq_features ($im, $section_start, $section_end);
	 
	 
	 binmode STDOUT;
     }
     print $im->png();
     exit(0);
 }     

sub create_repeat_image {
    my $im = shift;
    my $section_num = shift;
    my $baseline = shift;
    my $height = 0.75 * $settings->{ELEMENT_VERTICAL_SPACING};
    
    my $baseline = $image_height - ($section_num * $settings->{ELEMENT_VERTICAL_SPACING});
    
    #draw line where arcs will end.
    $im->line(&coord_transform($settings, $start), $baseline, &coord_transform($settings, $end), $baseline, $black);
    
    my $y_coord = sub { my $row_num = shift;
			return ($image_height - ($section_num * $settings->{ELEMENT_VERTICAL_SPACING}) - ($row_num * $height));
		    };
    
    my @arcs;
    my @rects;
    foreach my $repeat (@repeats) {
	my $midpoint1 = ($repeat->{rep1_c1_tr} + $repeat->{rep2_c1_tr}) /2;
	my $width1 = abs (&coord_transform($settings, $repeat->{rep1_c1}) - &coord_transform($settings, $repeat->{rep2_c1})); 
	my $rownum1 = $repeat->{rep1_row};
	my $rep1_c1y = $y_coord->($rownum1);
	my $rep1_c2y = $y_coord->($rownum1) + $height;
	$repeat->{rep1_c1y} = $rep1_c1y;
	$repeat->{rep1_c2y} = $rep1_c2y;
	
	my ($arcColor) = ($repeat->{type} eq "F") ? $blue : $red;

	push (@rects, [($repeat->{rep1_c1_tr}, $rep1_c1y, $repeat->{rep1_c2_tr}, 
			$rep1_c2y, $green)]);
	push (@arcs, [$midpoint1,$baseline,
		      $width1,($baseline*2),
		      180,0,$arcColor]);
	
	
	my $midpoint2 = ($repeat->{rep1_c2_tr} + $repeat->{rep2_c2_tr}) /2;
	my $width2 = abs (&coord_transform($settings, $repeat->{rep1_c2}) - &coord_transform($settings, $repeat->{rep2_c2})); 
	my $rownum2 = $repeat->{rep2_row};
	my $rep2_c1y = $y_coord->($rownum2);
	my $rep2_c2y = $y_coord->($rownum2) + $height;
	$repeat->{rep2_c1y} = $rep2_c1y;
	$repeat->{rep2_c2y} = $rep2_c2y;
	
	push (@rects, [$repeat->{rep2_c1_tr}, $rep2_c1y, $repeat->{rep2_c2_tr}, 
		       $rep2_c2y, $green]);
	
	push (@arcs, [$midpoint2,$baseline,
		      $width2,($baseline*2),
		      180,0,$arcColor]);
    }
    
    ## draw the arcs first,
    foreach my $arc (@arcs) {
	$im->arc(@$arc);
    }
    ## now draw the rectangles
    foreach my $rect (@rects) {
	$im->filledRectangle(@$rect);
	pop @$rect;
	push (@$rect, $black);
	$im->rectangle(@$rect);
    }
}



sub coord_transform {
    my ($settings, $xcoord, $pure_flag) = @_;
    my ($min_element_length,$max_element_length, $draw_panel_size, $draw_panel_spacing) = ($settings->{SEQ_START},
											   $settings->{SEQ_STOP} - $settings->{SEQ_START}, 
											   $settings->{DRAW_PANEL_SIZE}, 
											   $settings->{DRAW_PANEL_SPACING});
    print "xcoord_in\t$xcoord\t" if $DEBUG;
    $xcoord = (($xcoord-$min_element_length)/$max_element_length * ($draw_panel_size - $draw_panel_spacing)) ;
    unless ($pure_flag) {
	$xcoord += $draw_panel_spacing / 2;
    }
    #$xcoord = int ($xcoord);
    print "xcoord_out\t$xcoord\n" if $DEBUG;
    return ($xcoord);
} 


sub create_ticker {
    my ($settings, $im, $ticker_section_number) = @_;
    my ($min_element_length,$max_element_length, $element_vspacing) = ($settings->{SEQ_START},$settings->{SEQ_STOP}, $settings->{ELEMENT_VERTICAL_SPACING});
    #print " ($min_element_length,$max_element_length, $element_vspacing) \n";
    my $curr_y = $image_height - $ticker_section_number * $element_vspacing;
    $im->line (&coord_transform($settings, $min_element_length), $curr_y, &coord_transform($settings, $max_element_length), $curr_y, $black);
    my $ticker_height_small = int (0.20 * $element_vspacing);
    my $ticker_height_large = int (0.50 * $element_vspacing);
    my $line_text_pointer = int (0.3 * $element_vspacing);
    my $value = int(($max_element_length-$min_element_length)/10);
    my $length = length($value);
    #print "value: $value\tLength: $length\n";
    my $lrg_interval = 10 ** ($length);
    #print "TICKER: $lrg_interval\n";
    if ($lrg_interval < 100) {$lrg_interval = 100;}
    my $sm_interval = int ($lrg_interval/10);
    my ($text_interval);
    if (int ((int ($max_element_length)-int($min_element_length))/$lrg_interval) <= 2) {
	$text_interval = 4 * $sm_interval;
    } else {
	$text_interval = $lrg_interval;
    }
    my($min_element_start) = $min_element_length - ($min_element_length % $sm_interval);
    my($max_element_start) = $max_element_length - ($max_element_length % $sm_interval) + $sm_interval;
    
    for (my $i = $min_element_start; $i <= $max_element_start; $i+= $sm_interval) {
	#print "$i\n";
	my $line_length = 0;
	if ($i%$sm_interval == 0) {
	    $line_length = $ticker_height_small;
	} 
	if ($i%$lrg_interval == 0) {
	    $line_length = $ticker_height_large;
	}
	if ( $i%$text_interval == 0) {
	    #add ticker text
	    my($label,$flabel);
	    $flabel = $i;
	    $im->string(gdSmallFont, &coord_transform ($settings, $i), ($curr_y + $element_vspacing), "$flabel", $black);
	    $im->line (&coord_transform ($settings, $i), $curr_y, &coord_transform($settings, $i), ($curr_y - $line_text_pointer), $black);
	}
	if ($line_length) {
	    $im->line (&coord_transform ($settings, $i), $curr_y, &coord_transform($settings, $i), ($curr_y + $line_length), $black);
	}
    }
    return ($curr_y + $element_vspacing);
}




sub assign_repeats_to_rows {
    my (@repeats) = @_;
    my @rows;
    my $max_row = 1;
    foreach my $repeat (@repeats) {
	foreach my $rep_type ("rep1", "rep2") {
	    my (@coords) = sort {$a<=>$b} ($repeat->{"${rep_type}_c1"}, $repeat->{"${rep_type}_c2"});
	    my $row_assignment;
	    my $current_row = 1;
	    while (!$row_assignment) {
		unless (ref $rows[$current_row]) {
		    $rows[$current_row] = [];
		}
		my @currently_placed_elements = @{$rows[$current_row]};
		my $row_ok=1;
		foreach my $element (@currently_placed_elements) {
		    my ($elem_c1, $elem_c2) = @$element;
		    if ($coords[0] <= $elem_c2 && $coords[1] >= $elem_c1) { #overlap
			$row_ok = 0;
			last;
		    }
		}
		if ($row_ok) {
		    $row_assignment = $current_row;
		    $repeat->{"${rep_type}_row"} = $row_assignment;
		    push (@{$rows[$row_assignment]}, [@coords]);
		} else {
		    $current_row++;
		}
	    }
	    if ($row_assignment > $max_row) { $max_row = $row_assignment;}
	}
    }
    return ($max_row);
}


####
sub draw_seq_features {
    my ($im, $section_start, $section_end) = @_;
    
    my $y_coord = $image_height - $section_start*$settings->{ELEMENT_VERTICAL_SPACING};
    
    my $feature_height = 0.6 * $settings->{ELEMENT_VERTICAL_SPACING};

    foreach my $feature (@features_to_draw) {
	
	my ($y1, $y2) = ($y_coord - $feature_height, $y_coord);
	my ($feat_id, $feat_type, $end5, $end3) = @$feature;
	my $color = $black;
	if ($feat_type eq "TE") {
	    $color = $purple;
	} elsif ($feat_type eq "CDS") {
	    $color = $yellow;
	}
	
	my ($lend_x, $rend_x) = sort {$a<=>$b} (&coord_transform($settings, $end5), &coord_transform($settings, $end3));
	
	my $annot_text = "";
	
	
	## determine left and right coordinates
	if ($feat_type eq "CDS") {
	    ## get the parts:
	    my $query = "select end5, end3 from CDS_coords where orf_ID = ?";
	    my @results = &do_sql_2D($dbproc, $query, $feat_id);
	    foreach my $result (@results) {
		my ($end5, $end3) = @$result;
		my ($lend_x, $rend_x) = sort {$a<=>$b} (&coord_transform($settings, $end5), &coord_transform($settings, $end3));
		$im->filledRectangle($lend_x, $y1, $rend_x, $y2, $yellow);
	    }
	    	    
	} elsif ($feat_type eq "TE") {
	    $im->filledRectangle($lend_x, $y1, $rend_x, $y2, $purple);
	}
	
	$im->rectangle($lend_x-1, $y1-1, $rend_x+1, $y2+1, $black); ## add outline
	
	my $query = "select qualifier, annotText from FeatureAnnots where feat_id = ?";
	my @results = &do_sql_2D($dbproc, $query, $feat_id);
	my %featAnnots;
	foreach my $result (@results) {
	    $featAnnots{$result->[0]} = $result->[1];
	}
	

	$annot_text = $featAnnots{product} . " " . $featAnnots{gene} . " " . $featAnnots{note};
	$annot_text =~ s/\s+/ /g;
	$annot_text =~ s/^\s//;
	
	$im->string(gdTinyFont, $lend_x, $y1-1, "$feat_type: $annot_text", $black);

	$y_coord -= $settings->{ELEMENT_VERTICAL_SPACING};
    }
}

