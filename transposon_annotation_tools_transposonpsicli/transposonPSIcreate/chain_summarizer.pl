#!/usr/local/bin/perl


use lib ($ENV{EUK_MODULES}, $ENV{EGC_SCRIPTS});
use Egc_library;
use strict;
use DBI;
use Data::Dumper;
require "overlapping_nucs.ph";



my $usage = "usage: $0 chainsFile [swap_db_query]\n\n";
my $chainFile = $ARGV[0] or die $usage;
my $swap_db_query_flag = $ARGV[1];

my @chained_hits = &get_chained_hits($chainFile);

foreach my $chain (@chained_hits) {
	my ($hit_acc, $hit_lend, $hit_rend, $genome_acc, $genome_lend, $genome_rend, $pvalue, $orient) = ($chain->{hit_acc},
                                                                                                      $chain->{hit_lend},
                                                                                                      $chain->{hit_rend},
                                                                                                      $chain->{genome_acc},
                                                                                                      $chain->{genome_lend},
                                                                                                      $chain->{genome_rend},
                                                                                                      $chain->{pvalue},
                                                                                                      $chain->{orient}
                                                                                                      );
	if ($orient eq '-') {
        ($genome_lend, $genome_rend) = ($genome_rend, $genome_lend);
    }
    
    print "$genome_acc\t$genome_lend\t$genome_rend\t$hit_acc\t$hit_lend\t$hit_rend\t$pvalue\t$orient\n";
    
}

exit(0);


####
sub get_chained_hits {
    my $chainFile = shift;
    open (CHAINS, $chainFile) or die $!;
    
    my @chains;
    ## struct held in chains:
    #    lend, rend, orient, acc, pvalue
    #  pvalue is the smallest pvalue encountered for that chain.
    my $curr_chain;
    
    while (<CHAINS>) {
        if (/\w/) {
            my @x = split (/\t/);
            my ($hit, $genome_acc, $hit_lend, $hit_rend, $genome_end5, $genome_end3, $pvalue) = ($x[0], $x[5], $x[6], $x[7], $x[8], $x[9], $x[19]);
            
            if ($swap_db_query_flag) {
                ($hit, $hit_lend, $hit_rend, $genome_acc, $genome_end5, $genome_end3) = ($genome_acc, $genome_end5, $genome_end3, $hit, $hit_lend, $hit_rend);
            }
            

            my $orient = ($genome_end5 < $genome_end3) ? '+' : '-';
            my ($genome_lend, $genome_rend) = sort {$a<=>$b} ($genome_end5, $genome_end3);
            if ($curr_chain) {
                if ($curr_chain->{hit_acc} ne $hit) {
                    die "Error, $hit != " . Dumper ($curr_chain);
                }
                if ($curr_chain->{orient} ne $orient) {
                    die "Error, orient $orient != " . Dumper ($curr_chain);
                }
                if ($curr_chain->{genome_lend} > $genome_lend) {
                    $curr_chain->{genome_lend} = $genome_lend;
                }
                if ($curr_chain->{genome_rend} < $genome_rend) {
                    $curr_chain->{genome_rend} = $genome_rend;
                }
                if ($curr_chain->{hit_lend} > $hit_lend) {
                    $curr_chain->{hit_lend} = $hit_lend;
                }
                if ($curr_chain->{hit_rend} < $hit_rend) {
                    $curr_chain->{hit_rend} = $hit_rend;
                }
                
                if ($curr_chain->{pvalue} > $pvalue) {
                    $curr_chain->{pvalue} = $pvalue;
                }
            } else {
                $curr_chain = { 
                    hit_acc => $hit,
                    hit_lend => $hit_lend,
                    hit_rend => $hit_rend,
                    
                    genome_acc => $genome_acc,
                    genome_lend => $genome_lend,
                    genome_rend => $genome_rend,
                                        
                    pvalue => $pvalue,
                    orient => $orient,
                                        
                    };
                push (@chains, $curr_chain);
            }
        } else {
            # chain separator
            $curr_chain = undef;
        }
    }
    close CHAINS;
    return (@chains);
}


