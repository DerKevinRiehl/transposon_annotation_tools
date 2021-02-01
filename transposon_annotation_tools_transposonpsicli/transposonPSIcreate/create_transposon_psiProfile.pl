#!/usr/local/bin/perl

use strict;

my $transposon_db = "/usr/local/db/common/transposon_db.pep";
my $nr = "/usr/local/db/panda/AllGroup/AllGroup.niaa";

my $refSeq = $ARGV[0] or die "usage: $0 refSeq";
$refSeq =~ /^([^\.]+)\./;
my $coreSeqname = $1 or die "Sorry, couldn't extract core seqName from $refSeq\n";


## Create init profile:
my $cmd = "blastpgp -i $refSeq -d $nr -j 2 -C $coreSeqname.chkp";
system $cmd;
if ($?) {
    die "Error, cmd: $cmd\n died with ret($?)\n";
}

## Search profile against transposon_db for evaluation purposes:
my $cmd = "blastpgp -i $refSeq -d $transposon_db -R $coreSeqname.chkp -j 1 > $coreSeqname.psiblast.out";
system $cmd;
if ($?) {
    die "Error, cmd: $cmd\n died with ret($?)\n";
}

