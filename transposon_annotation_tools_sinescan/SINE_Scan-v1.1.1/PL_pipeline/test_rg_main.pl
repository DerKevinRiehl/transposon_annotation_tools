use strict;
use Bio::SimpleAlign;
use Bio::AlignIO;
use File::Basename;
if(@ARGV<6){
	print "$0 <TE seq> <genome> <min copy number of TE in genome> <output dir> <extend size> <end size for termial check>\n";
	print "Be sure that blastdb of genome is formated.\n";
	exit(0);
}

my $genome=$ARGV[1];
my $workdir=$ARGV[3];
my ($sineseq,$sinedir)=(basename($ARGV[0]),dirname($ARGV[0]));
my $minCp=$ARGV[2];
my $script=dirname($0);
my $sizeFlank=$ARGV[4];
my $sizeEnd=$ARGV[5];

print "Command\n$0 ".join(" ",@ARGV)."\n\n";
print "Create working dir: $workdir\n";


print "\nMake directory for each TE:\n";
my $filtersines=$ARGV[0];
system "perl ./PL_pipeline/makeDirforTE.pl $filtersines  $workdir";
open in,$filtersines or die "$!\n";
my %map=();
my %chr=();
my $order=0;
while(<in>){
	if($_=~/>/){
		$_=~s/>//;
		my @a=split(/\|/,$_);
		my $w="$a[0],$a[1]";
		$map{$w}=$order;
		$chr{$a[0]}.="$a[1],";
		$order++;
	}
}
close in;

open in,$genome or die "$!\n";
my %len=();
my $u;
while(<in>){
	chomp $_;
	if($_=~/>/){
		my @a=split(/\s+/,$_);
		$u=$a[0];
		$u=~s/>//;
		$len{$u}=0;
	}else{
		$len{$u}+=length($_);
	}
}
close in;

print "\nCheck copy number for TE clusters\n";
opendir DH,"$workdir" or die $!;
foreach my $name(sort {$a<=>$b} readdir DH){
	if(-d "$workdir/$name" && $name !~ /^\./){
###this is an iteration for an obvious boundary of Repeats#########
		my $convergence=0;####a condition to measure convergence of consensus sequence' length
		my $cycle=0;
		###first rank###
		my $CHR="";
		my $Position="";
		my $Direction='+';
		while(1){	
			print "\n-----work on cluster $name-----\n";
			my $path="$workdir/$name";
			print "Its path is: $path\n";
			print "TE sequence: $path/$name.sine.fa\n";
			print "Scan genome using a file above: $path/$name.sine.genome.bls\n";
			system "/home/maohlzj/ncbi-blast-2.2.31+/bin/blastn -task blastn -db $genome -query $path/$name.sine.fa -max_target_seqs 100000 -evalue 1e-10 -dust no -out $path/$name.sine.genome.bls -num_threads 1 -outfmt 6";
		
			open in,"$path/$name.sine.fa" or die "$!\n";
			my $query_length=0;
			while(<in>){
				chomp $_;
				if($_!~/^>/){
					$query_length+=length($_);
				}
			}
			close in;
#			system "rm $path/$name.sine.genome.bls";
			my %rank;
			my %Rank=();
			my $cnt=0;#estimated copy number of the candidate TE
			open fin,"<$path/$name.sine.genome.bls" or die $!;
			while(<fin>){
				my @x=split(/\t/,$_);
				my $dis=2*abs(5-$cycle);
				my $subject_length=abs($x[9]-$x[8])+1;
				if($x[9]-$x[8]<0){
					my $t=$x[8];
					$x[8]=$x[9];
					$x[9]=$t;
				}
				if($cycle == 0 || $x[2]*$x[3]>=80*$query_length && $subject_length>=0.8*$query_length && $subject_length<=$query_length/0.8){
					if($cycle>0){
						if(abs($x[6]-1)>3 || abs($x[7]-$query_length)>3){
							next;
						}
					}
					my $number=$x[2]*$x[3]/($x[9]-$x[8]+1);
					my $iden=$x[2]*$x[3]/100;
					$rank{$iden}{$_}=$number;
					++$cnt;
				}
			}
			close fin;
				
			print "In cycle $cycle,Find $cnt rough good hits with identical percentage more than 80%\n";
			if($cnt<$minCp){
				if($cycle == 0){
					print "Cluster $name has $cnt good hits in genome, less than the cutoff value $minCp: stop analyse this cluster!\n";
				}else{
					print "Cluster $name has $cycle round test and this round it has $cnt good hits in genome, less than the cutoff value $minCp and its boundaries perhaps is to be stable: stop analyse this cluster!\n";
				}
#				unlink glob "$path/* $path/.*";
#				rmdir $path;
				#system "rm $path/$name.sine.genome.* ";
				last;
			}
	#		system "rm $path/$name.sine.genome.solar";
			my $p=0;
			my $cc=0;
			my $flag=1;
			my @hits=();
			my @start=();
			my @stop=();
			foreach my $identity(reverse sort {$a<=>$b} keys %rank){###identical percentage of query
				if($flag==0){
					last;
				}
				foreach my $line(reverse sort {$rank{$identity}{$a}<=>$rank{$identity}{$b}} keys%{$rank{$identity}}){ ###identical percentage of subject
					if($p==35){###35?###
						$flag=0;
						last;
					}else{
						++$p;
						my @a=split(/\t/,$line);
						if($p == 1){
							$CHR=$a[1];
							$Position="$a[8]-$a[9]";
							if($a[8] > $a[9]){
								$Direction='-';
							}
						}
						push @hits,$line; 
						push @start,$a[6];
						push @stop,$a[7];
						++$cc;
					}
				}	
			}
			
			my $n=@hits;
			@start=sort{$a<=>$b} @start;
			@stop=sort{$a<=>$b} @stop;
			my $median_start=0;
			my $median_stop=0;
			if($n%2 == 0){
				$median_start=int(($start[$n/2-1]+$start[$n/2])/2);
				$median_stop=int(($stop[$n/2-1]+$stop[$n/2])/2);
			}else{
				$median_start=$start[$n/2];
				$median_stop=$stop[$n/2];
			}	
			my $deta_start=$median_start;
			my $deta_stop=abs($query_length-$median_stop)+1;
			my $preflank=$sizeFlank;
			my $suffixflank=$sizeFlank;
			if($deta_start > 10){
				$preflank=0;
			}
			if($deta_stop > 10){
				$suffixflank=0;
			}
			my @last_input=();
			foreach my $stronghit (@hits){
				my @x=split(/\s+/,$stronghit);
				if($x[9]-$x[8]<0){
					my $t=$x[8];
					$x[8]=$x[9];
					$x[9]=$t;
				}
				if($x[8] > $preflank and $x[9]+$suffixflank<$len{$x[1]}){
					push @last_input,$stronghit;
				}
			}
			$cnt=@last_input;
			print "Find $cnt  accurate good hits with 5 side preflank_length=$preflank bp sequence and 3 side suffixflank_length=$suffixflank bp sequence\n";
			if($cnt<$minCp){
				if($cycle == 0){
					print "Cluster $name has $cnt good hits in genome, less than the cutoff value $minCp: stop analyse this cluster!\n";
				}else{
					print "Cluster $name has $cycle round test and this round it has $cnt good hits in genome, less than the cutoff value $minCp and its boundaries perhaps is to be stable: stop analyse this cluster!\n";
				}
#				unlink glob "$path/* $path/.*";
#				rmdir $path;
				system "rm $path/$name.sine.genome.* ";
				last;
			}
			open fout,">$path/$name.sine.genome.filter" or die $!;
			foreach my $i (@last_input){
				print fout"$i";
			}	
			close fout;
			
			print "Extract $cc best hits: $path/$name.sine.extendseq\n";
			system "perl ./PL_pipeline/extendseq.pl $path/$name.sine.genome.filter $genome $preflank $suffixflank >$path/$name.sine.extendseq";
			if($cycle == 0){
				system "cat $path/$name.sine.extendseq >$path/$name.sine.origin.fa";
			}
			$cycle++;
#			system "rm $path/$name.sine.genome.solar.filter";

#####check end then output a good longer (with ends) seed sequence for next-round-searching ############
###step first: checkend #######
			my $New_name="$CHR,$Position,$Direction";
			my $scores=`perl ./PL_pipeline/test_fix.pl $path/$name.sine.extendseq $New_name`;
#			my @Position=split(/\s+/,$scores);
#			print "$Position[0]\t$Position[1]\n";
#			$scores=abs($Position[1]-$Position[0])+1;
			if($scores == 0){
				last;
			}
			my $deta=abs($scores-$convergence);
			$convergence=$scores;	
			if($deta < 4){
				my $BD=`perl ./PL_pipeline/CheckEnd.pl $path/$name.sine.extendseq 60 25`;
				print "$BD\n";
				#######Blast-Filter###########
				open fin,"<$path/$name.sine.genome.bls" or die $!;
				while(<fin>){
					my @x=split(/\t/,$_);
					my $dis=2*abs(5-$cycle);
					my $subject_length=abs($x[9]-$x[8])+1;
					if($x[9]-$x[8]<0){
						my $t=$x[8];
						$x[8]=$x[9];
						$x[9]=$t;
					}
					if($x[2]*$x[3]>=80*$query_length && $subject_length>=0.8*$query_length && $subject_length<=$query_length/0.8){
						if(abs($x[6]-1)<=3 && abs($x[7]-$query_length)<=3){
							my @a=split(/,/,$chr{$x[1]});
							foreach my $j (@a){
								my @A=split(/-/,$j);
								my $R=($x[8]-$A[1])*($x[9]-$A[0]);
								if($R < 0){
									my $B="$x[1],$j";
									my $DIR="$workdir/$map{$B}";
									if(-d $DIR && $map{$B} != $name){
										system "rm -rf $DIR";
										$chr{$x[1]}=~s/$j,//;
										last;
									}
								}		
							}		
						}
					}
				}
				close fin;
				print "The boundaries are stable and jump out the loop.\n";
				print "Output obvious boundaries in $path/$name.boundary.fa\n";	
				print "Give Repeat candidates with obvious ends\n";
				system "rm $path/$name.sine.genome.* ";
				last;
			}
			if($cycle > 3){
				my $BD=`perl ./PL_pipeline/CheckEnd.pl $path/$name.sine.extendseq 60 25`;
				print "$BD\n";
				print "The number of iteration is large enough and jump out the loop.\n";
				if($deta > 10){
					print "this cluster's boundaries perhaps are error, check it manually\n";
				}
				print "Output obvious boundaries in $path/$name.boundary.fa\n";	
				print "Give Repeat candidates with obvious ends\n";
				system "rm $path/$name.sine.genome.* ";
				last;
			}
###original seeds to genome_scan and use boundary to pairwise alignment,check its quality to determine whether it should be retained######	
		}
	}
}
system "perl ./PL_pipeline/clusterSeqs.pl $workdir/passList $filtersines";
closedir DH;
print "\nDone.\n";

