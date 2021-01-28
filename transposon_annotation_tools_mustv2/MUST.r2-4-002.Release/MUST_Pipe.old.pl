#!/usr/bin/perl -w
use strict;
use MIME::Base64;
use Bio::SeqIO;
use Bio::SearchIO;
my$tline="";my@tlist=();my$tempi=0;my$tempj=0;my$tempk=0;my$pline="";my@plist=();my$tid=0;my%tidx=();
my$egDirTemp="tempData";my$egScriptMUST="MUST";my$egScriptPWD="PairWiseDistance";my$egScriptFilter="filterRedundancy.new.pl";my$egLogFilter="log.filterRedundancy.new.log";
if(system("echo \"\" > \"$egLogFilter\"")){die"-Error when clearing the log file!\n";}my$egDirScript=$0;
@tlist=split(/\//,$egDirScript);
if(@tlist==1){$egDirScript=".";}else{$egDirScript=~s/\/$tlist[@tlist-1]$//g;}$|=1;
my$gMinTIR=8;
my$gMaxTIR=50;
my$gMinDR=2;
my$gMaxDR=30;
my$gMinMITE=100;
my$gMaxMITE=600;
my$gFixedFlanking=50;
my$gMutationRate=0.80;
my$egFileSeq="";
my$egFileOut="";
if((@ARGV==3)or(@ARGV==11)){}else{print"Error in syntax!\n";
&efSyntax();
die"-----------------------------------------------------------------------------\n";}$egFileSeq=$ARGV[0];$egFileOut=$ARGV[1];$egDirTemp=$ARGV[2];
if(-d$egDirTemp){}else{if(system("mkdir -p \"$egDirTemp\"")){die"-Error when creating the temporary directory!\n";}}if(@ARGV==11){$gMinTIR=$ARGV[3];$gMaxTIR=$ARGV[4];$gMinDR=$ARGV[5];$gMaxDR=$ARGV[6];$gMinMITE=$ARGV[7];$gMaxMITE=$ARGV[8];$gFixedFlanking=$ARGV[9];$gMutationRate=$ARGV[10];}my$egCutOffIdent=$gMutationRate;
my%mMITE2ID=();my@mID=();my@mGID=();my@mStart=();my@mEnd=();my@mStrand=();my@mLength=();my@mDR=();my@mDRIdent=();my@mDRLeft=();my@mDRRight=();my@mTIR=();my@mTIRIdent=();my@mTIRLeft=();my@mTIRRight=();my@mDRAT=();my@mTIRAT=();my@mMITEAT=();my@mFDRAT=();my@mPreSite=();my@mPostSite=();my@mMITE=();my@mCopyNum=();my@mCopies=();my@mTotalScore=();my@mFlag=();my$mNum=0;
my%mData=();
my$egCutOffCopyScore=10;
my$egLineLength=60;
my$egTimeOld=`date +%s`;
my$egTimeNew=`date +%s`;
my$egSeq="";
my$egID="";
my$lFile="";
my$lFSeq="";
my$lCmdLine="";
my$lFileDat="";
my$lFileTIR="";
my$lFileEdge="";
print"Scanning the nucleotide sequences for potential MITEs ... \n";
my$egSeqObj=new Bio::SeqIO(-file=>"$egFileSeq",-format=>"fasta");
while(my$egSeqO=$egSeqObj->next_seq){$tline=$egSeqO->seq;
$egSeq="\U$tline\E";
$egID=$egSeqO->id;
$lFile="";
$lFSeq="";
$lCmdLine="";
$lFileDat="";
$lFileTIR="";
print"Processing [$egID] ... ";
my$lFSeq=&efFormatSeq($egSeq,$egLineLength,1);
$lFSeq=~s/ +//g;
print" [FSeq:".length($egSeq).":".length($lFSeq)."]";
$lFile="$egDirTemp/temp_seq.raw";
open(lOut,">$lFile")or die"Error when saving the temporary sequence!\n";
print lOut"$lFSeq\n";
close(lOut);
print" [Save:FSeq]";
$lFileDat="$egDirTemp/temp_seq.dat";
$lFileTIR="$egDirTemp/temp_tir.dat";
$lCmdLine="\"$egDirScript\/$egScriptMUST\" \"$lFile\" \"$lFileDat.before_filter\" \"$lFileTIR\" $gMinTIR $gMaxTIR $gMinDR $gMaxDR $gMinMITE $gMaxMITE $gFixedFlanking $gMutationRate";
if(system("$lCmdLine >/dev/null")){die"-Error occurred when running the script \"$egDirScript\/$egScriptMUST\"!\n";}print" [MUST]";
$lCmdLine="\"$egDirScript/$egScriptFilter\" \"$lFileDat.before_filter\" \"$lFileDat\"";
if(system("$lCmdLine >> $egLogFilter")){die"-Error occurred when running the script \"$egDirScript\/$egScriptFilter\"!\n";}print" [Filter]";
$tempi=&efLoadMUST($lFileDat);
print" [MITE:$tempi]";
$egTimeNew=`date +%s`;
print" [Time:".($egTimeNew-$egTimeOld)." seconds]\n";
$egTimeOld=$egTimeNew;}my$egCmdLine="";
print"Removing redundancy in the predicted MITEs ... ";
for($tempi=0;$tempi<$mNum;$tempi++){$mTIR[$tempi]=int($mTIR[$tempi]);
$mDR[$tempi]=int($mDR[$tempi]);}$tempk=0;
for($tempi=0;$tempi<$mNum;$tempi++){if($mFlag[$tempi]==1){for($tempj=$tempi+1;$tempj<$mNum;$tempj++){if(($tempj<$mNum)and($mFlag[$tempj]==1)){if(($mGID[$tempi]eq$mGID[$tempj])and(&efIsOverlap($mStart[$tempi],$mEnd[$tempi],$mStart[$tempj],$mEnd[$tempj])==1)){if($mTIR[$tempi]>$mTIR[$tempj]){$mFlag[$tempj]=0;
$tempk++;}elsif($mTIR[$tempi]<$mTIR[$tempj]){$mFlag[$tempi]=0;
$tempk++;
last;}else{if($mDR[$tempi]>$mDR[$tempj]){$mFlag[$tempj]=0;
$tempk++;}elsif($mDR[$tempi]<$mDR[$tempj]){$mFlag[$tempi]=0;
$tempk++;
last;}else{$mFlag[$tempj]=0;
$tempk++;}}}}}}else{}}print" [RuleBasedRemoving:$tempk]";
$tempj=0;
for($tempi=0;$tempi<$mNum;$tempi++){if($mFlag[$tempi]==1){$tempj++;}}$egTimeNew=`date +%s`;
print" [MITE:$tempj] [done] [Time:".($egTimeNew-$egTimeOld)." seconds]\n";
$egTimeOld=$egTimeNew;
print"Clustering ... ";
my@mClusterID=();
my%mC2IDs=();
for($tempi=0;$tempi<$mNum;$tempi++){$mClusterID[$tempi]=-1;}my$mCluster=0;
my$mMITEWithCluster=0;
my%tInCluster=();
my$tFirstMITE=-1;
my@tInClusterList=();
my@tNeighborList=();
my$egClusterMode=1;
my$egFileMITESeq="$egDirTemp/temp-mite-seq.fasta";
open(efO,">$egFileMITESeq")or die"E!\n";
for($tempi=0;$tempi<$mNum;$tempi++){print efO">$tempi\n$mMITE[$tempi]\n";}close(efO);
system("formatdb -p F -i \"$egFileMITESeq\"");
system("megablast -F F -d \"$egFileMITESeq\" -i \"$egFileMITESeq\" -m 9 -o \"$egDirTemp/temp-mite-seq.all-vs-all.megablast.txt\"");
print" [BLAST]";
my%idxIdent=();
my@mCC=();
my@mFF=();
open(efI,"$egDirTemp/temp-mite-seq.all-vs-all.megablast.txt")or die"EE!\n---[$egDirTemp/temp-mite-seq.all-vs-all.megablast.txt]---\n";
while(<efI>){$tline=$_;$tline=~s/[\r\n]//g;@tlist=split(/\t/,$tline);
if($tline=~/^#/){}elsif(@tlist>=12){$idxIdent{"$tlist[0]|vs|$tlist[1]"}=$tlist[2]*$tlist[3]/100;
$idxIdent{"$tlist[1]|vs|$tlist[0]"}=$tlist[2]*$tlist[3]/100;}}close(efI);
for($tempi=0;$tempi<$mNum;$tempi++){$mCC[$tempi]="";$mFF[$tempi]=0;
$mClusterID[$tempi]=0;}for($tempi=0;$tempi<$mNum;$tempi++){for($tempj=$tempi+1;$tempj<$mNum;$tempj++){$tid=$idxIdent{"$tempi|vs|$tempj"};
if(!(defined$tid)){$tid=$idxIdent{"$tempj|vs|$tempi"};}if(defined$tid){if(($tid>=$egCutOffIdent*$mLength[$tempi])and($tid>=$egCutOffIdent*$mLength[$tempj])){$mCC[$tempi].=" $tempj";
$mCC[$tempj].=" $tempi";}}}}my$egClusterNumber=1;
for($tempi=0;$tempi<$mNum;$tempi++){my$kCID=$egClusterNumber;
$mCC[$tempi]=~s/^\s+//g;
$mCC[$tempi]=~s/\s+$//g;
@tlist=split(/\s+/,$mCC[$tempi]);
for($tempj=0;$tempj<@tlist;$tempj++){if(($mClusterID[$tlist[$tempj]]>0)and($mClusterID[$tlist[$tempj]]<$kCID)){$kCID=$mClusterID[$tempj];}}for($tempj=0;$tempj<@tlist;$tempj++){if($mClusterID[$tlist[$tempj]]==0){$mClusterID[$tlist[$tempj]]=$kCID;}}if($mClusterID[$tempi]==0){$mClusterID[$tempi]=$kCID;}if($kCID==$egClusterNumber){$egClusterNumber++;}}$mCluster=$egClusterNumber-1;
for($tempi=0;$tempi<$mNum;$tempi++){$tid=$mC2IDs{"$mClusterID[$tempi]"};
if(defined$tid){$tid.=";$tempi";}else{$tid="$tempi";}$mC2IDs{"$mClusterID[$tempi]"}=$tid;}$egTimeNew=`date +%s`;
print" [done] [ClusterNum:$mCluster] [Time:".($egTimeNew-$egTimeOld)." seconds]\n";
$egTimeOld=$egTimeNew;
print"Deciding the strand information of each valid MITE ... ";
for($tempi=0;$tempi<$mNum;$tempi++){if(($mFlag[$tempi]==1)and($mStrand[$tempi]eq"")){$mStrand[$tempi]="+";
for($tempj=$tempi+1;$tempj<$mNum;$tempj++){if(($mFlag[$tempj]==1)and($mClusterID[$tempj]==$mClusterID[$tempi])and($mStrand[$tempj]eq"")){if(&efIsSameStrand($mMITE[$tempi],$mMITE[$tempj])==1){$mStrand[$tempj]="+";}else{$mStrand[$tempj]="-";}}}}}$egTimeNew=`date +%s`;
print" [done] [Time:".($egTimeNew-$egTimeOld)." seconds]\n";
$egTimeOld=$egTimeNew;
print"Saving the data ... ";
$tempj=0;
open(efOut,">$egFileOut\.before_filter")or die"Error when saving the predicted MITEs!\n";
print efOut"#GID	ID	Cluster	Start	End	Strand	Length	DR	DRIdent	DRLeft	DRRight	TIR	TIRIdent	TIRLeft	TIRRight	DRAT	TIRAT	MITEAT	FDRAT	PreSite	PostSite	MITE	TotalScore\n";
my%aCID=();
for($tempi=0;$tempi<$mNum;$tempi++){if(($mClusterID[$tempi]>0)and($mFlag[$tempi]==1)){$mCopies[$tempi]=~s/\s/,/g;
print efOut"$mGID[$tempi]	$mID[$tempi]	$mClusterID[$tempi]	$mStart[$tempi]	$mEnd[$tempi]	$mStrand[$tempi]	$mLength[$tempi]	$mDR[$tempi]	$mDRIdent[$tempi]	$mDRLeft[$tempi]	$mDRRight[$tempi]	$mTIR[$tempi]	$mTIRIdent[$tempi]	$mTIRLeft[$tempi]	$mTIRRight[$tempi]	$mDRAT[$tempi]	$mTIRAT[$tempi]	$mMITEAT[$tempi]	$mFDRAT[$tempi]	$mPreSite[$tempi]	$mPostSite[$tempi]	$mMITE[$tempi]	$mTotalScore[$tempi]\n";
$tempj++;$aCID{"$mClusterID[$tempi]"}=1;}}close(efOut);
my$egCutOffCopyNum=3;
if(system("$egDirScript/filterCopyNum.pl \"$egFileOut.before_filter\" \"$egFileOut.filter-1\" $egCutOffCopyNum >/dev/null")or system("$egDirScript/filterSegDup-2.pl \"$egFileOut.filter-1\" \"$egFileOut\" 4 3 >/dev/null")or system("rm *.log")){die"-Error when filtering!\n";}$egTimeNew=`date +%s`;
print" [done] [Valid MITE Copy:$tempj] [MITEs:".(scalar keys%aCID)."] [Time:".($egTimeNew-$egTimeOld)." seconds]\n";
$egTimeOld=$egTimeNew;
my$qMITENum=0;
my$qClusterNum=0;
%tidx=();
open(efI,"$egFileOut")or die"E!\n";
while(<efI>){$tline=$_;$tline=~s/[\r\n]//g;@tlist=split(/\t/,$tline);
if($tline=~/^#/){}elsif(@tlist>=5){$qMITENum++;
$tidx{"$tlist[2]"}=1;}}$qClusterNum=scalar keys%tidx;
close(efI);
print"There are $qMITENum MITEs detected in $qClusterNum clusters.\n";
$|=0;
sub efLoadMUST{my($tFileDat)=@_;
%mData=();
my$tStartID=$mNum;
my$tNowNum=0;
open(tIn,"$tFileDat")or die"Error when loading the Dat!\n";
while(<tIn>){$tline=$_;
$tline=~s/[\r\n]//g;
if($tline=~/^\/\/$/){$mData{"SEQID"}=$egID;
my@ttlist=();
my$tNowID=&efValidString($mData{"ID"});
if($tNowID=~/^\d+$/){$mMITE2ID{"$egID $tNowID"}=$mNum;
$mID[$mNum]=&efValidString($tNowID);
$mGID[$mNum]=&efValidString($egID);
@ttlist=split(/-/,&efValidString($mData{"POSITION"}),2);
$mStart[$mNum]=&efValidString($ttlist[0]);
$mEnd[$mNum]=&efValidString($ttlist[1]);
$mStrand[$mNum]="";
$mLength[$mNum]=abs($mEnd[$mNum]-$mStart[$mNum])+1;
$mDR[$mNum]=&efValidString($mData{"DR"});
$mDRIdent[$mNum]=&efValidString($mData{"DRIDENT"});
$mDRLeft[$mNum]=&efValidString($mData{"DRLEFT"});
$mDRRight[$mNum]=&efValidString($mData{"DRRIGHT"});
$mTIR[$mNum]=&efValidString($mData{"TIR"});
$mTIRIdent[$mNum]=&efValidString($mData{"TIRIDENT"});
$mTIRLeft[$mNum]=&efValidString($mData{"TIRLEFT"});
$mTIRRight[$mNum]=&efValidString($mData{"TIRRIGHT"});
$mDRAT[$mNum]=&efValidString($mData{"DRAT"});
$mTIRAT[$mNum]=&efValidString($mData{"TIRAT"});
$mMITEAT[$mNum]=&efValidString($mData{"MITEAT"});
$mFDRAT[$mNum]=&efValidString($mData{"FDRAT"});
$mPreSite[$mNum]=&efValidString($mData{"PRESITE"});
$mPostSite[$mNum]=&efValidString($mData{"POSTSITE"});
$mMITE[$mNum]=&efValidString($mData{"MITE"});
$mTotalScore[$mNum]=0;
$mCopyNum[$mNum]=1;
$mCopies[$mNum]="$mNum";
$mFlag[$mNum]=1;
$tNowNum++;
$mNum++;}}elsif($tline=~/^([a-zA-Z]+)\t(.*)$/){my$tID=&efValidString($1);
my$tDat=&efValidString($2);
if(($tID ne"")and($tDat ne"")){$mData{"$tID"}=$tDat;}}}close(tIn);
return$tNowNum;}sub efValidString{my($tStr)=@_;
if((defined$tStr)and($tStr ne"")){return$tStr;}else{return"";}}sub efSyntax{my$egSyntax=decode_base64("TUlURSBVbmNvdmVyaW5nIFN5c1RlbSAgICB2ZXJzaW9uIDEuMAogICAgICAgICAgICAgICAgIEZl
bmdmZW5nIFpob3UgKGMpIDIwMDgtMDYtMjAKLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0t
LS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLQogICAu
L01VU1RfUGlwZS5wbCA8aW5wdXQuZmFzdGE+IDxvdXRwdXQuTUlURS5kYXQ+IDxEaXJUZW1wPiBb
b3B0aW9uc10KTm90ZToKICAgSWYgYW55IG9wdGlvbnMgYXJlIHByb3ZpZGVkLCB0aGUgdXNlciBo
YXZlIHRvIGlucHV0IHRoZSB2YWx1ZXMgZm9yIGFsbCB0aGUgb3B0aW9ucy4KT3B0aW9uIGxpc3Q6
CiAgIDxNSU5fVElSX2xlbmd0aD4gICAgICAgICBbOF0KICAgPE1BWF9USVJfbGVuZ3RoPiAgICAg
ICAgIFs1MF0KICAgPE1JTl9EUl9sZW5ndGg+ICAgICAgICAgIFsyXQogICA8TUFYX0RSX2xlbmd0
aD4gICAgICAgICAgWzMwXQogICA8TUlOX01JVEVfbGVuZ3RoPiAgICAgICAgWzEwMF0KICAgPE1B
WF9NSVRFX2xlbmd0aD4gICAgICAgIFs2MDBdCiAgIDxGSVhFRF9GTEFOS0lOR19sZW5ndGg+ICBb
NTBdCiAgIDxNdXRhdGlvbl9SYXRlPiAgICAgICAgICBbMC43XQo="
);
$egSyntax="MITE Uncovering SysTem    version 2.2.001
---------------------(c) Fengfeng Zhou (FengfengZhou\@gmail.com) 2013-11-20
---------------------------------------------------------------------------------
./MUST_Pipe.pl <input.fasta> <output.MITE.dat> <DirTemp> [options]
Note:
If any options are provided, the user have to input the values for all the options.
Option list:
-<MIN_TIR_length>         [8]
-<MAX_TIR_length>         [50]
-<MIN_DR_length>          [2]
-<MAX_DR_length>          [30]
-<MIN_MITE_length>        [100]
-<MAX_MITE_length>        [600]
-<FIXED_FLANKING_length>  [50]
-<Mutation_Rate>          [0.80]
-----------------------------------------------------------------------------
";
print$egSyntax;}sub efFormatSeq{my($fSeq,$tBlock,$tNumber)=@_;
my$tPos=0;
my$tRSeq="";
while(1){my$tStart=$tBlock*$tPos;
my$tEnd=$tBlock*($tPos+1)-1;
if($tStart>=length($fSeq)){last;}if($tEnd>=length($fSeq)){$tEnd=length($fSeq)-1;}$tRSeq=$tRSeq.substr($fSeq,$tStart,abs($tEnd-$tStart)+1)." ";
if(($tPos>=0)and($tPos%$tNumber==0)){$tRSeq=$tRSeq."\n     ";}$tPos++;}$tRSeq=~s/^ +//g;
$tRSeq=~s/ +$//g;
$tRSeq="     ".$tRSeq;
return$tRSeq;}sub efIsOverlap{my($s1,$e1,$s2,$e2)=@_;
my$t=0;
if($s1>$e1){$t=$s1;
$s1=$e1;
$e1=$t;}if($s2>$e2){$t=$s2;
$s2=$e2;
$e2=$t;}if((($s2>=$s1)and($s2<=$e1))or(($s1>=$s2)and($s1<=$e2))){return 1;}else{return 0;}}sub efIsSameStrand{my($iSeq1,$iSeq2)=@_;
if((length($iSeq1)>0)and(length($iSeq2)>0)){my$iF1="$egDirTemp/seq1.fasta";
open(iO,">$iF1")or die"Error when saving sequence 1\n";
print iO ">Query\n$iSeq1\n";
close(iO);
my$iF2="$egDirTemp/seq2.fasta";
open(iO,">$iF2")or die"Error when saving sequence 2\n";
print iO">DB\n$iSeq2\n";
close(iO);
my$iFO="$egDirTemp/seq.bl2seq.dat";
my$iCmdLine="bl2seq -F F -p blastn -i \"$iF1\" -j \"$iF2\" -o $iFO";
if(system($iCmdLine)){print"-Error when running bl2seq!\n";
return-1;}else{my$iSearchIO=new Bio::SearchIO(-file=>"$iFO",-format=>"blast");
my$iFlag=1;
while(my$iResult=$iSearchIO->next_result){my$iNameQ=$iResult->query_name();
while(my$iHit=$iResult->next_hit){my$iNameH=$iHit->name();
while(my$iHSP=$iHit->next_hsp){my$iStrandQ=$iHSP->strand('query');
my$iStrandH=$iHSP->strand('hit');
if(($iStrandQ eq$iStrandH)and($iStrandQ ne"")){$iFlag=1;}else{$iFlag=0;}goto lbGotIt;}}}lbGotIt:return$iFlag;}}else{return-1;}}sub efPrintData{my(%pIdx)=@_;
my@pKeys=keys%pIdx;
for(my$pi=0;$pi<@pKeys;$pi++){print"--[$pKeys[$pi]] => [".$pIdx{"$pKeys[$pi]"}."]---\n";}print"----------------------------------------------------------------------\n";}
