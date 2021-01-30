#!/bin/bash
MUSTPTH1=$(which makeblastdb)
MUSTPTH2=$(which blastn)
MUSTPTH3=$(which muscle)
MUSTPTH4=$(which stretcher)
MUSTPTH5=$(which cd-hit)
MUSTPTH6=$(which bedtools)

echo "DO INSTALLATION"
echo "perl $PREFIX/bin/SINE_Scan-v1.1.1/SINE_Scan_Installer.pl -d $PREFIX/bin/SINE_Scan-v1.1.1/ -a $PREFIX/bin/SINE_Scan-v1.1.1/sine_finder.py -f $MUSTPTH1 -b $MUSTPTH2 -M $MUSTPTH3 -e $MUSTPTH4 -c $MUSTPTH5 -S $PREFIX/bin/SINE_Scan-v1.1.1/SINEBase/SineDatabase.fasta -R $PREFIX/bin/SINE_Scan-v1.1.1/RNABase/RNAsbase.fasta -l $MUSTPTH6"
perl $PREFIX/bin/SINE_Scan-v1.1.1/SINE_Scan_Installer.pl -d $PREFIX/bin/SINE_Scan-v1.1.1/ -a $PREFIX/bin/SINE_Scan-v1.1.1/sine_finder.py -f $MUSTPTH1 -b $MUSTPTH2 -M $MUSTPTH3 -e $MUSTPTH4 -c $MUSTPTH5 -S $PREFIX/bin/SINE_Scan-v1.1.1/SINEBase/SineDatabase.fasta -R $PREFIX/bin/SINE_Scan-v1.1.1/RNABase/RNAsbase.fasta -l $MUSTPTH6

