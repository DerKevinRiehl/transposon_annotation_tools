#!/bin/bash
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/bin/MUST.r2-4-002.Release
cp -r $RECIPE_DIR/MUST.r2-4-002.Release $PREFIX/bin
cp $RECIPE_DIR/mustv2 $PREFIX/bin
chmod +x $PREFIX/bin/mustv2

#conda install -y -c bioconda/label/cf201901 blat
#conda install -y -c bioconda/label/cf201901 blast
#conda install -y -c conda-forge/label/cf202003 perl==5.26.2
#conda install -y -c bioconda/label/cf201901 perl-bioperl-core==1.007002
