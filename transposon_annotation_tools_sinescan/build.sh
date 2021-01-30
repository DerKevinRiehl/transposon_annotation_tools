#!/bin/bash
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/bin/SINE_Scan-v1.1.1
cp -r $RECIPE_DIR/SINE_Scan-v1.1.1 $PREFIX/bin
cp $RECIPE_DIR/sinescan $PREFIX/bin
chmod +x $PREFIX/bin/sinescan
