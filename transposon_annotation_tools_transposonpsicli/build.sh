#!/bin/bash
mkdir -p $PREFIX/share
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/share/transposonPSIcli

cp -r $RECIPE_DIR/* $PREFIX/share/transposonPSIcli
cp $RECIPE_DIR/transposonPSI $PREFIX/bin
chmod +x $PREFIX/bin/transposonPSI

