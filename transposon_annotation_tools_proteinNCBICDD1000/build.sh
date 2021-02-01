#!/bin/bash
mkdir -p $PREFIX/share
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/share/proteinNCBICDD1000

cp -r $RECIPE_DIR/* $PREFIX/share/proteinNCBICDD1000
cp $RECIPE_DIR/proteinNCBICDD1000 $PREFIX/bin
chmod +x $PREFIX/bin/proteinNCBICDD1000

