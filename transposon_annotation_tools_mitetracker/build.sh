#!/bin/bash
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/mitetrackerLIB
cp -r $RECIPE_DIR/mitetrackerLIB/ $PREFIX/bin
cp $RECIPE_DIR/mitetracker $PREFIX/bin

chmod +x $PREFIX/bin/mitetracker

