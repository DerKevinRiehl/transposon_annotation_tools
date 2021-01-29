#!/bin/bash
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/bin/transposon_annotation_tools_mitefinderii
cp $RECIPE_DIR/bin/* $PREFIX/bin/transposon_annotation_tools_mitefinderii
cp $RECIPE_DIR/*.txt $PREFIX/bin/transposon_annotation_tools_mitefinderii
cp $RECIPE_DIR/miteFinderII $PREFIX/bin
chmod +x $PREFIX/bin/miteFinderII

