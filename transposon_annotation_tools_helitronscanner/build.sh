#!/bin/bash
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/bin/transposon_annotation_tools_helitronscanner
cp -r $RECIPE_DIR/helitronscannerRES $PREFIX/bin/transposon_annotation_tools_helitronscanner
cp $RECIPE_DIR/helitronscanner $PREFIX/bin
chmod +x $PREFIX/bin/helitronscanner
