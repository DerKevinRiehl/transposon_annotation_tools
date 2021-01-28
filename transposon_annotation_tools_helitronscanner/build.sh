#!/bin/bash
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/share/transposon_annotation_tools_helitronscanner
cp $RECIPE_DIR/*.jar $PREFIX/share/transposon_annotation_tools_helitronscanner
cp $RECIPE_DIR/*.lcvs $PREFIX/share/transposon_annotation_tools_helitronscanner
cp $RECIPE_DIR/helitronscanner $PREFIX/bin
chmod +x $PREFIX/bin/helitronscanner

