mkdir -p $PREFIX\bin
mkdir -p $PREFIX\bin\transposon_annotation_tools_helitronscanner
mkdir -p $PREFIX\bin\transposon_annotation_tools_helitronscanner\helitronScannerRES

Xcopy /E /I $RECIPE_DIR\helitronscannerRES $PREFIX\bin\transposon_annotation_tools_helitronscanner\helitronScannerRES
copy $RECIPE_DIR\helitronscanner.bat $PREFIX\bin
copy $RECIPE_DIR\helitronscanner.py $PREFIX\bin

