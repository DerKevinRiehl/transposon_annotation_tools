echo "@@@"
echo "%PREFIX%"
echo "%RECIPE_DIR%"
echo "@@@"

mkdir "%PREFIX%\Scripts"
copy "%RECIPE_DIR%\seqkit.exe" "%PREFIX%\Scripts\seqkit.exe"
echo "copied 3 files"