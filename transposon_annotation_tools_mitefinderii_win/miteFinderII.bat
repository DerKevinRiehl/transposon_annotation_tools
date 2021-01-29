SET mypath=%~dp0
echo %mypath:~0,-1%
echo "%mypath:~0,-1%\transposon_annotation_tools_mitefinderii\miteFinder_windows_x64.exe -pattern_scoring %mypath:~0,-1%\transposon_annotation_tools_mitefinderii\pattern_scoring.txt %*"
%mypath:~0,-1%\transposon_annotation_tools_mitefinderii\miteFinder_windows_x64.exe -pattern_scoring %mypath:~0,-1%\transposon_annotation_tools_mitefinderii\pattern_scoring.txt %*
