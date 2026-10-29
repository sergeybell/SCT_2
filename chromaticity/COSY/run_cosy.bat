@echo off
rem COSY must run in a folder containing compiled cosy.bin (C:\GIT_REPOS\SCT_2\COSY\src)
pushd C:\GIT_REPOS\SCT_2\COSY\src
cosy.exe %~dp0chrom_ov321.fox > %~dp0output\cosy_ov321.log 2>&1
cosy.exe %~dp0chrom_ov330.fox > %~dp0output\cosy_ov330.log 2>&1
popd
type %~dp0output\cosy_ov321.txt %~dp0output\cosy_ov330.txt
