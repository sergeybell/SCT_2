@echo off
rem COSY 10.2: cosy.exe from C:\GIT_REPOS\SCT_2\COSY\src (its cosy.bin is INCLUDEd by absolute path)
cd /d %~dp0
set COSY=C:\GIT_REPOS\SCT_2\COSY\src\cosy.exe
%COSY% lattice_md.fox > output\lattice_md.log 2>&1
%COSY% ov321.fox > output\ov321.log 2>&1
%COSY% ov330.fox > output\ov330.log 2>&1
%COSY% ov321_user_method.fox > output\ov321_user_method.log 2>&1
%COSY% twiss.fox > output\twiss.log 2>&1
if not "%1"=="" %COSY% rf_test.fox > output\rf_test.log 2>&1
python ov330_map_analysis.py

