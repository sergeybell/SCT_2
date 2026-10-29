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
rem fringe fields: FR 1/2 need SYSCA.DAT in the working dir (made once by sysca\mk_sysca.fox = COSY's CRSYSCA, ~10 min)
if not exist sysca\SYSCA.DAT (cd sysca & %COSY% mk_sysca.fox > mk_sysca.log 2>&1 & cd ..)
copy /y sysca\SYSCA.DAT . > nul
%COSY% lattice_md_fr.fox > output\lattice_md_fr.log 2>&1
%COSY% ov321_fr.fox > output\ov321_fr.log 2>&1
%COSY% ov330_fr.fox > output\ov330_fr.log 2>&1
%COSY% rf_ts.fox > output\rf_ts.log 2>&1
python ov330_map_analysis.py
python fr_analysis.py
python rf_ts_analysis.py
del /q SYSCA.DAT RKLOG.DAT 2> nul

