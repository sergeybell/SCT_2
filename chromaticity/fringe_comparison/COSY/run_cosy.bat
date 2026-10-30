@echo off
rem COSY 10.2 (cosy.exe from C:\GIT_REPOS\SCT_2\COSY\src). FR 3 needs no SYSCA.DAT.
cd /d %~dp0
if not exist output mkdir output
set COSY=C:\GIT_REPOS\SCT_2\COSY\src\cosy.exe
%COSY% lattice_cmp.fox > output\lattice_cmp.log 2>&1
%COSY% chrom_cmp.fox > output\chrom_cmp.log 2>&1
