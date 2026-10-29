@echo off
cd /d %~dp0
set MADX=C:\GIT_REPOS\NICA\MADX\2.Programms\Twiss\madx.exe
%MADX% chrom_plain.madx > output\madx_plain.log 2>&1
%MADX% chrom_ptc.madx > output\madx_ptc.log 2>&1
type output\madx_plain_results.txt
