@echo off
cd /d %~dp0
julia --project=C:\GIT_REPOS\BELL\Julia chrom_bell.jl > output\bell_stdout.log 2>&1
type output\bell_stdout.log
