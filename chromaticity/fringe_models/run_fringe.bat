@echo off
cd /d %~dp0
if not exist output mkdir output
julia --project=. run.jl > output\fringe_stdout.log 2>&1
type output\fringe_stdout.log
