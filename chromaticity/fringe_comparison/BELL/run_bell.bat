@echo off
cd /d %~dp0
if not exist output mkdir output
julia --project=. run.jl > output\bell_stdout.log 2>&1
type output\bell_stdout.log
