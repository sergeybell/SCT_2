@echo off
rem OptiM is GUI-only: open magnetic_2_md.opt and press Ctrl+T (Twiss), then copy the summary to output\optim_reference.txt
start "" "C:\Program Files\OptiMX-x86-64\optimx64.exe" "%~dp0magnetic_2_md.opt"
