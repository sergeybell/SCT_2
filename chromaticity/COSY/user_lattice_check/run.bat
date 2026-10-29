@echo off
cd /d %~dp0
set COSY=C:\GIT_REPOS\SCT_2\COSY\src\cosy.exe
%COSY% ov321_user_lattice.fox > user_lattice.out 2>&1
%COSY% magnetic_2_notilt.fox > notilt_compile.out 2>&1
%COSY% ov321_user_lattice_notilt.fox > user_lattice_notilt.out 2>&1
findstr RES user_lattice.out user_lattice_notilt.out
