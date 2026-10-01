@echo off
rem Usage:
rem   run_cosy.bat mapping                 job from COSY\jobs
rem   run_cosy.bat lattice magnetic_2      compile COSY\structures\magnetic_2 (lattice + maps)
rem   run_cosy.bat lattice magnetic_2 mapping
setlocal
if "%~1"=="" (
    echo Usage: run_cosy.bat [lattice STEM] JOB [JOB ...]
    exit /b 1
)
pushd "%~dp0..\src"
if not exist tmp mkdir tmp

if /I "%~1"=="lattice" (
    if not exist "..\dat\%~2" mkdir "..\dat\%~2"
    call :run "..\structures\%~2\%~2.fox" || goto :fail
    call :run "..\structures\%~2\%~2_maps.fox" || goto :fail
    shift
    shift
)

:loop
if "%~1"=="" goto :done
call :run "..\jobs\%~n1.fox" || goto :fail
shift
goto :loop

:run
echo ========================================
echo RUNNING: %~1
echo ========================================
cosy.exe %1
echo.
exit /b %errorlevel%

:done
popd
exit /b 0

:fail
popd
exit /b 1
