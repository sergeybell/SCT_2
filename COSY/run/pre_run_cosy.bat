@echo off
rem Compile base COSY modules (cosy, utilities, elements, header) in COSY\src.
setlocal
pushd "%~dp0..\src"
if not exist tmp mkdir tmp
for %%F in (cosy utilities elements header) do (
    echo ========================================
    echo RUNNING: %%F.fox
    echo ========================================
    cosy.exe %%F.fox || goto :fail
    echo.
)
popd
exit /b 0

:fail
popd
exit /b 1
