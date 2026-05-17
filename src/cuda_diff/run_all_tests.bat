@echo off
REM Regression runner for the CUDA differentiable substrate (Windows).
REM
REM Runs every test_*.py in this directory sequentially. Each test is
REM self-contained: builds its own fixture, invokes sib_cuda.exe directly,
REM validates outputs in-process. Exit 0 iff all PASS.
REM
REM Mirrors run_all_tests.sh -- same loop, batch syntax.
REM
REM Prerequisite: sib_cuda.exe built. Run "cmd /c build.bat" first.

setlocal enabledelayedexpansion

set "HERE=%~dp0"
pushd "%HERE%..\.."

set "PY=python"
where %PY% >nul 2>&1
if errorlevel 1 set "PY=python3"

set /a n_pass=0
set /a n_fail=0
set "failed="

for %%t in ("%HERE%test_*.py") do (
    echo === %%~nt ===
    %PY% "%%t"
    if errorlevel 1 (
        echo [FAIL] %%~nt ^(exit !errorlevel!^)
        echo.
        set /a n_fail+=1
        set "failed=!failed! %%~nt"
    ) else (
        echo [PASS] %%~nt
        echo.
        set /a n_pass+=1
    )
)

echo === Summary ===
echo   passed: !n_pass!
echo   failed: !n_fail!
if !n_fail! gtr 0 (
    echo   failed tests:!failed!
    popd
    exit /b 1
)
echo [OVERALL PASS] all !n_pass! CUDA tests
popd
exit /b 0
