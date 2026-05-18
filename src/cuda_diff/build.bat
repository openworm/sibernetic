@echo off
REM Windows CUDA build for the differentiable Sibernetic substrate.
REM Mirrors src/metal_diff/build.sh but targets nvcc + MSVC.
REM
REM Requirements:
REM   - CUDA Toolkit 12.0+ (tested on 12.4 in CI, 13.2 locally)
REM   - MSVC via Visual Studio Build Tools or full VS install
REM     (cl.exe must be reachable after running vcvarsall.bat)
REM
REM VS lookup order (override by setting VCVARS before invocation):
REM   1. %VCVARS% if set
REM   2. vswhere.exe (works for VS 2017/2019/2022/Insiders/preview)
REM   3. Hardcoded fallback paths for VS 2022/2019 Community + BuildTools
REM
REM Architecture: sm_75 (Turing/RTX 20-series) as the floor. Bump via
REM   set CUDA_ARCH=sm_86  (Ampere / RTX 30-series)
REM   set CUDA_ARCH=sm_89  (Ada / RTX 40-series / L4)
REM
REM Layout: six translation units linked via nvcc -rdc=true. Device-link
REM is performed by nvcc automatically in the final invocation.
REM     cuda_common.cu         host I/O helpers
REM     shaders.cu             all __global__ kernel definitions
REM     ops_kernels_m6.cu      M6 fwd/bwd + M7.1 bwd standalone drivers
REM     ops_xpbd_step.cu       run_xpbd_step (all-in-one orchestrator)
REM     ops_xpbd_full.cu       xpbd_full / xpbd_step_diff / xpbd_step_full_diff
REM     sib_cuda.cu            main() + CLI dispatcher

setlocal
if "%CUDA_ARCH%"=="" set CUDA_ARCH=sm_75

REM -- Locate vcvarsall.bat ----------------------------------------------
REM Priority 1: caller-provided override.
if not "%VCVARS%"=="" goto :have_vcvars

REM Priority 2: vswhere.exe (Microsoft's recommended way; ships with
REM every VS 2017+ install and on most build agents). Asks for the
REM latest install that has the VC++ tools component, writing the
REM resolved install path to a temp file. Implementation notes:
REM   - %ProgramFiles(x86)% contains a literal `)` that wrecks
REM     parenthesized cmd blocks; we use `goto` instead of an `if (...)`
REM     body to dodge the parser issue.
REM   - We redirect vswhere's stdout to a temp file rather than capture
REM     it via `for /f` backticks; the latter strips inner quotes when
REM     the executable path contains spaces, causing a misleading
REM     "'vswhere.exe' is not recognized" error.
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" set "VSWHERE=%ProgramFiles%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" goto :try_fallback_paths

set "VSWHERE_OUT=%TEMP%\sib_cuda_vswhere.txt"
"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath >"%VSWHERE_OUT%" 2>nul
for /f "usebackq delims=" %%I in ("%VSWHERE_OUT%") do (
    if exist "%%I\VC\Auxiliary\Build\vcvarsall.bat" set "VCVARS=%%I\VC\Auxiliary\Build\vcvarsall.bat"
)
del "%VSWHERE_OUT%" >nul 2>nul
if not "%VCVARS%"=="" goto :have_vcvars

:try_fallback_paths

REM Priority 3: hardcoded fallbacks. Cover both x64 Program Files (legacy
REM BuildTools layout) and x86 Program Files (Community / older).
for %%E in (Community BuildTools Professional Enterprise) do (
    for %%R in ("%ProgramFiles%" "%ProgramFiles(x86)%") do (
        for %%V in (2022 2019 2017) do (
            if exist "%%~R\Microsoft Visual Studio\%%V\%%E\VC\Auxiliary\Build\vcvarsall.bat" (
                set "VCVARS=%%~R\Microsoft Visual Studio\%%V\%%E\VC\Auxiliary\Build\vcvarsall.bat"
                goto :have_vcvars
            )
        )
    )
)

echo ERROR: Could not find vcvarsall.bat. Tried:
echo   - %%VCVARS%% override (unset)
echo   - vswhere.exe at %ProgramFiles(x86)%%\Microsoft Visual Studio\Installer\vswhere.exe
echo   - VS 2022/2019/2017 Community/BuildTools/Professional/Enterprise
echo     under both Program Files and Program Files (x86).
echo.
echo Install Visual Studio Build Tools (Desktop development with C++)
echo or set VCVARS=^<path-to-vcvarsall.bat^> before running this script.
exit /b 1

:have_vcvars
echo Using vcvarsall: %VCVARS%
call "%VCVARS%" x64 >nul
if errorlevel 1 (
    echo ERROR: vcvarsall.bat failed for x64.
    exit /b 1
)

set HERE=%~dp0
cd /d "%HERE%"

nvcc -std=c++17 -O2 -arch=%CUDA_ARCH% -rdc=true ^
    cuda_common.cu ^
    shaders.cu ^
    ops_kernels_m6.cu ^
    ops_xpbd_step.cu ^
    ops_xpbd_full.cu ^
    ops_pair_spring.cu ^
    sib_cuda.cu ^
    -o sib_cuda.exe
if errorlevel 1 exit /b 1

echo Built %HERE%sib_cuda.exe (arch=%CUDA_ARCH%)
endlocal
