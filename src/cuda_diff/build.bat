@echo off
REM Windows CUDA build for the differentiable Sibernetic substrate.
REM Mirrors src/metal_diff/build.sh but targets nvcc + MSVC.
REM
REM Requirements:
REM   - CUDA Toolkit 12.0+ (tested on 13.2)
REM   - MSVC via Visual Studio Build Tools (cl.exe must be reachable
REM     after running vcvarsall.bat)
REM
REM Architecture: sm_75 (Turing/RTX 20-series) as the floor. Bump via
REM   set CUDA_ARCH=sm_86  (Ampere / RTX 30-series)
REM   set CUDA_ARCH=sm_89  (Ada / RTX 40-series / L4)
REM
REM Layout (post #15 refactor): six translation units linked via
REM   nvcc -rdc=true. Device-link is performed by nvcc automatically in
REM   the final invocation.
REM     cuda_common.cu         host I/O helpers
REM     shaders.cu             all __global__ kernel definitions
REM     ops_kernels_m6.cu      M6 fwd/bwd + M7.1 bwd standalone drivers
REM     ops_xpbd_step.cu       run_xpbd_step (all-in-one orchestrator)
REM     ops_xpbd_full.cu       xpbd_full / xpbd_step_diff / xpbd_step_full_diff
REM     sib_cuda.cu            main() + CLI dispatcher

setlocal
if "%CUDA_ARCH%"=="" set CUDA_ARCH=sm_75
if "%VCVARS%"=="" set VCVARS=C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvarsall.bat

if not exist "%VCVARS%" (
    echo ERROR: vcvarsall.bat not found at "%VCVARS%"
    echo Set VCVARS=^<path-to-vcvarsall.bat^> before running.
    exit /b 1
)
call "%VCVARS%" x64 >nul

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
