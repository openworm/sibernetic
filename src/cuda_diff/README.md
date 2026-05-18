# Sibernetic CUDA differentiable substrate

CUDA port of `src/metal_diff/`, implementing the work plan in
`src/cuda/README.md`. Mirrors the Metal substrate file-for-file so the
shared Python harnesses (`load_config.py`, `sgd_true.py`,
`dump_*_trajectory.py`) work against either backend by swapping only the
binary path.

## Layout

```
src/cuda_diff/
  build.sh / build.bat        nvcc -arch=sm_75 -rdc=true
  cuda_common.{h,cu}          CudaCtx, allocate_pool, build_static_grid
  shaders.{h,cu}              all __global__ kernels (port of shaders.metal)
  ops.h                       op-dispatcher declarations
  ops_kernels_m6.cu           M6 host drivers (dist_*, wpoly6, density, ...)
  ops_xpbd_step.cu            M7 imperative XPBD pipeline driver
  ops_xpbd_full.cu            differentiable forward + backward drivers
  ops_pair_spring.cu          pair-force + spring-bond drivers
  sib_cuda.cu                 main() + op dispatcher
  load_config.py              shared with metal_diff (no changes)
  sgd_true.py                 shared with metal_diff (BIN path only)
  dump_cuda_trajectory.py     xpbd_step chunked dump
  dump_cuda_full_trajectory.py xpbd_full_fwd one-shot dump
  test_*.py                   per-kernel and end-to-end tests
  run_all_tests.{sh,bat}      regression-runner entry points
```

## Build

Windows (primary):  `cmd /c src\cuda_diff\build.bat`
Linux (unverified): `./src/cuda_diff/build.sh`

Requires CUDA Toolkit 12.0+ (CI exercises 12.4; locally tested up to
13.2), VS Build Tools 2017/2019/2022 on Windows (`build.bat` probes for
them via `vswhere.exe`), and an sm_75 (Turing) or newer GPU. Outputs
`sib_cuda` / `sib_cuda.exe` next to the sources.

## Run tests

    src/cuda_diff/run_all_tests.sh         # Linux/macOS
    cmd /c src\cuda_diff\run_all_tests.bat # Windows

Both iterate over `test_*.py` sequentially and exit 0 iff every test
PASSes. No pytest dependency; each test is a standalone script that
builds its fixture, calls `sib_cuda` via subprocess, and validates
outputs in-process.

## Spec coverage (per src/cuda/README.md)

* Phase 1 (single-kernel parity)        complete
* Phase 2 (M6 atomic ops)               complete
* Phase 3 (imperative xpbd_step)        complete
* Phase 4 (differentiable XPBD)         complete, FD-validated K=2..1000
* Phase 5 (cross-backend regression)    fwd parity scripted: CUDA and
  metal-native substrates are wired into scripts/cross_backend_regression.py
  via a paired --substrate flag, so demo1 can be driven OpenCL / Metal /
  CUDA from one entry point. Metal-vs-CUDA gradient parity within 1 %
  is still blocked on Apple-silicon hardware in this dev environment.

## Known limitations

* Membranes (M10) are out-of-scope per the spec scope rule. Orchestrators
  accept the membrane CLI args (matching Metal's positional contract)
  but no-op them and emit a stderr warning when `n_membranes > 0`.
* Only sm_75 has been exercised on real hardware (RTX 2070, both Windows
  and Linux). sm_70 / sm_86 / sm_89 should compile but are untested.
* Linux build + runtime verified on Ubuntu 24.04 + CUDA 12.6 via WSL2
  against the same RTX 2070: `build.sh` produces a working ELF binary
  and all 13 tests in `run_all_tests.sh` PASS with the same FD
  tolerances as the Windows runs. Native (non-WSL) Linux is untested
  but the WSL path exercises the same nvcc + glibc + driver stack.
* `test_demo1_edge_cases.py` flags spring_K=1e5 as unstable at demo1's
  dt/sim_scale combination -- physics divergence, not a kernel bug. See
  the V5 commit message for details.

## Upstream reference

The line-by-line counterpart to every file here is in `src/metal_diff/`.
When in doubt about a kernel's contract or a host driver's argv layout,
look at the matching `.mm` / `.metal` first - the CUDA file is intended
to be a mechanical translation.

Validation history lives in the git log for `ow-native-gpu-0.1.0`
(search commit messages for `V1` .. `V5`).
