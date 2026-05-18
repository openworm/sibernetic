# Sibernetic CUDA differentiable substrate

CUDA port of `src/metal_diff/`, implementing the work plan in
`src/cuda/README.md`. Mirrors the Metal substrate's kernel set and the
integrated XPBD CLI surface (`xpbd_full_fwd` / `xpbd_full_bwd` /
`xpbd_step`) so the shared Python harnesses (`load_config.py`,
`sgd_true.py`, `dump_*_trajectory.py`) work against either backend by
swapping only the binary path. Unit-level CLI ops (the standalone
`pair_forces_*`, `spring_bonds_*`, density-chain primitives) follow
similar argv conventions but were renamed where the underlying kernel
gained a "_grid" qualifier; see the **Op-name map** below. Membranes
(M10) are intentionally out-of-scope — see **Known limitations**.

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
* Phase 5 (cross-backend regression)    fwd parity scripted: OpenCL,
  metal-native, and CUDA substrates are wired into
  scripts/cross_backend_regression.py via a paired --substrate flag, so
  demo1 can be driven through any of them from one entry point. Each
  substrate's analytic backward gradients are independently FD-validated
  in-tree (CUDA: this directory's `test_xpbd_full_bwd_*` suite, all PASS
  on RTX 2070; Metal: src/metal_diff/test_*.py on M-series). The
  cross-substrate Metal-vs-CUDA gradient comparison (~1 % tolerance) is
  deferred because it requires simultaneous Apple-silicon + NVIDIA
  access which isn't available in this dev environment; it is not a
  defect in either substrate.

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

## Op-name map (CUDA ↔ Metal)

Integrated XPBD ops (`xpbd_full_fwd`, `xpbd_full_bwd`, `xpbd_step`) have
identical names and argv layouts in both substrates, so the shared
Python harnesses are substrate-agnostic. Standalone kernel ops differ:

| CUDA op (`sib_cuda <op>`) | Metal op (`sib_metal <op>`) | Notes |
|---|---|---|
| `pair_forces_grid_fwd` / `_bwd` | `pair_forces_fwd` / `_bwd` | CUDA name reflects the spatial-grid acceleration in the kernel body |
| `spring_bonds_force` / `_bwd` | `spring_bonds_fwd` / `_bwd` | "force" = the force-based (not constraint-based) bond model used by both |
| `visc_K_partial`, `spring_K_partial` | (host-only on Metal) | CUDA exposes these as standalone ops for FD validation |
| `apply_ext_accel`, `apply_ext_accel_bwd` | (inlined in Metal) | CUDA exposes these as standalone ops for the integrated chain |
| `dist_active_static[_bwd]`, `dist_active_active[_bwd]`, `wpoly6_inplace[_bwd]`, `rowsum_density[_bwd]`, `density_constraint_grad[_bwd]`, `solve_density_constraint_bwd` | identical names | M6 / M9 chain primitives |

## Test-suite mapping

CUDA ships 13 tests vs Metal's 19 because several Metal unit tests were
consolidated into broader end-to-end checks on the CUDA side:

| CUDA test | Metal-side coverage |
|---|---|
| `test_m6_cuda.py`, `test_m6_bwd_cuda.py` | `test_dist.py` + `test_dens_grad.py` + `test_grad.py` + `test_solve_dens_bwd.py` |
| `test_xpbd_full_fwd_cuda.py` | `test_xpbd.py` + `test_xpbd_full.py` + `test_xpbd_full_spring.py` + `test_xpbd_full_visc.py` |
| `test_xpbd_full_bwd_cuda.py`, `test_xpbd_full_bwd_longchain_cuda.py` | `test_constraint_grad_bwd.py` + `test_param_grads.py` + `test_dens_alpha_grad.py` |
| `test_pair_forces_grid_cuda.py` | `test_pair_forces_bwd.py` |
| `test_visc_K_partial_cuda.py`, `test_m7_1_bwd_cuda.py` | per-parameter FD checks |
| `test_cube_drop_cuda.py` | demo1 integration check |
| `test_conservation_laws.py`, `test_fuzz_xpbd_full.py`, `test_demo1_edge_cases.py` | CUDA-only audits (no Metal counterpart) |

Membrane tests (`test_membrane_correction.py`, `test_xpbd_full_membrane.py`)
have no CUDA counterpart — M10 is out-of-scope.

## Upstream reference

The line-by-line kernel-level counterpart to every file here is in
`src/metal_diff/`. When in doubt about a kernel's math or the integrated
XPBD argv layout, look at the matching `.mm` / `.metal` first — the
CUDA kernel implementations are intended as mechanical ports. Unit-level
CLI argv layouts may differ at the surface; see the Op-name map above.

Validation history lives in the git log for `ow-native-gpu-0.1.0`
(search commit messages for `V1` .. `V5`).
