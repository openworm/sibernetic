# Sibernetic CUDA backend — original scaffolding plan

> **Status: historical. The work described below was completed and
> lives in [`src/cuda_diff/`](../cuda_diff/). The single `sphFluid.cu`
> in this directory is a leftover from an earlier OpenCL-mirroring plan
> and is no longer the starting point.**
>
> For the as-built native CUDA differentiable substrate, build status,
> per-phase test suite, and known limitations, see
> [`src/cuda_diff/README.md`](../cuda_diff/README.md). The rest of this
> file is preserved verbatim as the original work plan that motivated
> the `src/cuda_diff/` layout.

## Why this exists

The repo already has two GPU paths — OpenCL (gold-standard PCISPH on
Linux/older Macs) and native Metal (`src/metal_diff/sib_metal`, our
differentiable XPBD substrate on Apple Silicon). NVIDIA hardware is
covered by OpenCL today, but Apple killed OpenCL on Apple Silicon, so
long-term we want each platform's *vendor-backed* API:

| Platform | Backend | Status |
|---|---|---|
| Apple Silicon | native Metal (`src/metal_diff/sib_metal`) | working differentiable XPBD substrate |
| NVIDIA | native CUDA (`src/cuda_diff/sib_cuda`) | complete; FD-validated through K=1000; see `src/cuda_diff/README.md` |
| Linux server | OpenCL via NVIDIA runtime | parity baseline; do not invest |

## Target architecture: mirror `src/metal_diff/`

The strategic direction (set 2026-05) was that the CUDA substrate
should mirror the metal_diff differentiable substrate file-for-file,
NOT the abstract `owSolver` virtual interface from earlier plans.
Reasons:

- The metal_diff substrate has worked out the differentiable contract:
  every forward kernel pairs with a backward kernel, shared types live
  in a single header (`metal_common.h`), the run-op dispatcher is a
  thin extern table, and the build is a single shell script.
- A CUDA port that follows that pattern gets analytic gradients on
  NVIDIA hardware essentially for free.
- OpenCL's `owOpenCLSolver` predates the differentiable design and is
  forward-only. Mirroring it would re-create that limitation.

### As-built layout

The `src/cuda_diff/` substrate followed the suggested layout almost
exactly:

```
src/cuda_diff/
├── build.sh / build.bat        # nvcc compile + link (sm_75 floor)
├── cuda_common.{h,cu}          # CudaCtx, allocate_pool, build_static_grid
│                                 (port of metal_common.h / .mm)
├── shaders.cu                  # all __global__ kernel definitions
│                                 (port of shaders.metal — single file)
├── shaders.h                   # kernel prototypes
├── ops.h                       # op-dispatcher declarations
├── ops_kernels_m6.cu           # M6 fwd/bwd standalone drivers
├── ops_xpbd_step.cu            # M7 imperative xpbd_step orchestrator
├── ops_xpbd_full.cu            # differentiable xpbd_full_fwd / _bwd
├── ops_pair_spring.cu          # pair_forces + spring_bonds (fwd/bwd)
├── sib_cuda.cu                 # main + op dispatcher
├── load_config.py              # shared with metal_diff
├── dump_cuda_trajectory.py     # xpbd_step chunked dump
├── dump_cuda_full_trajectory.py # xpbd_full_fwd one-shot dump
├── sgd_true.py                 # shared with metal_diff (BIN path differs)
└── test_*.py                   # 13 FD-validated tests
```

The math doesn't change. XPBD's constraint formulation, the kernel
signatures, the tested backward derivations all carried over directly.

## Per-kernel translation rules (mostly mechanical)

| MSL (Metal Shading Language) | CUDA C++ |
|---|---|
| `kernel void foo(...)` | `__global__ void foo(...)` |
| `device float *buf [[buffer(0)]]` | `float *buf` (positional arg) |
| `[[thread_position_in_grid]]` | `blockIdx.x * blockDim.x + threadIdx.x` |
| `threadgroup float partials[256]` | `__shared__ float partials[256]` |
| `threadgroup_barrier(mem_flags::mem_threadgroup)` | `__syncthreads()` |
| `dispatchThreads:MTLSizeMake(N,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)` | `kernel<<<(N+255)/256, 256>>>(...)` |
| `MTLBuffer` | `cudaMalloc` + raw pointer |
| `[buffer contents]` (host pointer) | `cudaMemcpyDeviceToHost` into staging buffer |

OpenCL math intrinsics map similarly: `fabs` → `fabsf`, `sqrt` →
`sqrtf`, `pow` → `powf`, `dot()` → manual `(a.x*b.x + a.y*b.y + a.z*b.z)`
since CUDA `float3` doesn't have built-in `dot`.

## Phased implementation plan (as executed)

Sequenced so each phase produced a working artifact:

### Phase 1: single-kernel parity (1 day) — complete
Picked `wpoly6_inplace` (simplest M6 kernel), ported to CUDA, ran
against an FD reference, verified bit-equality with the Metal output.
Got the nvcc build working before any further kernels.

### Phase 2: M6 atomic ops (~3 days) — complete
Ported the rest of M6 (`dist_active_static`, `dist_active_active`,
`rowsum_density`, `density_constraint_grad`). Each got an FD test
mirroring `test_dist.py` / `test_dens_grad.py` / `test_grad.py`.

### Phase 3: imperative `xpbd_step` (~2 days) — complete
Ported the M7 imperative pipeline (predict → density solve → distance
constraint → pair forces → floor → update_velocity). Cube-drop smoke
test (`test_cube_drop_cuda.py`) produces the same output as the Metal
version within float32 noise.

### Phase 4: differentiable pipeline (~3 days) — complete
Ported `xpbd_full_fwd` / `xpbd_full_bwd`, added the per-step
gradient-clip env var (`BWD_CLIP_NORM`), reran `sgd_true.py` on demo1
(converges to the same optima as on Metal).

### Phase 5: parity sweep + cross-backend regression (~1 day) — partial
Added the CUDA backend to `scripts/cross_backend_regression.py`.
Forward parity scripted: OpenCL / Metal / CUDA all run demo1 from one
entry point. Metal-vs-CUDA gradient parity within 1% is blocked on
Apple-silicon hardware in this dev environment (see Phase-5 note in
`src/cuda_diff/README.md`).

## Why not just use OpenCL on NVIDIA?

It actually works fine — a remote L4 GPU running NVIDIA's OpenCL
runtime measures at 86 sec for a 1-sec demo1 sim with cube physics
intact. However:

- Apple killed OpenCL on Apple Silicon, so it is not a *single-vendor*
  cross-platform path forward.
- The 2015 AMD APP SDK Sibernetic historically links against is
  abandoned upstream.
- NVIDIA's OpenCL is maintained but receives little ongoing investment;
  CUDA gets all the new features (graph capture, cooperative groups,
  FP16/BF16 compute, sparse warp-level primitives).
- For *differentiable* physics on NVIDIA hardware specifically, CUDA
  graph capture and cooperative-group reductions speed up the per-step
  backward by 2–3× over OpenCL on the same hardware.

OpenCL on NVIDIA stays as the **parity baseline** in
`scripts/cross_backend_regression.py`: native CUDA outputs must match
OpenCL within tolerance, just as the Metal substrate's outputs do
today.

## Stub file in this directory

`sphFluid.cu` here is a placeholder from an earlier plan that mirrored
the OpenCL kernels in a single .cu file. It is *not* the right starting
point under the as-built strategy in `src/cuda_diff/`. The file is
kept for reference only and is not part of any build target.
