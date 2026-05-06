# Sibernetic CUDA backend (skeleton)

> **Status: scaffolding only.** The single `sphFluid.cu` placeholder in
> this directory plus this README sketch the structure for a native
> CUDA differentiable substrate that mirrors `src/metal_diff/` rather
> than the legacy OpenCL `owSolver` abstract base. The actual port is
> deferred — it is ~2 weeks of focused CUDA work.

## Why this exists

The repo already has two GPU paths — OpenCL (gold-standard PCISPH on
Linux/older Macs) and native Metal (`src/metal_diff/sib_metal`, our
differentiable XPBD substrate on Apple Silicon). NVIDIA hardware is
covered by OpenCL today, but Apple killed OpenCL on Apple Silicon, so
long-term we want each platform's *vendor-backed* API:

| Platform | Backend | Status |
|---|---|---|
| Apple Silicon | native Metal (`src/metal_diff/sib_metal`) | ✅ working differentiable XPBD substrate |
| NVIDIA | **native CUDA (this directory)** | 🔧 **scaffolded; not implemented** |
| Linux server | OpenCL via NVIDIA runtime | parity baseline; do not invest |

## Target architecture: mirror `src/metal_diff/`

The current strategic direction (2026-05) is that the CUDA substrate
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

### Suggested layout

```
src/cuda_diff/
├── build.sh                    # nvcc compile + link
├── cuda_common.{h,cu}          # CudaCtx, allocate_pool, build_static_grid
│                                 (port of metal_common.h / .mm)
├── ops_kernels_m6.cu           # M6 kernels: dist_*, wpoly6, rowsum, density_grad
├── ops_xpbd_step.cu            # M7 imperative pipeline
├── ops_xpbd_full.cu            # differentiable forward + backward
├── ops_pair_spring.cu          # pair forces + spring bonds
├── shaders.cu                  # all __global__ kernel definitions
│                                 (port of shaders.metal — single file)
├── sib_cuda.cu                 # main + op dispatcher
├── load_config.py              # shared with metal_diff (no changes)
├── dump_cuda_trajectory.py     # port of dump_metal_trajectory.py
├── sgd_true.py                 # shared with metal_diff (just change BIN path)
└── test_*.py                   # FD-validated per-kernel tests, mirroring
                                  metal_diff/test_*.py
```

The math doesn't change. XPBD's constraint formulation, the kernel
signatures, the tested backward derivations all carry over directly.

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

## Phased implementation plan

Sequenced so each phase produces a working artifact:

### Phase 1: single-kernel parity (1 day)
Pick `wpoly6_inplace` (simplest M6 kernel), port to CUDA, run against
an FD reference, verify bit-equality with the Metal output. Get the
build (plain nvcc or CMake) working before any further kernels.

### Phase 2: M6 atomic ops (~3 days)
Port the rest of M6 (`dist_active_static`, `dist_active_active`,
`rowsum_density`, `density_constraint_grad`). Each gets an FD test
mirroring `test_dist.py` / `test_dens_grad.py` / `test_grad.py`.

### Phase 3: imperative `xpbd_step` (~2 days)
Port the M7 imperative pipeline (predict → density solve → distance
constraint → pair forces → floor → update_velocity). Cube-drop smoke
test (`test_xpbd.py` analog) should produce the same output as the
Metal version within float32 noise.

### Phase 4: differentiable pipeline (~3 days)
Port `xpbd_full_fwd` / `xpbd_full_bwd`, add the per-step gradient-clip
env var (`BWD_CLIP_NORM`), rerun `sgd_true.py` on demo1 (should
converge to the same optima as on Metal).

### Phase 5: parity sweep + cross-backend regression (~1 day)
Add the CUDA backend to `scripts/cross_backend_regression.py` and run
all three substrates (OpenCL, Metal, CUDA) against demo1. Gradients
should agree between Metal and CUDA to within 1 % (the OpenCL path
stays forward-only).

### Total: ~2 weeks of focused work for a competent CUDA developer.

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
`scripts/cross_backend_regression.py`: when the native CUDA backend
lands, its outputs must match OpenCL within tolerance, just as the
Metal substrate's outputs do today.

## Stub file in this directory

`sphFluid.cu` here is a placeholder from an earlier plan that mirrored
the OpenCL kernels in a single .cu file. It is *not* the right starting
point under the current strategy — start fresh from `src/metal_diff/`
following the layout above. The file is kept for reference only.
