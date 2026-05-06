# Sibernetic CUDA backend (skeleton)

> **Status: scaffolding only.** Files in this directory and `inc/owCudaSolver.h` lay out the structure for a native CUDA backend that mirrors PR #222's native Metal port. The actual CUDA kernel ports (translating `src/sphFluid.cl` to `src/cuda/sphFluid.cu`) are deferred — they are weeks of focused work and depend on PR #222's `owSolver` abstract base landing first.

## Why this exists

The `ow-native-gpu-0.1.0` line is built on the strategic decision (see `DEVELOPMENT_LOG.md`) to use **two vendor-backed GPU backends** instead of relying on cross-platform Python compilers like Taichi (whose maintenance has slowed) or PyTorch (whose per-kernel-launch overhead makes it 21× slower than OpenCL on the same hardware):

| Platform | Backend | Status |
|---|---|---|
| Apple Silicon | native Metal (`src/metal_diff/sib_metal`) | working, differentiable XPBD substrate |
| NVIDIA | **native CUDA (this directory)** | **scaffolded; not implemented** |
| Linux server | OpenCL via NVIDIA runtime (existing) | parity baseline; do not invest |

## Structure (mirrors PR #222 / Metal)

```
src/cuda/
├── README.md          ← you are here
├── sphFluid.cu        ← all CUDA __global__ kernels (port of src/sphFluid.cl)
├── CudaContext.cpp/h  ← CUDA device init, stream, memory pools (TODO)
└── kernels/           ← one .cuh per kernel descriptor (mirrors PR #222's src/kernels/) (TODO)

inc/
└── owCudaSolver.h     ← public C++ interface, mirrors owOpenCLSolver

src/
├── owCudaSolver.cpp   ← bridge from owSolver virtual interface to .cu kernels (TODO)
└── backend/
    └── CudaBackend.cpp/h  ← (TODO) CUDA runtime API wrapper, equivalent to PR #222's MetalBackend
```

## Implementation plan

Sequenced so each step produces a working artifact:

### Phase 0: Wait for PR #222 to land
Reason: PR #222 introduces `inc/owSolver.h` (the abstract base both backends implement) and the `src/kernels/` descriptor pattern. Building the CUDA backend on the pre-PR-#222 OpenCL-only structure means refactoring once #222 lands. Wait until #222 merges.

### Phase 1: Port `sphFluid.cl` → `sphFluid.cu` (literal translation)
Translate every OpenCL kernel in `src/sphFluid.cl` (1515 lines) to CUDA. Mostly mechanical:
- `__kernel void` → `__global__ void`
- `__global float4 *buf` → `float4 *buf` with explicit pointer args
- `get_global_id(0)` → `blockIdx.x * blockDim.x + threadIdx.x`
- `barrier(CLK_LOCAL_MEM_FENCE)` → `__syncthreads()`
- `__local` → `__shared__`
- OpenCL math intrinsics → CUDA intrinsics (`fabs` → `fabsf`, etc.)

Estimated: 2-3 days for a careful translation, with a parity test against the OpenCL output at each kernel.

### Phase 2: Implement `CudaBackend.cpp` (host-side dispatch)
Wraps cuBLAS-style CUDA runtime calls:
- `cudaMalloc` / `cudaFree` for buffers
- `cudaMemcpy` for host-device transfers
- `<<<grid, block>>>` kernel launches
- `cudaStreamSynchronize` for ordering

Estimated: 2-3 days.

### Phase 3: Implement `owCudaSolver.cpp`
Bridge between `owSolver` virtual interface (PR #222's abstraction) and `CudaBackend.cpp`'s kernel launches. Mirrors `owMetalSolver.cpp` from PR #222 line-by-line.

Estimated: 1-2 days.

### Phase 4: Wire into the build
- Add `nvcc` to the Linux makefile path
- Add CUDA backend selection to `owConfigProperty.cpp` (`backend=cuda`)
- Update `Dockerfile` for sibernetic-runner to install CUDA toolkit (already has CUDA runtime via nvidia/cuda image)
- Add `backend=cuda` to the cross-backend regression script

Estimated: 1 day.

### Phase 5: Cross-backend parity validation
Run `scripts/cross_backend_regression.py --backend cuda --backend opencl --local-binary <PR222-Metal>`. All three should produce demo1 cube-stability metrics within the existing tolerance bands (extent retention ≥ 80%, mean_y fell ≥ 50%).

Estimated: 1 day of measurement + tuning.

### Total estimated effort: ~2 weeks of focused work for a competent CUDA developer.

## Reference files in PR #222 to model from

When PR #222 lands, the matching CUDA files would mirror these structurally:

| Metal file | CUDA equivalent |
|---|---|
| `inc/owMetalSolver.h` | `inc/owCudaSolver.h` |
| `src/owMetalSolver.cpp` | `src/owCudaSolver.cpp` |
| `src/owMetalPrivateImpl.cpp` | `src/owCudaPrivateImpl.cpp` (if needed) |
| `src/backend/MetalBackend.{cpp,h}` | `src/backend/CudaBackend.{cpp,h}` |
| `src/metal/sphFluid.metal` | `src/cuda/sphFluid.cu` |
| `src/kernels/*.h` | `src/kernels/*.h` (already shared with Metal — same descriptors) |

The Metal/CUDA divergence is **only** in the kernel language (MSL vs CUDA C++) and the host-side runtime API (Metal C++ vs CUDA Runtime). The algorithm specification, kernel descriptors, and abstract solver interface are shared.

## Why not just use OpenCL on NVIDIA?

It actually works fine — Cloud Run + L4 + NVIDIA's OpenCL runtime measures at 86 sec for a 1-sec demo1 sim, with cube physics intact. However:
- Apple killed OpenCL on Apple Silicon, so it's not a path forward for cross-platform dev
- The 2015 AMD APP SDK we historically link against is abandoned
- NVIDIA's OpenCL is still maintained but not actively invested in
- For long-term maintainability we want vendor-backed APIs (CUDA on NVIDIA, Metal on Apple)

OpenCL on NVIDIA stays as the **parity baseline** in the cross-backend regression: when we add the native CUDA backend, its outputs must match OpenCL within tolerance.
