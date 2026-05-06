# metal_diff — hand-written differentiable Metal for Sibernetic

This directory holds hand-written Metal compute kernels for the native
differentiable Sibernetic substrate. Apple Silicon GPU only for now;
native CUDA equivalent is a separate follow-up.

## Why hand-written

Earlier substrate-prototyping attempts (PyTorch on MPS, Taichi-Metal)
either lacked the autodiff coverage we needed or had algorithmic-scale
issues (the coordinate-space bug in the Taichi imperative solver, MPS
backward limitations) that didn't map cleanly to Sibernetic's particle
counts. The decision (2026-05-03) was to write the GPU kernels directly
in Metal Shading Language and pair each forward kernel with a
hand-derived backward kernel. Both legacy backends have since been
removed from the tree.

## What's here

| File | Purpose |
|---|---|
| `sib_metal.mm` | Single Objective-C++ host program. Embeds Metal shaders inline, parses `op` argument, runs the corresponding kernel. |
| `build.sh` | `clang++ -framework Metal -framework Foundation` → `sib_metal` |
| `test_dist.py` | Numpy-reference validator + steady-state benchmark for each kernel. |

## Currently supported ops

### Steady-state perf summary (demo1 scale, 343 active × 17,498 static)

| Op | ms/iter | Cumulative |
|---|---|---|
| `dist_active_static` (M6.0) | 0.358 | 0.358 |
| `wpoly6_inplace` (M6.1) | 0.508 | 0.866 |
| `rowsum_density` (M6.2) | 0.370 | 1.236 |
| `dist_active_active` (M6.3) | 0.179 | 1.415 |
| `density_constraint_grad` (M6.4) | 0.426 | **1.841** |

For reference: OpenCL's full demo1 step on L4 GPU is **1.61 ms** (sort + PCISPH + neighbor lookup + everything). M6 forward kernels at 1.84 ms cover full SPH density + density-constraint gradient — the entire forward pipeline that XPBD's density constraint projection needs.

For one full XPBD step with N inner projection iters (rough estimate):
- Setup (one-time per step): M6.0 + M6.1 + M6.2 + M6.3 = 1.42 ms
- Per inner iter: M6.4 = 0.43 ms
- For N=3 iters: ~2.7 ms/step

Slower than OpenCL's 1.61 ms/step on L4 but in the same ballpark. Optimization targets if needed: distance reuse across XPBD iters (small displacements), kernel fusion, fp16 distance matrix.

### `dist_active_static`

Pairwise squared distance matrix between active (dynamic) and static
(boundary) particles.

```
sib_metal dist_active_static <n_active> <n_static> \
    <active.bin> <static.bin> <out.bin> [iters]
```

- Inputs are little-endian fp32 binary files: `[n*3]` for positions,
  `[n_active * n_static]` for the output matrix (row-major).
- The optional `iters` re-runs the kernel `iters` times against the same
  buffers and prints `ms/iter` on stderr (steady-state, excludes Metal
  startup).

**Why active×static, not all×all:** demo1 has ~343 dynamic + ~17,498
static (boundary) particles. Static-static interactions never fire
because the boundary doesn't move. Computing only active×static drops
the matrix from 17K×17K (1.16 GB) to 343×17K (23 MB).

**Steady-state perf** (Apple Silicon, demo1 scale):

```
n_active=343, n_static=17498  →  0.45 ms/iter
```

For comparison, the full OpenCL step on the same scenario (L4 GPU) is
1.6 ms; the SPH neighbor lookup alone is 0.08 ms. Our distance kernel
is doing more work (full N×M materialization vs sparse list lookup) but
in a form that's directly usable for matrix-formulated SPH downstream.

### `wpoly6_inplace`

Applies Müller 2003 Wpoly6 SPH smoothing kernel elementwise on a squared distance buffer, in-place.

```
sib_metal wpoly6_inplace <n_total> <h> <inout.bin> [iters]
```

Input/output is the same `[n_total]` fp32 buffer. The kernel formula `W = (315/(64πh⁹)) · (h² - r²)³` for `r < h` else 0; the host computes `poly6_const` and `h²`.

### `rowsum_density`

Per-row reduction of an `[n_rows, n_cols]` matrix scaled by particle mass.

```
sib_metal rowsum_density <n_rows> <n_cols> <mass> <W.bin> <density.bin> [iters]
```

Uses threadgroup tree reduction (256 threads/row) — much faster than naive per-row serial sum AND more numerically accurate due to fewer accumulation steps.

### `dist_active_active`

Squared distance matrix between active particles (cube-cube, liquid-liquid).

```
sib_metal dist_active_active <n_active> <active.bin> <out.bin> [iters]
```

Same pattern as `dist_active_static` with one buffer used twice. Diagonal entries are exactly 0 (self-distance); downstream consumers skip i==j.

### `density_constraint_grad`

XPBD's density constraint gradient — fused: sums Wspiky-grad contributions over all neighbors (active + static), produces per-active-particle ∇C (3-vector).

```
sib_metal density_constraint_grad <n_active> <n_static> <h> <mass> <rho_rest> \
    <active.bin> <static.bin> <r2_aa.bin> <r2_as.bin> <gradC.bin> [iters]
```

Uses Wpoly6 for density (separately, M6.1+M6.2) but Wspiky for the gradient — Macklin 2013 PBD-Fluids convention to avoid pressure-clustering pathology. Threadgroup reduction over float3 partials.

## Roadmap (revised after gradient-chain Option 3 + XPBD switch)

- **M6.0** ✓ active×static squared distance kernel
- **M6.1** ✓ Wpoly6 elementwise on distance matrix
- **M6.2** ✓ density via threadgroup-reduced row-sum
- **M6.3** ✓ active×active distance kernel
- **M6.4** ✓ density constraint gradient + denom_helper (fused
            Wspiky-grad + reduction; outputs both ∇C and Σ|∇W|² for
            the proper XPBD denominator)
- **DROPPED** Wvisc-lap kernel — XPBD handles viscosity differently
- **DROPPED** semi-implicit Euler integrator — XPBD owns the integration
- **DROPPED** PyTorch wrapper — Option 3 (paired forward/backward step
              functions) doesn't need general autograd

- **M7.A** ✓ XPBD forward orchestration: predict_positions,
            solve_density_constraint (with proper denominator),
            solve_floor_constraint, update_velocities, add_inplace
            utility, plus xpbd_step driver that runs full step end-to-end
- **M7.B** TODO: distance-constraint solver (elastic bonds) for cube
            cohesion + worm muscle
- **M7.C** TODO: backward kernels paired with each forward kernel
            (gradient chain Option 3 — no general autograd)
- **M7.D** TODO: differentiability validation (numerical-vs-autograd
            gradient check)
- **M7.E** TODO: parameter-learning demo (start with miscalibrated
            rho_rest, learn the right value via SGD)

### `xpbd_step`

Full XPBD timestep on hand-written Metal.

```
sib_metal xpbd_step \
    <n_active> <n_static> <h> <mass> <rho_rest> <dt> <gravity_y> \
    <floor_y> <alpha_density> <n_iters> \
    <pos_active.bin> <vel_active.bin> <pos_static.bin> [bench_steps]
```

Outputs written to `/tmp/xpbd_pos_out.bin` and `/tmp/xpbd_vel_out.bin`.
For multi-step benchmarking, pass `bench_steps` > 1; per-step wall time
prints to stderr.

**Pipeline per step:**
1. `predict_positions` — apply gravity to vel, integrate to x_pred
2. For each XPBD iter:
   a–b. `dist_active_active` + `dist_active_static` (recompute distances
        — particles moved on previous iter)
   c–d. `wpoly6_inplace` on r²_aa and r²_as (via Metal blit copy first)
   e–g. `rowsum_density` × 2 + `add_inplace` to combine
   h.  `density_constraint_grad` (also outputs denom_helper)
   i.  `solve_density_constraint`
   j.  `solve_floor_constraint`
3. `update_velocities` — recover v from x_pred change

**Steady-state perf:**

| Scenario | n_active | n_static | iters | ms/step |
|---|---|---|---|---|
| Small cube + floor plate | 64 | 100 | 3 | 0.67 |
| demo1-like (random) | 343 | 17,498 | 3 | 2.67 |

Reference: OpenCL full demo1 step on L4 = 1.61 ms (with PCISPH). We're
~65% slower on demo1-comparable workload — first-pass implementation.
Optimization headroom: distance reuse across XPBD inner iters (small
displacement assumption), kernel fusion, spatial grid for >2K particles.

**Bug found and fixed during M7.A:** initial implementation had a CPU
memcpy between distance computation and Wpoly6 evaluation. The memcpy
ran on the CPU before the GPU command buffer executed → Wpoly6 read
stale/uninitialized data → density constraint fired with garbage values
→ particles launched at >400 m/s on step 1. Fix: use Metal `blitCommandEncoder`
for GPU-side copy, properly ordered within the command buffer. Same
class of bug we saw in the imperative Taichi solver (different cause).
Lesson: never CPU-touch a Metal shared buffer between queued encoders
without `waitUntilCompleted` first.

## Build

```bash
bash src/metal_diff/build.sh
```

Requires the Apple Metal toolchain (`xcrun metal`, ships with Xcode
Command Line Tools). Verified on Darwin 25.3.0 / Apple Silicon.

## Test

```bash
.venv/bin/python src/metal_diff/test_dist.py
```

Runs two cases — small (correctness vs numpy) and demo1-realistic
(correctness + wall time + numpy comparison).
