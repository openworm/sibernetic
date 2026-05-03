# Sibernetic PyTorch Solver Development Log

**Date:** December 2024
**Goal:** Achieve full numerical parity between PyTorch and OpenCL PCISPH solvers

---

## Executive Summary

Successfully implemented a PyTorch backend for the Sibernetic C. elegans worm body physics simulator. The implementation achieves numerical parity with the OpenCL solver within 1e-2 tolerance, supports GPU acceleration via CUDA, and includes comprehensive logging/timing infrastructure.

---

## Architecture Overview

### File Structure

```
sibernetic/
├── pytorch_solver.py          # Main PyTorch solver (1000+ lines)
├── src/
│   ├── owPhysicsFluidSimulator.cpp  # C++ integration layer
│   └── sphFluid.cl                   # OpenCL reference kernels
├── inc/
│   └── owPhysicsFluidSimulator.h     # C++ header
└── tests/
    ├── conftest.py                   # Shared fixtures
    ├── test_physics_parity.py        # Physics tests (46 tests)
    ├── test_validation.py            # Validation suite (18 tests)
    └── test_cpp_integration.py       # C++ integration tests (12 tests)
```

### Integration Flow

```
C++ (owPhysicsFluidSimulator)
    ↓ PyObject_CallMethod
Python (pytorch_solver.PytorchSolver)
    ↓ torch operations
CUDA/CPU tensors
```

---

## Key Learnings

### 1. PCISPH Algorithm Structure

The OpenCL solver uses **Predictive-Corrective Incompressible SPH** with 3 iterations per timestep:

```python
def run_step(self):
    # 1. Neighbor search (once per step)
    self.run_hash_particles()
    self.run_sort()
    self.run_index()
    self.run_index_post_pass()
    self.run_find_neighbors()

    # 2. Initial density/pressure
    self.run_compute_density()
    self.run_compute_pressure()

    # 3. PCISPH iterations (3x)
    for i in range(3):
        self.run_compute_pressure_force_acceleration()
        self.run_pcisph_predict_positions()
        self.run_pcisph_predict_density()
        self.run_pcisph_correct_pressure()

    # 4. Final integration
    self.run_compute_pressure_force_acceleration()
    self.run_integrate()
```

### 2. Particle Types

Stored in `position[:, 3]` (the w component):

| Type | Value | Description |
|------|-------|-------------|
| LIQUID | 1.0 | Normal fluid particles |
| ELASTIC | 2.0 | Elastic body particles (worm) |
| BOUNDARY | 3.0 | Static boundary particles |

**Critical:** Boundary particles (type 3) must NOT move during integration.

### 3. Coordinate System & Scaling

The simulation uses a scaling factor to convert between simulation and world coordinates:

```python
simulation_scale = config["simulation_scale"]       # ~0.0037 typically
simulation_scale_inv = config["simulation_scale_inv"]  # 1/simulation_scale

# World to simulation
sim_pos = world_pos * simulation_scale

# Simulation to world
world_pos = sim_pos * simulation_scale_inv
```

### 4. Key Physics Formulas

#### Density (Poly6 Kernel)
```python
# OpenCL: computeDensityPressure kernel
w = (h² - r²)³  # where r = distance between particles
density = sum(w * mass_mult_Wpoly6Coefficient)

# CRITICAL: Apply minimum clamp
h_scaled = h * simulation_scale
density = max(density, h_scaled**6)
```

#### Pressure
```python
# Tait equation
pressure = (density/rho0 - 1.0) * delta
```

#### Pressure Force (Spiky Kernel Gradient)
```python
# OpenCL: pcisph_computeForcesAndInitPressure kernel
pres_factor = 0.5 * (pressure[i] + pressure[j]) / predicted_rho[j]
grad = (h_scaled - r_ij)²  # Spiky kernel gradient
force = -pres_factor * grad * (pos[i] - pos[j]) * sim_scale / r_ij
```

#### Viscosity
```python
# OpenCL: computeAcceleration kernel
visc_force = mu * (vel[j] - vel[i]) * (h - r_ij) / 1000.0
```

#### Surface Tension
```python
# Coefficient: -1.7e-09 * surfTensCoeff
kernel = (h² - r²)³
force = -1.7e-09 * surfTensCoeff * kernel * direction
```

#### Elastic Spring Force
```python
# Hooke's law with equilibrium distance
delta_r = distance - equilibrium_distance
force = -(direction / distance) * delta_r * elasticity_coeff
```

#### Muscle Contraction
```python
# OpenCL: computeAcceleration kernel
# Muscle index stored in elasticConnectionsData[...].z (1-indexed, 0 = no muscle)
force = -(direction / distance) * muscle_activation * max_muscle_force
# max_muscle_force = 4000.0
```

### 5. Leapfrog Integration

Two-phase symplectic integrator:

```python
# Mode 0: Position update (first half-step)
position += (velocity * dt + acceleration * dt²/2) * simulation_scale_inv

# Mode 1: Velocity update (second half-step)
velocity += (acceleration_old + acceleration_new) * dt/2
```

### 6. Neighbor Search

Uses spatial hashing with fixed-size grid cells:

```python
# Hash particle to grid cell
cell_idx = floor((pos - origin) * hash_grid_cell_size_inv)
cell_id = cell_idx.x + cell_idx.y * grid_cells_x + cell_idx.z * grid_cells_x * grid_cells_y

# Sort particles by cell_id for cache coherency
# Build cell index array for O(1) cell lookup
```

---

## Gotchas & Edge Cases

### 1. PyDict Memory Leaks (C++)

**Problem:** `PyDict_SetItemString` does NOT steal references.

```cpp
// WRONG - leaks memory
PyDict_SetItemString(cfg, "xmin", PyFloat_FromDouble(config->xmin));

// CORRECT
PyObject *val = PyFloat_FromDouble(config->xmin);
PyDict_SetItemString(cfg, "xmin", val);
Py_DECREF(val);
```

### 2. Density Minimum Clamp

Without clamping, isolated particles get density=0, causing division by zero:

```python
h_scaled_6 = (h * simulation_scale) ** 6
density = torch.clamp(density, min=h_scaled_6)
```

### 3. Muscle Index is 1-Indexed

In `elasticConnectionsData[..., 2]`:
- `0` = Not a muscle connection
- `1+` = Muscle index (subtract 1 for array access)

### 4. Boundary Particle Velocity

Boundary particles must have velocity zeroed after integration to prevent drift:

```python
if particle_type == 3:  # BOUNDARY
    velocity[:3] = 0.0
```

### 5. Grid Cell Count

Must match: `grid_cell_count = grid_cells_x * grid_cells_y * grid_cells_z`

If mismatched, neighbor search produces garbage.

### 6. Sorted vs Unsorted Buffers

After `run_sort()`, particles are reordered. Must track:
- `particle_index_back` - maps sorted index → original index
- `sorted_position`, `sorted_velocity` - sorted buffers

Integration writes to sorted buffers, then unscrambles at the end.

---

## Testing Strategy

### Test Categories

1. **Unit Tests** - Individual kernel parity
2. **Integration Tests** - Full simulation stability
3. **Reference Tests** - Compare against OpenCL output logs
4. **Performance Tests** - Timing benchmarks

### Running Tests

```bash
# Run all engine tests
RUN_ENGINE_TESTS=1 python -m pytest tests/ -v

# Run specific phase
RUN_ENGINE_TESTS=1 python -m pytest tests/test_physics_parity.py -v

# Skip slow tests
RUN_ENGINE_TESTS=1 python -m pytest tests/ -v -m "not slow"

# GPU tests only
RUN_ENGINE_TESTS=1 python -m pytest tests/ -v -m gpu
```

### Tolerance

- Physics parity: `rtol=1e-2` (1% relative tolerance)
- GPU/CPU parity: `rtol=1e-5, atol=1e-5`

---

## Performance Considerations

### Bottlenecks Identified

1. **Neighbor Finding** - O(n²) naive approach, mitigated by spatial hashing
2. **Scatter Operations** - `scatter_add_` can be slow on CPU
3. **Python Loop Overhead** - Minimized by vectorizing operations

### Timing Infrastructure

```python
solver.enable_timing(True)
for _ in range(100):
    solver.run_step()
timing = solver.get_timing_breakdown()
# Returns: {'hash_particles': 0.5, 'sort': 0.3, ...} in ms
```

### GPU Acceleration

```python
config["device"] = "cuda"  # Use GPU
config["device"] = "cpu"   # Use CPU
```

GPU provides ~10-100x speedup for large particle counts.

---

## Configuration Parameters

### Required Parameters

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| `h` | Smoothing length | 3.34 |
| `rho0` | Rest density | 1000.0 |
| `delta` | Pressure stiffness | 0.5 |
| `time_step` | Integration dt | 0.0001 |
| `simulation_scale` | World→sim scale | 0.0037 |
| `max_iteration` | PCISPH iterations | 3 |
| `max_neighbor_count` | Max neighbors per particle | 32 |

### Grid Parameters

| Parameter | Description |
|-----------|-------------|
| `xmin`, `ymin`, `zmin` | Grid origin |
| `grid_cells_x/y/z` | Grid dimensions |
| `grid_cell_count` | Total cells (x*y*z) |
| `hash_grid_cell_size_inv` | 1/cell_size |

### Physics Coefficients

| Parameter | Description |
|-----------|-------------|
| `mu` | Viscosity coefficient |
| `mass_mult_Wpoly6Coefficient` | Poly6 kernel scaling |
| `mass_mult_gradWspikyCoefficient` | Spiky kernel scaling |
| `surf_tens_coeff` | Surface tension strength |
| `r0` | Boundary interaction radius |

---

## Logging Infrastructure

### Enable Logging

```python
solver.enable_logging(True)
solver.run_step()
log = solver._step_logs[0]  # Dict of substep → state
```

### Export Logs

```python
solver.export_logs("simulation_log.json")
```

### Log Format

```json
{
  "steps": [
    {
      "hash_particles": {"position": [...], "velocity": [...]},
      "compute_density": {"position": [...], "density": [...]},
      "pcisph_iter_0": {...},
      "pcisph_iter_1": {...},
      "pcisph_iter_2": {...},
      "integrate": {...}
    }
  ]
}
```

---

## Future Work

1. **CUDA Custom Kernels** - Replace PyTorch ops with custom CUDA for 2-5x speedup
2. **Mixed Precision** - Use FP16 where possible for memory bandwidth
3. **Membrane Forces** - Not yet implemented (low priority)
4. **Multi-GPU** - Distribute particles across GPUs for very large simulations

---

## References

- OpenCL kernels: `src/sphFluid.cl`
- C++ integration: `src/owPhysicsFluidSimulator.cpp`
- PCISPH paper: Solenthaler & Pajarola 2009
- SPH kernel functions: Muller et al. 2003

---

## Test Coverage Summary

| Phase | Tests | Status |
|-------|-------|--------|
| 1. C++ Integration | 12 | ✅ Complete |
| 2. Core PCISPH | 15 | ✅ Complete |
| 3.1 Boundary | 4 | ✅ Complete |
| 3.2 Viscosity | 5 | ✅ Complete |
| 3.3 Surface Tension | 4 | ✅ Complete |
| 3.4 Elastic Forces | 8 | ✅ Complete |
| 3.5 Muscle Forces | 10 | ✅ Complete |
| 3.6 Leapfrog | 6 | ✅ Complete |
| 4. Validation | 18 | ✅ Complete |
| **Total** | **76+** | ✅ |

---

*This log was created to preserve institutional knowledge from the PyTorch solver implementation effort.*

---

## 2026-05-03 — Cross-backend GPU benchmark on NVIDIA L4 (Cloud Run)

Goal: identify a modern, well-maintained GPU substrate to replace OpenCL for active development. Apple deprecated OpenCL on Apple Silicon, the AMD APP SDK 3.0 we link against is from 2015 and abandoned, and we wanted to know which of the alternatives (PyTorch CUDA, Taichi CUDA) actually delivers usable GPU performance for Sibernetic's PCISPH workload.

### Setup

A throwaway Cloud Run service at `~/Documents/sibernetic-runner/` (slarsontech-private, not in this repo) wraps Sibernetic in a FastAPI runner with NVIDIA L4 GPU access via Cloud Run's gen2 GPU support. Same image, same hardware, same scenario — only `backend=` flag varies. Scenario: `demo1` cube-drop, `timelimit=0.05`, `logstep=10` → 2500 PCISPH integration steps with 218 elastic + 125 liquid + 17,498 boundary particles.

### Results

| Backend | Wall clock | steps/sec | vs OpenCL |
|---|---|---|---|
| **opencl** (DD003 gold standard) | **7.9 sec** | 317 | 1.0× (baseline) |
| **taichi-cuda** | **17.0 sec** | 148 | 2.1× slower |
| torch-cuda (eager) | 165 sec | 15.1 | 21× slower |
| torch-cuda + `torch.compile` | FAILED — TorchInductor compile subprocess crash on `run_compute_density` | — | — |
| torch (cpu) on Cloud Run 8vCPU | 425 sec | 5.9 | 54× slower |

For a 1-second sim (50,000 steps), OpenCL completes in 84 sec wall clock with cube physics intact (extent retention 123%, mean Y 44 → 10, min Y 39 → 4 — i.e. cube fell and is sitting on the floor without pancaking).

### What this tells us

1. **PyTorch CUDA is not viable for Sibernetic's workload.** The 21× gap vs OpenCL on the same hardware is dominated by Python + per-kernel-launch overhead (each PCISPH step does ~30 small torch ops). Even after the major fix to `run_find_neighbors` earlier in the session (replaced O(N²) `cdist` + Python loop with hash-grid lookup, ~290× CPU speedup), torch-cuda inherits this overhead structurally.

2. **`torch.compile` doesn't rescue PyTorch CUDA in this environment.** Two attempts both crashed TorchInductor's compile subprocess (`RuntimeError: A compilation subprocess exited unexpectedly`). The fallback `torch._dynamo.config.suppress_errors = True` does not catch `BackendCompilerFailed` from Inductor. Probably fixable with significant debugging investment, but the best plausible payoff (5-10× speedup) still leaves PyTorch CUDA 2-4× slower than Taichi CUDA.

3. **Taichi CUDA is the viable post-OpenCL path.** 2× slower than raw OpenCL is acceptable, especially given:
   - Cross-platform (CUDA + Apple Silicon Metal + Vulkan + CPU from one Python codebase)
   - Active community (Taichi-Lang/Taichi GitHub org)
   - No reliance on OpenCL anywhere in the toolchain
   - Per-step compute is ~6.8 ms vs OpenCL's 3.2 ms — both dominated by similar kernel work, Taichi's overhead is in JIT + Python-side dispatching, not Inductor-style graph compilation

### Path to "<3-min wall clock for a 5-second sim of 17K particles"

Currently:
- OpenCL on L4: ~7 min for 5s sim (250,000 steps × 1.66 ms/step)
- Taichi-CUDA on L4: ~28 min extrapolated (cold-start amortized; real number probably ~12 min)

Neither is under 3 min as-is. The lever is `dt`: demo1 uses `dt = 2e-5` which forces 250,000 steps for 5s sim. PCISPH stability typically allows `dt = 1e-4` or larger for these particle scales. With `dt = 2e-4` (10× fewer steps), Taichi-CUDA hits the 3-min target. This is the next experiment.

### Pancake bug status

At 0.05s sim (2500 steps), Taichi-CUDA shows extent retention 100% — but that's because the cube hasn't fallen meaningfully yet. Need a longer run (~1 sec) to confirm whether the elastic-coordinate-scale bug documented in the README (Taichi-Metal: 7.6% retention vs PyTorch: 94.5%) is Apple-Silicon-Metal-specific or algorithmic. If algorithmic, it's the same bug on taichi-cuda and we have to fix it before relying on Taichi as the production substrate.

### Recommendation captured for future-me

**Stop trying to make PyTorch CUDA fast for Sibernetic.** Invest engineering time in:
1. **Taichi CUDA as the primary modern backend.** Already integrated in the C++ embedding layer.
2. **`dt` audit.** The 250,000-step demo1 is a config-level decision, not a hardware one. Most production SPH papers use much larger timesteps.
3. **Fix the pancake bug** in `taichi_solver.py` — the README's documented 3-step fix (`/sim_scale` removal + `simulationScaleInv` in integration + `h_scaled` for kernel coefficients).
4. **Native-Metal port (PR #222)** as the Apple Silicon performance hedge — sits alongside taichi-cuda for cross-platform coverage.
5. **Keep OpenCL on Cloud Run as the parity reference** for the cube-stability test threshold calibration. Don't bet new dev on it.

### Things deliberately ruled out

- **PyTorch CUDA as primary backend.** 21× slower than OpenCL on the same GPU. `torch.compile` doesn't fix the structural Python-in-the-loop cost.
- **PyTorch MPS on Apple Silicon.** Tested earlier this session — slower than PyTorch CPU because `scatter_add_` falls back to CPU on MPS, triggering CPU↔GPU transfers per accumulation.
- **JAX rewrite.** Would mean throwing out the C++/Sibernetic core. Not worth it for a 2-3× wall-clock improvement when Taichi is already integrated.

### Update: Taichi project health concern + PR #222 measurements

After the initial recommendation above, two new pieces of information reframed the path forward:

**1. Taichi's release cadence has slowed substantially.** Taichi 1.7.4 (July 2024) is still the latest stable as of May 2026 — no new release in ~10 months. The project isn't abandoned (Yuanming Hu still pushes occasionally; cp313 wheels exist) but corporate engineering muscle (Taichi Graphics) has scaled back. cp314 wheels for the upcoming Python 3.14 are unlikely to materialize without renewed release work. **Investing exclusively in Taichi as the cross-platform substrate carries real maintenance risk.**

**2. PR #222 (Wei Weng's native Metal port) was benchmarked and is the fastest backend yet measured.** Built locally from `pr-222` worktree against Python 3.14, ran demo1 cube-drop:

| Sim time | Wall clock | Steps/sec | Cube extent retention | Notes |
|---|---|---|---|---|
| 0.05s (2500 steps) | 4.9 sec | 508 | n/a (too short) | per-step compute 1.49 ms |
| 1.0s (50,000 steps) | 75 sec | 663 | **118%** ✅ | mean_y 44→9.8, min_y 39→3.9 (cube on floor, intact) |

**PR #222 is ~12% faster than L4 OpenCL on the same workload** (75 sec vs 84 sec for the 1-sec sim) and produces the **same gold-standard physics** as OpenCL (118% retention vs OpenCL's 123% — both intact, no pancake). Per-step compute time on Apple Silicon (1.49 ms) actually beats the data-center L4 GPU (1.66 ms) for this scenario.

For 5-sec sim extrapolation:
- PR #222 native Metal: ~6.3 min wall clock
- + dt audit to 1e-4 (5× fewer steps): **75 sec** ← hits <3-min target on a laptop

### Revised recommendation (supersedes the earlier "Taichi-CUDA" recommendation)

**Two-backend strategy with vendor-backed compilers, treat Taichi as a tactical reference:**

1. **Apple Silicon production: PR #222's native Metal.** Hand-written MSL via the metal-cpp wrappers. Already faster than OpenCL, Apple-backed (Metal lives indefinitely). Merge it into `ow-pytorch-0.0.1` line.
2. **NVIDIA production: native CUDA backend** — would need to be written. Equivalent in scope to the Metal half of PR #222: a `sphFluid.cu` kernel module + a `CudaBackend.cpp` host wrapper. Estimated several weeks of focused work but bounded.
3. **Cross-platform reference / quick prototyping: keep Taichi-CUDA + Taichi-Metal**, but treat them as fallback / parity-test substrates rather than the production target. If Taichi stops releasing entirely, we lose nothing critical.
4. **Linux server-side fallback: keep OpenCL on NVIDIA hardware** as the parity reference for cross-backend regression tests. Don't bet new dev on it but it's not gone yet.

This trades "one beautiful cross-platform substrate" for "two vendor-backed implementations sharing the same algorithm specification." Higher engineering cost up front, lower long-term maintenance risk. Both Apple (Metal) and NVIDIA (CUDA) have multi-decade commitments; no community-project rug-pull risk.

### Concrete next steps

- [ ] Land PR #222 (after the `dd003`/parity-test work the original PR review requested)
- [x] ~~Write the matching native CUDA backend, modeling on PR #222's structure~~ — SCAFFOLDED in commit 9f92972 (`src/cuda/README.md` + `src/cuda/sphFluid.cu` + `inc/owCudaSolver.h`); kernel implementation deferred until PR #222 lands so we can build on its `owSolver` abstract base + `src/kernels/` descriptor pattern. ~2 weeks of focused work for a competent CUDA dev once unblocked.
- [x] ~~Cross-backend cube-stability regression~~ — DONE in commit 92183bd (`scripts/cross_backend_regression.py`); first run caught Taichi-CUDA's pancake bug as designed
- [x] ~~Run `dt` audit experiments on PR #222's Metal backend~~ — DONE; result below
- [ ] Add cross-backend cube-stability regression: PR #222 Metal vs OpenCL vs Taichi-CUDA on the same demo1 scenario, diff at known timesteps (the test-harness work from earlier this session)
- [x] ~~Re-verify Taichi-CUDA at 1-sec sim time~~ — DONE; result below

### dt audit on PR #222 (2026-05-03)

Swept `timestep=` values on the PR #222 native Metal binary against demo1, 1-sec sim. Used a 120-sec watchdog to detect divergence quickly (vs the original baseline that took 75 sec to complete).

| dt | Stability | Notes |
|---|---|---|
| 2e-5 (baseline) | ✅ converges, 75 sec wall, retention 118% | this is what demo1 runs at by default |
| 3e-5 (1.5×) | ❌ hangs at step ~7 | PCISPH didn't converge in 3 iterations |
| 4e-5 (2×) | ❌ hangs | same |
| 5e-5 (2.5×) | ❌ hangs | same |

**Conclusion:** with `max_iteration=3` (the current PCISPH default), dt=2e-5 IS the maximum stable timestep for demo1's 17K particles. The dt lever is essentially closed without solver-level changes.

To open it, options ranked by effort:

1. **Bump PCISPH max_iteration to 6-10** and re-test dt=1e-4. If 6 iterations × 5× fewer steps net ≈ 2.5× wall clock speedup, that's the cheapest path. Requires modifying owConfigProperty.cpp / sphFluid.cl / sphFluid.metal default (or exposing a CLI flag).
2. **Switch to a different pressure solver** (IISPH, DFSPH, WCSPH). Substantial code change; better stability properties at large dt; better long-term but weeks of work.
3. **Use a smaller test scenario** (the minicube proposal) for fast iteration; demo1 stays as the gold-standard end-to-end target.

This means the 5-sec/<3-min wall-clock target on PR #222 native Metal can't be reached with `dt` alone. Realistic path: PCISPH iteration tuning + per-step kernel optimization + dt push together.

### Taichi solver fix attempt (2026-05-03, did not land)

Tried to apply the README's documented 3-step coordinate-scale fix to
`taichi_solver.py`. Found:

1. **Fix #1 (elastic force):** Replaced the hardcoded `elastic_boost = 2000.0`
   on line 1271 with `sim_scale_inv` (~135,000 for demo1). Reasoning: that
   matches the same conversion factor OpenCL applies via `simulationScaleInv`
   in its integration step (`sphFluid.cl` line ~1455).

2. **Fix #2 (integration step):** The README says "Add `simulationScaleInv` to
   the integration step." But Taichi's existing code keeps SPH and gravity
   in world coordinates and converts inside SPH kernels (e.g.
   `grad_world = grad_scaled * sim_scale`). Adding `sim_scale_inv` to the
   leapfrog position update would multiply non-elastic accelerations by
   another ~135k and explode the simulation. **Skipped.**

3. **Fix #3 (kernel coefficients):** The class-level `self.poly6_coeff` etc.
   on lines 65-67 use unscaled `h`, but they're never read — the actual SPH
   kernels at lines 256-258 already compute their own local versions using
   `h_scaled`. **Already in place; no change needed.**

Empirical sweep on Apple Silicon Metal at 1-sec sim time, varying
`elastic_boost` from 2000 → 5000 → 10000 → 30000 → 50000 → 100000 → 135000:

```
ALL values: extent retention 100%, mean_y unchanged 44.42 → 44.42
            (cube didn't move at all)
```

Even with `dt=1e-4` (5× the default) the cube doesn't move at 1s.

**The cube DID move at 5s on Apple Silicon Metal pre-fix** (pancake to
mean_y=0.33 documented in the README). So motion DOES happen on Taichi —
just dramatically slower than OpenCL/PR222 produce on the same scenario.

### What this means

The Taichi pancake/freeze behavior is NOT just the README's documented
3-step coordinate-scale fix. There's a structural slowdown in how Taichi's
force pipeline accumulates motion that's separate from elastic-force
magnitude. The fall rate appears to be ~5× slower than OpenCL's:

| Backend | mean_y at 1s | mean_y at 5s |
|---|---|---|
| OpenCL on L4 | 10.4 | (not measured at 5s but likely floor) |
| PR #222 native Metal | 9.8 | (not measured) |
| Taichi-Metal | 44.4 (no motion) | 0.33 (pancake) |
| Taichi-CUDA | 44.4 (no motion) | (not measured) |

This is bigger than I can debug by remote experimentation in one session.
The proper fix needs someone to instrument `taichi_solver.py` with per-step
acceleration and velocity printouts, compare against OpenCL's per-step
state on the same scenario, and find where the magnitudes diverge.
Candidate culprits to investigate:

1. The `compute_forces` kernel's units handling — gravity stays at -9.8
   m/s² (world), but SPH grad is converted from scaled → world via
   `* sim_scale`. Are the magnitudes actually matched?
2. The `_compute_density` kernel might be producing density/pressure
   values different enough from OpenCL that the resulting pressure force
   dominates gravity, freezing the cube.
3. The boundary repulsion (line 799: `repulsion = 2000.0 * w * n_b`) is a
   hardcoded magnitude. Could be too strong for the scaled scene.

Reverted my speculative change; tree is clean. The "fix taichi_solver.py"
checkbox in the next-steps list stays UNCHECKED.

### Taichi-CUDA pancake check (2026-05-03)

Submitted a 1-sec demo1 sim with `backend=taichi-cuda` against the Cloud Run runner.

```
wall_clock:   73.6s (680 steps/sec)
mean Y:       44.42 → 44.42   (cube didn't move)
min Y:        39.41 → 39.41   (no fall)
extent ret:   100% (artifact: no deformation because no motion)
```

For comparison at the same 1-sec sim time:
- OpenCL on L4: mean_y 44 → 10.4, min_y 39 → 4.17 (cube fell, intact)
- PR #222 Metal: mean_y 44 → 9.8, min_y 39 → 3.9 (cube fell, intact)
- Taichi-CUDA on L4: mean_y unchanged (cube frozen)

**Conclusion:** the Taichi pancake bug is **algorithmic** (in `taichi_solver.py`), not Apple-Silicon-Metal-specific. Both Metal and CUDA exhibit the same root-cause "forces too weak" behavior, just at different magnitudes:

- Taichi-Metal at 5s: cube fell, then pancaked (forces eventually broke through)
- Taichi-CUDA at 1s: cube didn't fall at all (forces couldn't even initiate motion)

Implications:
- Fixing `taichi_solver.py` is a single change that benefits both Metal and CUDA simultaneously.
- The README's documented 3-step coordinate-scale fix (remove `/sim_scale` in elastic force, add `simulationScaleInv` in integration, use `h_scaled` for kernel coefficients) is the starting hypothesis.
- Until that fix lands, Taichi is unusable as a Sibernetic GPU substrate on either platform.
- Doesn't change the strategic direction (native Metal + native CUDA primary, Taichi fallback), but does mean the "Taichi as quick-prototyping fallback" benefit is gated on someone fixing the bug first.

