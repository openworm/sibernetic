# M7: XPBD plan for differentiable Sibernetic

## Why XPBD over original PBD or PCISPH

Decided 2026-05-03 (energy-conservation conversation). XPBD (Macklin et
al., 2016) gives PBD's benefits — no implicit pressure loop, naturally
differentiable, fast per step — without PBD's stiffness-via-iteration-
count pathology. As iteration count grows, XPBD converges to the
implicit-Euler force-based solution; energy behavior approaches that of
a respectable force-based integrator. PCISPH retained as a future
fallback if biomechanics validation requires it; for now XPBD is
correct enough and differentiable enough.

## Constraints to enforce per timestep

| Constraint | Particles | C(x) = 0 means | Compliance α |
|---|---|---|---|
| Density (incompressibility) | Liquid | ρ_i / ρ_rest - 1 = 0 | ≈ 0 (near-rigid) |
| Distance (elastic bonds) | Elastic-bonded pairs | \|\|p_i - p_j\|\| - rest_len = 0 | 1/elasticity_coefficient (from Sibernetic config) |
| Floor (boundary) | Any (above floor) | p_i.y - floor_y >= 0 | 0 (hard) |

`α = 1/k` where `k` is stiffness. Sibernetic's `elasticity_coefficient = 3e8`
maps directly: α_distance = 1/3e8 ≈ 3.3e-9. Density constraint compliance
controls how compressible the fluid is — start at α_density ≈ 1e-7 for
near-incompressible behavior, tunable.

## Standard XPBD update (per timestep)

```
# 1. Predict positions under external forces (gravity, viscous drag)
v_pred = v + dt * F_ext / m
x_pred = x + dt * v_pred

# 2. Initialize Lagrange multipliers (one per constraint)
λ = zeros(n_constraints)

# 3. Iterate constraint projections (3-10 iterations typical)
for iter in range(n_iters):
    for c in constraints:
        C    = c.value(x_pred)              # constraint violation
        ∇C   = c.gradient(x_pred)           # gradient w.r.t. positions
        # Δλ formula derived from implicit Euler stationary equations:
        Δλ   = -(C + (α / dt²) * λ) / (∇Cᵀ M⁻¹ ∇C + α / dt²)
        Δx   = M⁻¹ ∇C Δλ
        x_pred += Δx
        λ += Δλ

# 4. Recover velocity from position change
v = (x_pred - x) / dt
x = x_pred
```

The `α / dt²` term is the XPBD modification — this is what removes the
timestep-dependent stiffness pathology of original PBD.

## Mapping to M6's matrix-formulated SPH

For the **density constraint**:

```
C_density(x) = (ρ_i(x) / ρ_rest) - 1
ρ_i(x) = Σ_j m_j W_poly6(||p_i - p_j||², h)   ← from M6.1+M6.2
∇C_density / ∂p_i = (1/ρ_rest) Σ_j m_j ∇W_poly6
∇W_poly6 = -gradient of poly6 kernel
```

The matrix structure from M6 (D², W matrix, density vector) is exactly
what XPBD's density constraint needs to evaluate `C` and `∇C`. We reuse
the same kernels.

For the **distance constraint**:

```
C_dist(p_i, p_j) = ||p_i - p_j|| - rest_len
∇C_dist / ∂p_i = (p_i - p_j) / ||p_i - p_j||
∇C_dist / ∂p_j = -(p_i - p_j) / ||p_i - p_j||
```

Per-bond computation. We need an "elastic_bonds" buffer (i, j, rest_len,
α) loaded from Sibernetic's connection format. Use a small Metal kernel
that iterates over all bonds in parallel.

For the **floor constraint** (and any other inequality):

```
C_floor(p_i) = max(0, floor_y - p_i.y)   ← only active when violated
∇C_floor / ∂p_i.y = -1 (when violated)
```

Trivial per-particle kernel.

## Metal kernels to add (M7.x ladder)

| Kernel | Purpose | Inputs | Outputs |
|---|---|---|---|
| M7.0 `predict_positions` | Apply external forces | x, v, F_ext, dt, m | x_pred, v_pred |
| M7.1 `solve_density_constraint` | Project density violation | x_pred, ρ, ∇W matrix, ρ_rest, α, dt, λ | Δx accumulated, Δλ accumulated |
| M7.2 `solve_distance_constraints` | Project bond stretches | x_pred, bonds, rest_lens, α, dt, λ | Δx, Δλ |
| M7.3 `solve_floor_constraints` | Project floor penetration | x_pred, floor_y | Δx |
| M7.4 `update_velocities` | v = (x_pred - x)/dt | x_pred, x_old, dt | v_new |
| M7.5 outer loop driver | Run M7.0 → M7.1×N → M7.2×N → M7.3 → M7.4 | All buffers | Updated x, v |

## Differentiability (M6.6 + M7.x backward)

Each forward kernel above gets a hand-derived backward kernel that
computes `∂L/∂inputs` from `∂L/∂outputs`. XPBD is naturally
differentiable because it has no implicit loop — each iteration is an
explicit operation with a closed-form derivative.

The backward pass for the constraint projection:
- Input gradient: `∂L/∂x_post_projection`
- Output gradient: `∂L/∂x_pre_projection`, `∂L/∂α` (compliance), `∂L/∂rest_len`
- Computation: chain rule through `Δx = M⁻¹ ∇C Δλ` and the `Δλ` formula

Material parameters (compliance, rest length, density rest) become
learnable. This is exactly the "fix bugs by gradient descent" application
— start with miscalibrated rest_density and let autodiff find the value
that minimizes a target metric (e.g., cube extent retention).

## Order of operations

M6 must finish first (need M6.1 = Wpoly6, M6.2 = density, M6.3 =
active-active variant for liquid-liquid, M6.4 = ∇W kernel for force
gradient, M6.5 = integrator infrastructure). Then M7 builds on M6 by
adding constraint projection kernels and the outer iteration loop.

Estimated kernel count: M6 = 8 forward kernels + 8 backward = 16. M7
adds 5 more forward + 5 backward = 10. Total ≈ 26 hand-written kernels
to reach a fully differentiable XPBD demo1 substrate. Each is small
(20-80 lines of Metal); the overall project is significant but tractable.

## Validation milestones for M7

- M7.A: distance constraint only (1D chain) — verify XPBD recovers
  Hooke's-law equilibrium under gravity (this is M3 redone, but in XPBD
  rather than force-based)
- M7.B: density constraint only (liquid in box) — verify fluid is
  incompressible
- M7.C: full demo1 (cube + liquid + boundary, all 3 constraints) —
  cube doesn't pancake, fluid stays bounded, performance comparable to
  OpenCL baseline
- M7.D: differentiability proof — start from miscalibrated parameters,
  use autodiff to learn the true Sibernetic config values
