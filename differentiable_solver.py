"""Differentiable PCISPH solver — strategic pivot from imperative Taichi.

This file is the differentiable-Taichi successor to taichi_solver.py.
The goal is gradient-based optimization use cases: parameter tuning,
muscle activation learning for the worm, etc.

Strategic context: native Metal (PR #222) and a future native CUDA backend
serve the high-performance non-differentiable case. This file's reason to
exist is autodiff. If autodiff is not needed, use OpenCL or native backends.

MILESTONES (build in order; each one validates before the next starts):
  M1 (this file): differentiable gravity-only fall — single particle drops
       under gravity for N steps; verify dY_final/dV_init = N*dt.
  M2: Multi-particle gravity (1000 particles, gradient flows correctly)
  M3: Elastic spring force (PCISPH-free; springs only)
  M4: SPH density + viscosity (no PCISPH iteration)
  M5: Full PCISPH (deferred; the implicit pressure correction loop is
       non-trivially differentiable — may need fixed-point unroll or
       implicit-function-theorem trick)

Backend: Taichi 1.7.4 autodiff has been verified working on both CPU
and Metal for the physics-simulation patterns used in M1–M3 (single
particle, multi-particle, spring chain). All milestones default to CPU
for cross-platform reproducibility but accept an `arch` parameter.
"""
import taichi as ti


def run_milestone1(n_steps=100, dt=0.01, init_velocity=5.0, gravity=-9.8,
                   arch=None):
    """Differentiable gravity-only fall.

    Returns:
        (final_y, dY_dV0) where dY_dV0 is the gradient of the final Y
        position with respect to the initial velocity.

    For semi-implicit Euler with constant gravity:
        v(t+1) = v(t) + g*dt
        x(t+1) = x(t) + v(t+1)*dt
    Unrolling N steps:
        y_final = y_0 + sum_{i=1..N} v(i) * dt
                = y_0 + N*v_0*dt + g*dt^2 * (N(N+1)/2)
    So  dy_final/dv_0 = N*dt   (analytical reference)
    """
    ti.init(arch=arch if arch is not None else ti.cpu)

    # Time-indexed fields — Taichi's reverse-mode autodiff requires a
    # distinct slot per timestep so the tape can reconstruct dependencies.
    pos = ti.field(dtype=ti.f32, shape=n_steps + 1, needs_grad=True)
    vel = ti.field(dtype=ti.f32, shape=n_steps + 1, needs_grad=True)
    init_vel_param = ti.field(dtype=ti.f32, shape=(), needs_grad=True)
    loss = ti.field(dtype=ti.f32, shape=(), needs_grad=True)

    @ti.kernel
    def initialize():
        pos[0] = 100.0
        vel[0] = init_vel_param[None]

    @ti.kernel
    def integrate(t: ti.i32):
        vel[t + 1] = vel[t] + gravity * dt
        pos[t + 1] = pos[t] + vel[t + 1] * dt

    @ti.kernel
    def compute_loss():
        loss[None] = pos[n_steps]

    init_vel_param[None] = init_velocity

    with ti.ad.Tape(loss=loss):
        initialize()
        for t in range(n_steps):
            integrate(t)
        compute_loss()

    return float(loss[None]), float(init_vel_param.grad[None])


def run_milestone2(
    n_particles=1000,
    n_steps=100,
    dt=0.01,
    gravity=-9.8,
    target_final_y=50.0,
    n_optimization_iters=50,
    learning_rate=100.0,
):
    """Multi-particle differentiable gravity fall + gradient-based optimization.

    N particles each fall under gravity for n_steps. Each particle has its
    own learnable initial velocity. The loss is the squared distance between
    each particle's final Y and a target. We use Taichi's autodiff to
    compute per-particle gradients and run gradient descent.

    Returns:
        (initial_loss, final_loss) — should converge if autodiff is plumbed
        correctly across the particle dimension.
    """
    ti.init(arch=ti.cpu)

    # 2D fields: time-axis × particle-axis. Taichi's autodiff handles this
    # cleanly as long as each timestep writes to a separate slot.
    pos = ti.field(dtype=ti.f32, shape=(n_steps + 1, n_particles), needs_grad=True)
    vel = ti.field(dtype=ti.f32, shape=(n_steps + 1, n_particles), needs_grad=True)
    init_vel_param = ti.field(dtype=ti.f32, shape=n_particles, needs_grad=True)
    loss = ti.field(dtype=ti.f32, shape=(), needs_grad=True)

    @ti.kernel
    def initialize():
        for p in range(n_particles):
            pos[0, p] = 100.0
            vel[0, p] = init_vel_param[p]

    @ti.kernel
    def integrate(t: ti.i32):
        for p in range(n_particles):
            vel[t + 1, p] = vel[t, p] + gravity * dt
            pos[t + 1, p] = pos[t, p] + vel[t + 1, p] * dt

    @ti.kernel
    def compute_loss():
        # Mean squared error vs target — averaged so loss magnitude is
        # independent of n_particles.
        for p in range(n_particles):
            diff = pos[n_steps, p] - target_final_y
            loss[None] += diff * diff / n_particles

    @ti.kernel
    def gradient_step():
        for p in range(n_particles):
            init_vel_param[p] -= learning_rate * init_vel_param.grad[p]

    # Start every particle at zero velocity — they will all undershoot the
    # target and need a positive boost.
    init_vel_param.fill(0.0)

    losses = []
    for it in range(n_optimization_iters):
        loss[None] = 0.0
        with ti.ad.Tape(loss=loss):
            initialize()
            for t in range(n_steps):
                integrate(t)
            compute_loss()
        losses.append(float(loss[None]))
        gradient_step()

    return losses[0], losses[-1]


def run_milestone3(
    n_particles=10,
    n_steps=600,
    dt=0.005,
    gravity=-9.8,
    rest_length=1.0,
    init_stiffness=50.0,
    target_total_length=12.0,
    n_optimization_iters=40,
    learning_rate=2000.0,
):
    """Differentiable spring chain: gradient through particle interaction.

    A chain of N particles is hung from particle 0 (pinned) under gravity.
    Adjacent particles are connected by springs with learnable stiffness.
    The loss measures (chain_total_length - target). We use autodiff to
    learn the stiffness that produces the desired hanging length.

    Validates that gradients flow correctly through:
      - particle-particle force coupling
      - the integrator (semi-implicit Euler)
      - a learnable physical parameter (stiffness, not just initial state)

    Returns:
        (initial_loss, final_loss, learned_stiffness)
    """
    ti.init(arch=ti.cpu)

    pos = ti.field(dtype=ti.f32, shape=(n_steps + 1, n_particles), needs_grad=True)
    vel = ti.field(dtype=ti.f32, shape=(n_steps + 1, n_particles), needs_grad=True)
    stiffness = ti.field(dtype=ti.f32, shape=(), needs_grad=True)
    loss = ti.field(dtype=ti.f32, shape=(), needs_grad=True)

    @ti.kernel
    def initialize():
        for p in range(n_particles):
            # Hang particles vertically below origin, spaced at rest_length.
            pos[0, p] = -float(p) * rest_length
            vel[0, p] = 0.0

    @ti.kernel
    def integrate(t: ti.i32):
        # Single for-loop covers all particles, with p==0 special-cased
        # inline so Taichi autodiff sees a pure loop kernel (mixed loops
        # and bare statements break reverse-mode tape generation).
        for p in range(n_particles):
            if p == 0:
                # Pinned at the top.
                pos[t + 1, p] = pos[t, p]
                vel[t + 1, p] = 0.0
            else:
                # Spring from above-neighbor.
                above_dy = pos[t, p - 1] - pos[t, p]
                above_force = stiffness[None] * (above_dy - rest_length)
                # Spring from below-neighbor.
                below_force = 0.0
                if p < n_particles - 1:
                    below_dy = pos[t, p + 1] - pos[t, p]
                    below_force = stiffness[None] * (below_dy - (-rest_length))
                damping = -2.0 * vel[t, p]
                accel = above_force + below_force + gravity + damping
                vel[t + 1, p] = vel[t, p] + accel * dt
                pos[t + 1, p] = pos[t, p] + vel[t + 1, p] * dt

    @ti.kernel
    def compute_loss():
        # Total chain length = pos[0] - pos[N-1] (top minus bottom; both
        # are negative in our coordinates, so length is positive).
        chain_length = pos[n_steps, 0] - pos[n_steps, n_particles - 1]
        diff = chain_length - target_total_length
        loss[None] = diff * diff

    @ti.kernel
    def gradient_step():
        stiffness[None] -= learning_rate * stiffness.grad[None]

    stiffness[None] = init_stiffness
    initial_loss = None
    for it in range(n_optimization_iters):
        loss[None] = 0.0
        with ti.ad.Tape(loss=loss):
            initialize()
            for t in range(n_steps):
                integrate(t)
            compute_loss()
        if initial_loss is None:
            initial_loss = float(loss[None])
        gradient_step()
        # Keep stiffness positive so the chain doesn't go unstable.
        if stiffness[None] < 1.0:
            stiffness[None] = 1.0

    final_loss = float(loss[None])
    return initial_loss, final_loss, float(stiffness[None])


def run_milestone4(
    cube_n=4,
    n_steps=600,
    dt=0.005,
    gravity=-9.8,
    floor_y=0.0,
    cube_init_y=8.0,
    rest_length=1.0,
    stiffness=200.0,
    damping=2.0,
    floor_repulsion=500.0,
    floor_thickness=0.5,
    arch=None,
    learn=False,
    n_optimization_iters=30,
    learning_rate=2000.0,
    verbose=True,
):
    """3D elastic cube dropped on a floor — non-pancaking target.

    Builds a cube_n × cube_n × cube_n lattice (default 64 particles)
    connected by springs to axis-aligned neighbors. Drops the cube under
    gravity onto a floor at y=floor_y with a soft repulsive boundary.
    All physics is differentiable.

    Two modes:
      learn=False: run with fixed parameters; report cube extent retention
                   (1.0 = perfectly intact, 0.0 = fully pancaked).
      learn=True:  use autodiff to optimize (stiffness, floor_repulsion)
                   pair to maximize cube extent retention. This is the
                   "fix bugs by gradient descent" application — instead
                   of hand-tuning the magic constants in taichi_solver.py,
                   we let the solver learn them.

    Returns:
        learn=False: (initial_extent, final_extent, retention_ratio,
                      wall_seconds)
        learn=True:  (initial_retention, final_retention, learned_stiffness,
                      learned_repulsion)
    """
    ti.init(arch=arch if arch is not None else ti.cpu)

    n_particles = cube_n ** 3

    pos = ti.Vector.field(3, dtype=ti.f32,
                          shape=(n_steps + 1, n_particles), needs_grad=True)
    vel = ti.Vector.field(3, dtype=ti.f32,
                          shape=(n_steps + 1, n_particles), needs_grad=True)
    # Learnable parameters.
    k_stiffness = ti.field(dtype=ti.f32, shape=(), needs_grad=True)
    k_repulsion = ti.field(dtype=ti.f32, shape=(), needs_grad=True)
    loss = ti.field(dtype=ti.f32, shape=(), needs_grad=True)

    # Per-particle initial positions (computed once, host-side).
    import numpy as np
    init_pos_np = np.zeros((n_particles, 3), dtype=np.float32)
    for ix in range(cube_n):
        for iy in range(cube_n):
            for iz in range(cube_n):
                p = ix * cube_n * cube_n + iy * cube_n + iz
                init_pos_np[p, 0] = ix * rest_length
                init_pos_np[p, 1] = cube_init_y + iy * rest_length
                init_pos_np[p, 2] = iz * rest_length
    init_pos = ti.Vector.field(3, dtype=ti.f32, shape=n_particles)
    init_pos.from_numpy(init_pos_np)

    # Spring connectivity: each particle has up to 6 axis-aligned
    # neighbors (-x, +x, -y, +y, -z, +z). We pre-compute neighbor IDs
    # and write -1 for boundary slots that have no neighbor.
    neighbors_np = -np.ones((n_particles, 6), dtype=np.int32)
    for ix in range(cube_n):
        for iy in range(cube_n):
            for iz in range(cube_n):
                p = ix * cube_n * cube_n + iy * cube_n + iz
                if ix > 0:
                    neighbors_np[p, 0] = (ix - 1) * cube_n * cube_n + iy * cube_n + iz
                if ix < cube_n - 1:
                    neighbors_np[p, 1] = (ix + 1) * cube_n * cube_n + iy * cube_n + iz
                if iy > 0:
                    neighbors_np[p, 2] = ix * cube_n * cube_n + (iy - 1) * cube_n + iz
                if iy < cube_n - 1:
                    neighbors_np[p, 3] = ix * cube_n * cube_n + (iy + 1) * cube_n + iz
                if iz > 0:
                    neighbors_np[p, 4] = ix * cube_n * cube_n + iy * cube_n + (iz - 1)
                if iz < cube_n - 1:
                    neighbors_np[p, 5] = ix * cube_n * cube_n + iy * cube_n + (iz + 1)
    neighbors = ti.field(dtype=ti.i32, shape=(n_particles, 6))
    neighbors.from_numpy(neighbors_np)

    @ti.kernel
    def initialize():
        for p in range(n_particles):
            pos[0, p] = init_pos[p]
            vel[0, p] = ti.Vector([0.0, 0.0, 0.0])

    @ti.kernel
    def integrate(t: ti.i32):
        for p in range(n_particles):
            # Spring contributions from up to 6 neighbors.
            spring_force = ti.Vector([0.0, 0.0, 0.0])
            for c in range(6):
                j = neighbors[p, c]
                if j >= 0:
                    delta = pos[t, j] - pos[t, p]
                    dist = delta.norm() + 1e-6
                    direction = delta / dist
                    # Spring force on p toward j when stretched.
                    spring_force += k_stiffness[None] * (dist - rest_length) * direction

            # Floor repulsion — soft penalty when y < floor_y + floor_thickness.
            floor_force = ti.Vector([0.0, 0.0, 0.0])
            penetration = (floor_y + floor_thickness) - pos[t, p][1]
            if penetration > 0.0:
                floor_force[1] = k_repulsion[None] * penetration

            # Linear damping.
            damping_force = -damping * vel[t, p]

            accel = spring_force + floor_force + damping_force
            accel[1] += gravity

            vel[t + 1, p] = vel[t, p] + accel * dt
            pos[t + 1, p] = pos[t, p] + vel[t + 1, p] * dt

    # Compute extent helpers (max - min along each axis at a given time).
    extent = ti.Vector.field(3, dtype=ti.f32, shape=(n_steps + 1), needs_grad=True)
    pmin = ti.Vector.field(3, dtype=ti.f32, shape=(n_steps + 1), needs_grad=True)
    pmax = ti.Vector.field(3, dtype=ti.f32, shape=(n_steps + 1), needs_grad=True)

    # Three split kernels — autodiff requires each kernel to be either
    # all-loop or all-bare-statement, never both. See M3 docstring.
    @ti.kernel
    def init_extent_seed(t: ti.i32):
        pmin[t] = pos[t, 0]
        pmax[t] = pos[t, 0]

    @ti.kernel
    def reduce_extent(t: ti.i32):
        for p in range(n_particles):
            ti.atomic_min(pmin[t][0], pos[t, p][0])
            ti.atomic_min(pmin[t][1], pos[t, p][1])
            ti.atomic_min(pmin[t][2], pos[t, p][2])
            ti.atomic_max(pmax[t][0], pos[t, p][0])
            ti.atomic_max(pmax[t][1], pos[t, p][1])
            ti.atomic_max(pmax[t][2], pos[t, p][2])

    @ti.kernel
    def finalize_extent(t: ti.i32):
        extent[t] = pmax[t] - pmin[t]

    def compute_extent_at(t):
        init_extent_seed(t)
        reduce_extent(t)
        finalize_extent(t)

    @ti.kernel
    def compute_loss():
        # Smooth loss: sum of squared bond-length deviations at the final
        # timestep. A perfectly intact cube has all bonds at rest_length
        # and produces zero loss; pancake → bonds crushed → high loss.
        # Avoids min/max reductions which break autodiff (non-smooth
        # gradient through Taichi atomic_min/max).
        for p in range(n_particles):
            for c in range(6):
                j = neighbors[p, c]
                if j >= 0:
                    delta = pos[n_steps, j] - pos[n_steps, p]
                    dist = delta.norm() + 1e-6
                    dev = dist - rest_length
                    # Each bond counted twice (from both endpoints).
                    loss[None] += 0.5 * dev * dev / float(n_particles)

    @ti.kernel
    def gradient_step():
        k_stiffness[None] -= learning_rate * k_stiffness.grad[None]
        k_repulsion[None] -= learning_rate * k_repulsion.grad[None]

    k_stiffness[None] = stiffness
    k_repulsion[None] = floor_repulsion

    import time

    if not learn:
        # Single forward pass with timing.
        wall_start = time.perf_counter()
        initialize()
        for t in range(n_steps):
            integrate(t)
        compute_extent_at(0)
        compute_extent_at(n_steps)
        ti.sync()
        wall_seconds = time.perf_counter() - wall_start

        e0 = extent[0].to_numpy()
        ef = extent[n_steps].to_numpy()
        # Mean Y position at start vs end — confirms the cube actually fell.
        pos_init_np = pos.to_numpy()[0]
        pos_final_np = pos.to_numpy()[n_steps]
        mean_y_init = float(pos_init_np[:, 1].mean())
        mean_y_final = float(pos_final_np[:, 1].mean())
        min_y_final = float(pos_final_np[:, 1].min())
        # Retention = mean of (final / initial) along each axis, clamped to [0, 1].
        ratios = [
            min(1.0, max(0.0, ef[i] / max(e0[i], 1e-9)))
            for i in range(3)
        ]
        retention = sum(ratios) / 3.0

        if verbose:
            print(f"  Initial extent: ({e0[0]:.2f}, {e0[1]:.2f}, {e0[2]:.2f})")
            print(f"  Final extent:   ({ef[0]:.2f}, {ef[1]:.2f}, {ef[2]:.2f})")
            print(f"  Mean Y:         {mean_y_init:.2f} → {mean_y_final:.2f}  (min Y final = {min_y_final:.2f})")
            print(f"  Retention:      x={ratios[0]:.1%} y={ratios[1]:.1%} z={ratios[2]:.1%}  mean={retention:.1%}")
            print(f"  Wall time:      {wall_seconds:.3f}s for {n_steps} steps "
                  f"({n_steps/wall_seconds:.0f} steps/s)")

        return mean_y_init, mean_y_final, retention, wall_seconds

    # learn=True: optimize (stiffness, repulsion) via gradient descent.
    losses = []
    for it in range(n_optimization_iters):
        loss[None] = 0.0
        # Reset min/max accumulators each iter.
        pmin.fill(0.0)
        pmax.fill(0.0)
        with ti.ad.Tape(loss=loss):
            initialize()
            for t in range(n_steps):
                integrate(t)
            compute_extent_at(n_steps)
            compute_loss()
        losses.append(float(loss[None]))
        gradient_step()
        # Keep parameters positive.
        if k_stiffness[None] < 1.0:
            k_stiffness[None] = 1.0
        if k_repulsion[None] < 1.0:
            k_repulsion[None] = 1.0

    initial_loss = losses[0]
    final_loss = losses[-1]
    return (initial_loss, final_loss,
            float(k_stiffness[None]), float(k_repulsion[None]))


def run_milestone5(
    n_liquid=125,
    n_steps=600,
    dt=0.005,
    gravity=-9.8,
    # Calibrate rest_density to the actual initial density measured from
    # this layout (h=2.0, spacing=h*0.6, particle_mass=1.0). Setting
    # rest=initial means pressure ≈ 0 at t=0 → particles fall under
    # gravity rather than exploding outward from negative pressure.
    rest_density=0.20,
    pressure_k=10.0,
    viscosity_mu=0.5,
    smoothing_h=2.0,
    particle_mass=1.0,
    box_min_y=0.0,
    floor_repulsion=500.0,
    floor_thickness=0.5,
    arch=None,
    verbose=True,
):
    """Differentiable SPH liquid in a box — first apples-to-apples block.

    Implements the standard Müller 2003 SPH formulation:
      density:   ρ_i = Σ_j m_j W_poly6(|r_ij|, h)
      pressure:  p_i = k (ρ_i - ρ_rest)         (equation of state)
      pressure F: F_p_i = -m_i Σ_j m_j (p_i + p_j) / (2 ρ_j) ∇W_spiky(r_ij, h)
      viscosity: F_v_i = μ m_i Σ_j m_j (v_j - v_i) / ρ_j ∇²W_visc(|r_ij|, h)

    Uses ALL-PAIRS neighbor search (O(N²)). Adequate for N≤2000 on Metal.
    For demo1-scale (~17K) we'd need a spatial-grid neighbor search; that's
    out of scope for M5 and will be a separate spatial-acceleration work
    item (M6 likely).

    Differentiability: all SPH kernels and the integrator are smooth
    functions of position/velocity, so autodiff should propagate cleanly
    on CPU. Metal autodiff has the same caveat as M4 (forward fast, train
    on CPU).

    Returns (mean_y_init, mean_y_final, wall_seconds, mean_density_final).
    """
    ti.init(arch=arch if arch is not None else ti.cpu)

    import math
    import time
    import numpy as np

    # SPH kernel constants (Müller 2003).
    pi = math.pi
    h = smoothing_h
    h2 = h * h
    poly6_const = 315.0 / (64.0 * pi * h ** 9)
    spiky_grad_const = -45.0 / (pi * h ** 6)
    visc_lap_const = 45.0 / (pi * h ** 6)

    pos = ti.Vector.field(3, dtype=ti.f32,
                          shape=(n_steps + 1, n_liquid), needs_grad=True)
    vel = ti.Vector.field(3, dtype=ti.f32,
                          shape=(n_steps + 1, n_liquid), needs_grad=True)
    density = ti.field(dtype=ti.f32, shape=n_liquid, needs_grad=True)
    pressure = ti.field(dtype=ti.f32, shape=n_liquid, needs_grad=True)
    accel = ti.Vector.field(3, dtype=ti.f32, shape=n_liquid, needs_grad=True)

    # Initial liquid layout: 5×5×5 cluster (or as many as n_liquid allows)
    init_pos_np = np.zeros((n_liquid, 3), dtype=np.float32)
    side = int(np.ceil(n_liquid ** (1 / 3)))
    spacing = h * 0.6  # close enough to overlap kernel support
    cluster_y = box_min_y + 5.0  # well above floor
    for p in range(n_liquid):
        ix = p % side
        iy = (p // side) % side
        iz = p // (side * side)
        init_pos_np[p, 0] = (ix - side / 2) * spacing
        init_pos_np[p, 1] = cluster_y + iy * spacing
        init_pos_np[p, 2] = (iz - side / 2) * spacing
    init_pos = ti.Vector.field(3, dtype=ti.f32, shape=n_liquid)
    init_pos.from_numpy(init_pos_np)

    @ti.kernel
    def initialize():
        for p in range(n_liquid):
            pos[0, p] = init_pos[p]
            vel[0, p] = ti.Vector([0.0, 0.0, 0.0])

    @ti.kernel
    def compute_density(t: ti.i32):
        for i in range(n_liquid):
            d = 0.0
            for j in range(n_liquid):
                r = (pos[t, i] - pos[t, j]).norm()
                if r < h:
                    diff = h2 - r * r
                    d += particle_mass * poly6_const * diff * diff * diff
            density[i] = d

    @ti.kernel
    def compute_pressure():
        for i in range(n_liquid):
            pressure[i] = pressure_k * (density[i] - rest_density)

    @ti.kernel
    def compute_accel(t: ti.i32):
        for i in range(n_liquid):
            f_pressure = ti.Vector([0.0, 0.0, 0.0])
            f_viscosity = ti.Vector([0.0, 0.0, 0.0])
            for j in range(n_liquid):
                if i != j:
                    rij = pos[t, i] - pos[t, j]
                    r = rij.norm() + 1e-6
                    if r < h:
                        # Pressure force (symmetric, Wspiky gradient).
                        h_minus_r = h - r
                        grad_w = spiky_grad_const * h_minus_r * h_minus_r
                        denom = 2.0 * density[j] + 1e-3
                        f_pressure += -particle_mass * (
                            (pressure[i] + pressure[j]) / denom
                        ) * grad_w * (rij / r)
                        # Viscosity force (Wviscosity laplacian).
                        lap_w = visc_lap_const * h_minus_r
                        f_viscosity += viscosity_mu * particle_mass * (
                            (vel[t, j] - vel[t, i]) / (density[j] + 1e-3)
                        ) * lap_w
            # Floor repulsion.
            f_floor = ti.Vector([0.0, 0.0, 0.0])
            penetration = (box_min_y + floor_thickness) - pos[t, i][1]
            if penetration > 0.0:
                f_floor[1] = floor_repulsion * penetration
            total = f_pressure + f_viscosity + f_floor
            total[1] += gravity * particle_mass
            accel[i] = total / particle_mass

    @ti.kernel
    def integrate(t: ti.i32):
        for p in range(n_liquid):
            vel[t + 1, p] = vel[t, p] + accel[p] * dt
            pos[t + 1, p] = pos[t, p] + vel[t + 1, p] * dt

    initialize()
    wall_start = time.perf_counter()
    for t in range(n_steps):
        compute_density(t)
        compute_pressure()
        compute_accel(t)
        integrate(t)
    ti.sync()
    wall_seconds = time.perf_counter() - wall_start

    pos_init_np = pos.to_numpy()[0]
    pos_final_np = pos.to_numpy()[n_steps]
    mean_y_init = float(pos_init_np[:, 1].mean())
    mean_y_final = float(pos_final_np[:, 1].mean())
    density_np = density.to_numpy()
    mean_density = float(density_np.mean())

    if verbose:
        print(f"  N={n_liquid}, n_steps={n_steps}, h={smoothing_h}")
        print(f"  Mean Y:           {mean_y_init:.2f} → {mean_y_final:.2f}")
        print(f"  Mean density:     {mean_density:.1f}  (rest={rest_density})")
        print(f"  Wall time:        {wall_seconds:.3f}s "
              f"({n_steps/wall_seconds:.0f} steps/s)")

    return mean_y_init, mean_y_final, wall_seconds, mean_density


if __name__ == "__main__":
    # ── Milestone 1 ──
    n_steps = 100
    dt = 0.01
    init_v = 5.0
    final_y, dY_dV = run_milestone1(
        n_steps=n_steps, dt=dt, init_velocity=init_v
    )

    expected_grad = n_steps * dt
    error = abs(dY_dV - expected_grad)

    print(f"=== Milestone 1: differentiable gravity-only fall ===")
    print(f"n_steps={n_steps}, dt={dt}, init_velocity={init_v}")
    print(f"Final Y:           {final_y:.4f}")
    print(f"dY/dV0 (autodiff): {dY_dV:.6f}")
    print(f"dY/dV0 (expected): {expected_grad:.6f}")
    print(f"Error:             {error:.2e}")
    m1_pass = error < 1e-3
    print(f"[{'PASS' if m1_pass else 'FAIL'}] M1: "
          f"{'autodiff matches analytical gradient' if m1_pass else 'diverges from analytical'}")

    # ── Milestone 2 ──
    print()
    print(f"=== Milestone 2: multi-particle gradient descent ===")
    initial_loss, final_loss = run_milestone2()
    print(f"Initial loss: {initial_loss:.4f}")
    print(f"Final loss:   {final_loss:.4f}")
    print(f"Reduction:    {initial_loss / max(final_loss, 1e-9):.1f}x")
    m2_pass = final_loss < initial_loss * 0.01  # 100x reduction
    print(f"[{'PASS' if m2_pass else 'FAIL'}] M2: "
          f"{'gradient descent converges across 1000 particles' if m2_pass else 'failed to converge'}")

    # ── Milestone 3 ──
    print()
    print(f"=== Milestone 3: spring chain (gradient through particle interaction) ===")
    initial_loss, final_loss, learned_k = run_milestone3()
    print(f"Initial loss:       {initial_loss:.4f}")
    print(f"Final loss:         {final_loss:.4f}")
    print(f"Learned stiffness:  {learned_k:.2f}")
    m3_pass = final_loss < initial_loss * 0.5  # 2x reduction is meaningful
    print(f"[{'PASS' if m3_pass else 'FAIL'}] M3: "
          f"{'gradients flow through spring coupling' if m3_pass else 'gradient descent stalled'}")

    # ── Milestone 4 ──
    print()
    print(f"=== Milestone 4: 3D elastic cube on Metal ===")
    print(f"Phase A: forward run with hand-tuned params on arch=metal")
    mean_y_init, mean_y_final, retention, wall = run_milestone4(arch=ti.metal, learn=False)
    m4a_pass = retention > 0.8 and mean_y_final < mean_y_init - 1.0  # cube fell AND retained shape
    print(f"[{'PASS' if m4a_pass else 'FAIL'}] M4-A: "
          f"{'cube fell intact and retained shape' if m4a_pass else 'cube either did not fall or pancaked'}")

    # ── Milestone 4 — performance scan ──
    print()
    print(f"Phase B: performance scan on arch=metal (find OpenCL-comparable size)")
    print(f"  {'cube_n':>8} {'n_particles':>12} {'wall':>10} {'steps/s':>10} {'retention':>10}")
    for cn in [4, 8, 12, 16, 20, 24]:
        my_init, my_final, ret, w = run_milestone4(
            cube_n=cn, n_steps=600, arch=ti.metal, learn=False, verbose=False
        )
        print(f"  {cn:>8d} {cn**3:>12d} {w:>9.3f}s {600/w:>10.0f} {ret:>9.1%}")

    # ── Milestone 4 — autodiff still works on this scenario ──
    print()
    print(f"Phase C: autodiff parameter tuning on 3D cube (arch=cpu)")
    print(f"  (Metal autodiff partially propagates gradients on the 3D")
    print(f"   spring kernel — Taichi 1.7.4 Metal backend limitation. CPU")
    print(f"   autodiff works fully. Production strategy: train on CPU,")
    print(f"   inference/forward on Metal.)")
    init_loss, final_loss, learned_k, learned_rep = run_milestone4(
        cube_n=4, n_steps=600, arch=ti.cpu, learn=True,
        cube_init_y=15.0, damping=0.0,
        stiffness=10.0, floor_repulsion=2000.0,
        n_optimization_iters=20, learning_rate=50.0,
    )
    print(f"  Initial loss: {init_loss:.4f}")
    print(f"  Final loss:   {final_loss:.4f}")
    print(f"  Learned: stiffness={learned_k:.1f}  floor_repulsion={learned_rep:.1f}")
    m4c_pass = final_loss < init_loss * 0.8
    print(f"[{'PASS' if m4c_pass else 'FAIL'}] M4-C: "
          f"{'autodiff tunes 3D cube parameters' if m4c_pass else 'autodiff failed in 3D — cube too robust'}")

    # ── Milestone 5: SPH liquid in a box (Metal performance scan) ──
    print()
    print(f"=== Milestone 5: differentiable SPH liquid (arch=metal) ===")
    print(f"  {'N':>6} {'wall':>10} {'steps/s':>10} {'mean Y init→final':>22} {'mean ρ':>10}")
    for n in [125, 500, 1000, 2000]:
        my_init, my_final, w, rho = run_milestone5(
            n_liquid=n, n_steps=300, arch=ti.metal, verbose=False
        )
        print(f"  {n:>6d} {w:>9.3f}s {300/w:>10.0f} "
              f"{my_init:>10.2f} → {my_final:>6.2f} {rho:>10.1f}")
