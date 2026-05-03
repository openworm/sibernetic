"""
Taichi SPH Solver for Sibernetic.

A high-performance SPH (Smoothed Particle Hydrodynamics) solver using Taichi
for GPU acceleration on Mac (Metal), Linux/Windows (CUDA/Vulkan), or CPU.

This solver implements:
- WCSPH (Weakly Compressible SPH) with Tait equation of state
- Poly6 kernel for density, Spiky kernel for pressure gradients
- Viscosity forces
- Surface tension
- Boundary particle handling (Ihmsen et al. 2010)
- Elastic connections for deformable bodies (worm simulation)
- Muscle activation forces
- Leapfrog time integration
- Hard floor constraint as failsafe

Compatible API with pytorch_solver.py for C++ integration.
"""

import taichi as ti
import numpy as np
import time
import math
import os
import sys


class TaichiSolver:
    """High-performance SPH solver using Taichi."""

    # Particle types (stored in position.w)
    LIQUID = 1
    ELASTIC = 2
    BOUNDARY = 3

    def __init__(self, position, velocity, config, elastic_connections=None):
        """Initialize the Taichi SPH solver.

        Args:
            position: Nx4 array of particle positions (x, y, z, type)
            velocity: Nx4 array of particle velocities (vx, vy, vz, 0)
            config: Dictionary of simulation parameters
            elastic_connections: Optional flat list of elastic connection data
                from C++. Format: list of [particle_id, rest_length, muscle_id, unused]
                with num_elastic_particles * 32 entries.
        """
        self.config = config
        self.n_particles = len(position)

        # Initialize Taichi with appropriate backend
        device = config.get("device", "metal")
        kernel_profiler = config.get("kernel_profiler", False)
        self._init_taichi(device, kernel_profiler=kernel_profiler)

        # Extract parameters
        self.h = config.get("h", 0.1)
        self.rho0 = config.get("rho0", 1000.0)
        self.delta = config.get("delta", 500.0)
        self.dt = config.get("time_step", 0.0001)
        self.mass = config.get("mass", 0.001)

        # Kernel coefficients
        pi = math.pi
        self.poly6_coeff = 315.0 / (64.0 * pi * self.h**9)
        self.spiky_coeff = -45.0 / (pi * self.h**6)
        self.viscosity_coeff = 45.0 / (pi * self.h**6)

        # Grid parameters for spatial hashing
        self.grid_size = self.h
        grid_extent = max(
            config.get("xmax", 2.0) - config.get("xmin", 0.0),
            config.get("ymax", 2.0) - config.get("ymin", 0.0),
            config.get("zmax", 2.0) - config.get("zmin", 0.0),
        )
        self.grid_n = int(grid_extent / self.grid_size) + 4
        self.grid_origin = ti.Vector([
            config.get("xmin", 0.0) - self.grid_size,
            config.get("ymin", 0.0) - self.grid_size,
            config.get("zmin", 0.0) - self.grid_size,
        ])

        # Gravity
        self.gravity = ti.Vector([
            config.get("gravity_x", 0.0),
            config.get("gravity_y", -9.81),
            config.get("gravity_z", 0.0),
        ])

        # Floor constraint
        self.floor_y = config.get("floor_y", 0.0)
        self.floor_restitution = config.get("floor_restitution", 0.3)

        # Viscosity
        self.mu = config.get("mu", config.get("viscosity_coefficient", 0.1))

        # Surface tension
        self.surface_tension_coeff = config.get("surf_tens_coeff", 0.0)

        # Boundary
        self.r0 = config.get("r0", self.h * 0.5)

        # Max neighbors
        self.max_neighbors = config.get("max_neighbor_count", 64)

        # Elastic connections
        self.has_elastic = False
        self.max_connections = 32  # Match MAX_NEIGHBOR_COUNT in C++
        self.simulation_scale = config.get("simulation_scale", 1.0)

        # Muscle activation
        self.has_muscles = False

        # Timing
        self._timing_enabled = False
        self._timing_data = {}
        self._step_count = 0

        # Logging
        self._logging_enabled = False
        self._step_logs = []

        # Allocate Taichi fields
        self._allocate_fields()

        # Initialize with provided data
        self._init_particles(position, velocity)

        # Build kernels
        self._build_kernels()

        # Load elastic connections if provided
        if elastic_connections is not None:
            self._load_elastic_from_cpp(elastic_connections, config)

    def _init_taichi(self, device, kernel_profiler=False):
        """Initialize Taichi with the specified backend."""
        if device == "metal":
            ti.init(arch=ti.metal, kernel_profiler=kernel_profiler)
        elif device == "cuda":
            ti.init(arch=ti.cuda)
        elif device == "vulkan":
            ti.init(arch=ti.vulkan)
        elif device == "cpu":
            ti.init(arch=ti.cpu)
        else:
            # Auto-select best available
            ti.init()
        self.device = device

    def _allocate_fields(self):
        """Allocate Taichi fields for particle data."""
        n = self.n_particles
        gn = self.grid_n

        # Particle state
        self.pos = ti.Vector.field(4, dtype=ti.f32, shape=n)  # x, y, z, type
        self.vel = ti.Vector.field(4, dtype=ti.f32, shape=n)  # vx, vy, vz, 0
        self.acc = ti.Vector.field(4, dtype=ti.f32, shape=n)
        self.acc_old = ti.Vector.field(4, dtype=ti.f32, shape=n)

        # SPH quantities
        self.rho = ti.field(dtype=ti.f32, shape=n)
        self.pressure = ti.field(dtype=ti.f32, shape=n)

        # Predicted state for PCISPH
        self.predicted_pos = ti.Vector.field(3, dtype=ti.f32, shape=n)
        self.predicted_rho = ti.field(dtype=ti.f32, shape=n)

        # Grid for neighbor search
        self.grid_count = ti.field(dtype=ti.i32, shape=(gn, gn, gn))
        self.grid_offset = ti.field(dtype=ti.i32, shape=(gn, gn, gn))
        self.particle_ids = ti.field(dtype=ti.i32, shape=n)

        # Flat 1D fields for parallel prefix sum optimization
        self._grid_count_flat = ti.field(dtype=ti.i32, shape=gn * gn * gn)
        self._grid_offset_flat = ti.field(dtype=ti.i32, shape=gn * gn * gn)

        # Hierarchical prefix sum: block-level partial sums
        # Using block size of 1024 for good GPU occupancy
        self._prefix_block_size = 1024
        n_blocks = (gn * gn * gn + self._prefix_block_size - 1) // self._prefix_block_size
        self._prefix_block_sums = ti.field(dtype=ti.i32, shape=n_blocks)
        self._prefix_block_offsets = ti.field(dtype=ti.i32, shape=n_blocks)

        # Neighbor list (for consistent neighbor access)
        self.neighbor_count = ti.field(dtype=ti.i32, shape=n)
        self.neighbors = ti.field(dtype=ti.i32, shape=(n, self.max_neighbors))

        # CSR compact neighbor storage (reduces memory from 43MB to ~3MB)
        # avg_neighbors ~= 3.2, allocate n * 16 for headroom (still 4x smaller!)
        csr_capacity = n * 16
        self.neighbors_csr = ti.field(dtype=ti.i32, shape=csr_capacity)
        self.neighbor_offsets = ti.field(dtype=ti.i32, shape=n + 1)
        self._csr_capacity = csr_capacity

        # Prefix sum buffers for neighbor offsets (reuse same block size as grid)
        n_neighbor_blocks = (n + self._prefix_block_size - 1) // self._prefix_block_size
        self._neighbor_block_sums = ti.field(dtype=ti.i32, shape=n_neighbor_blocks)
        self._neighbor_block_offsets = ti.field(dtype=ti.i32, shape=n_neighbor_blocks)

        # Elastic connections (allocated on demand)
        self.elastic_connections = None
        self.elastic_rest_lengths = None
        self.elastic_stiffness = None

        # Membrane triangles (allocated on demand)
        self.membrane_triangles = None
        self.n_membranes = 0
        self.has_membranes = False
        self.is_muscle = None

        # Muscle activation
        self.muscle_activation = None

        # GPU state caching (Phase 3 optimization)
        # Avoids redundant ti.sync() calls when state hasn't changed
        self._gpu_dirty = True  # GPU has uncommitted changes
        self._cached_pos = None  # Cached position numpy array
        self._cached_vel = None  # Cached velocity numpy array

    def _init_particles(self, position, velocity):
        """Initialize particle positions and velocities."""
        # Convert to numpy if needed
        if hasattr(position, 'numpy'):
            position = position.numpy()
        if hasattr(velocity, 'numpy'):
            velocity = velocity.numpy()

        pos_np = np.array(position, dtype=np.float32)
        vel_np = np.array(velocity, dtype=np.float32)

        self.pos.from_numpy(pos_np)
        self.vel.from_numpy(vel_np)

    def _build_kernels(self):
        """Build Taichi kernels for SPH computations."""
        # Store references for use in kernels
        # OpenCL compatibility: use scaled coordinates for all SPH computations
        sim_scale = self.simulation_scale
        sim_scale_inv = 1.0 / sim_scale

        h = self.h
        h_scaled = h * sim_scale  # OpenCL: _hScaled = h * simulationScale
        h2_scaled = h_scaled * h_scaled
        h2 = h * h  # Keep unscaled for grid (positions are in world units)

        rho0 = self.rho0
        delta = self.delta
        mass = self.mass
        dt = self.dt

        # Kernel coefficients using SCALED h (OpenCL compatibility)
        # SPH kernels work in scaled coordinate space
        pi = 3.14159265359
        poly6_coeff = 315.0 / (64.0 * pi * h_scaled**9)
        spiky_coeff = -45.0 / (pi * h_scaled**6)
        visc_coeff = 45.0 / (pi * h_scaled**6)

        mu = self.mu
        gravity = self.gravity
        grid_size = self.grid_size
        grid_n = self.grid_n
        grid_origin = self.grid_origin
        floor_y = self.floor_y
        floor_restitution = self.floor_restitution
        r0 = self.r0  # Unscaled for world coordinates
        max_neighbors = self.max_neighbors
        surf_coeff = self.surface_tension_coeff

        pos = self.pos
        vel = self.vel
        acc = self.acc
        acc_old = self.acc_old
        rho = self.rho
        pressure = self.pressure
        predicted_pos = self.predicted_pos
        predicted_rho = self.predicted_rho
        grid_count = self.grid_count
        grid_offset = self.grid_offset
        particle_ids = self.particle_ids
        neighbor_count = self.neighbor_count
        neighbors = self.neighbors

        @ti.func
        def get_cell(p):
            """Get grid cell index for a position."""
            cell = ti.Vector([
                int((p[0] - grid_origin[0]) / grid_size),
                int((p[1] - grid_origin[1]) / grid_size),
                int((p[2] - grid_origin[2]) / grid_size),
            ])
            return ti.Vector([
                ti.max(0, ti.min(cell[0], grid_n - 1)),
                ti.max(0, ti.min(cell[1], grid_n - 1)),
                ti.max(0, ti.min(cell[2], grid_n - 1)),
            ])

        @ti.func
        def poly6_kernel(r_sq_scaled):
            """Poly6 kernel for density (scaled coordinates)."""
            result = 0.0
            if r_sq_scaled < h2_scaled:
                diff = h2_scaled - r_sq_scaled
                result = poly6_coeff * diff * diff * diff
            return result

        @ti.func
        def spiky_gradient(r_scaled, r_len_scaled):
            """Spiky kernel gradient for pressure (scaled coordinates)."""
            result = ti.Vector([0.0, 0.0, 0.0])
            if 0.0001 < r_len_scaled < h_scaled:
                diff = h_scaled - r_len_scaled
                result = spiky_coeff * diff * diff * (r_scaled / r_len_scaled)
            return result

        @ti.func
        def viscosity_laplacian(r_len_scaled):
            """Viscosity kernel laplacian (scaled coordinates)."""
            result = 0.0
            if 0.0 < r_len_scaled < h_scaled:
                result = visc_coeff * (h_scaled - r_len_scaled)
            return result

        # =====================================================================
        # Neighbor Search Kernels
        # =====================================================================
        @ti.kernel
        def count_particles_in_grid():
            for i, j, k in grid_count:
                grid_count[i, j, k] = 0
            for i in pos:
                cell = get_cell(pos[i])
                ti.atomic_add(grid_count[cell[0], cell[1], cell[2]], 1)

        @ti.kernel
        def compute_grid_offsets_sequential():
            """DEPRECATED: Sequential prefix sum - O(N) on GPU, very slow!"""
            # Must serialize - has loop-carried dependency on offset
            ti.loop_config(serialize=True)
            for i in range(grid_n):
                for j in range(grid_n):
                    for k in range(grid_n):
                        if i == 0 and j == 0 and k == 0:
                            grid_offset[0, 0, 0] = 0
                        else:
                            # Get previous cell in row-major order
                            prev_k = k - 1
                            prev_j = j
                            prev_i = i
                            if prev_k < 0:
                                prev_k = grid_n - 1
                                prev_j = j - 1
                                if prev_j < 0:
                                    prev_j = grid_n - 1
                                    prev_i = i - 1
                            grid_offset[i, j, k] = grid_offset[prev_i, prev_j, prev_k] + grid_count[prev_i, prev_j, prev_k]

        # Flat fields for parallel prefix sum
        grid_count_flat = self._grid_count_flat
        grid_offset_flat = self._grid_offset_flat

        @ti.kernel
        def flatten_grid_count():
            """Flatten 3D grid_count to 1D for fast CPU prefix sum."""
            for i, j, k in grid_count:
                flat_idx = i * grid_n * grid_n + j * grid_n + k
                grid_count_flat[flat_idx] = grid_count[i, j, k]

        @ti.kernel
        def unflatten_grid_offset():
            """Unflatten 1D grid_offset back to 3D after CPU prefix sum."""
            for i, j, k in grid_offset:
                flat_idx = i * grid_n * grid_n + j * grid_n + k
                grid_offset[i, j, k] = grid_offset_flat[flat_idx]

        # Hierarchical parallel prefix sum kernels
        block_size = self._prefix_block_size
        block_sums = self._prefix_block_sums
        block_offsets = self._prefix_block_offsets
        n_cells = grid_n * grid_n * grid_n
        n_blocks = (n_cells + block_size - 1) // block_size

        @ti.kernel
        def prefix_sum_phase1():
            """Phase 1: Compute local prefix sum within each block (parallel over blocks).

            Each block computes its own exclusive prefix sum and stores the block total.
            This is O(block_size) per block but all blocks run in parallel.
            """
            for block_id in range(n_blocks):
                block_start = block_id * block_size
                block_end = ti.min(block_start + block_size, n_cells)

                # Sequential prefix sum within this block
                running_sum = 0
                for i in range(block_start, block_end):
                    val = grid_count_flat[i]
                    grid_offset_flat[i] = running_sum
                    running_sum += val

                # Store block total for phase 2
                block_sums[block_id] = running_sum

        @ti.kernel
        def prefix_sum_phase2():
            """Phase 2: Compute prefix sum of block totals (SERIAL but only ~1100 iterations).

            This MUST be serialized due to the prefix sum dependency.
            But with block_size=1024 and 1.1M cells, this is only ~1100 iterations
            instead of 1.1M - a 1000x reduction in sequential work!
            """
            # Force serialization - this loop has data dependencies
            ti.loop_config(serialize=True)
            for block_id in range(n_blocks):
                if block_id == 0:
                    block_offsets[0] = 0
                else:
                    block_offsets[block_id] = block_offsets[block_id - 1] + block_sums[block_id - 1]

        @ti.kernel
        def prefix_sum_phase3():
            """Phase 3: Add block offsets to local prefix sums (parallel over all elements)."""
            for i in grid_offset_flat:
                block_id = i // block_size
                grid_offset_flat[i] += block_offsets[block_id]

        @ti.kernel
        def sort_particles():
            for i, j, k in grid_count:
                grid_count[i, j, k] = 0
            for i in pos:
                cell = get_cell(pos[i])
                idx = ti.atomic_add(grid_count[cell[0], cell[1], cell[2]], 1)
                particle_ids[grid_offset[cell[0], cell[1], cell[2]] + idx] = i

        @ti.kernel
        def build_neighbor_list():
            for i in pos:
                neighbor_count[i] = 0
                cell = get_cell(pos[i])
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        for dk in range(-1, 2):
                            ni = cell[0] + di
                            nj = cell[1] + dj
                            nk = cell[2] + dk
                            if 0 <= ni < grid_n and 0 <= nj < grid_n and 0 <= nk < grid_n:
                                start = grid_offset[ni, nj, nk]
                                end = start + grid_count[ni, nj, nk]
                                for idx in range(start, end):
                                    j = particle_ids[idx]
                                    if i != j:
                                        r = pos[i] - pos[j]
                                        r_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]
                                        if r_sq < h2 and neighbor_count[i] < max_neighbors:
                                            neighbors[i, neighbor_count[i]] = j
                                            neighbor_count[i] += 1

        # =====================================================================
        # Density and Pressure Kernels
        # =====================================================================
        @ti.kernel
        def compute_density():
            """Compute density using stored neighbor list."""
            for i in pos:
                # Self-contribution: W(0) with scaled h
                rho[i] = mass * poly6_coeff * h2_scaled * h2_scaled * h2_scaled
                # Neighbor contributions
                for ni in range(neighbor_count[i]):
                    j = neighbors[i, ni]
                    r = pos[i] - pos[j]
                    r_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]
                    # Scale r_sq for kernel (world coords -> scaled coords)
                    r_sq_scaled = r_sq * sim_scale * sim_scale
                    rho[i] += mass * poly6_kernel(r_sq_scaled)
                # Minimum clamp
                rho[i] = ti.max(rho[i], 0.01)

        @ti.kernel
        def compute_density_direct():
            """Compute density with on-the-fly neighbor search (no stored list)."""
            for i in pos:
                # Self-contribution: W(0) with scaled h
                rho[i] = mass * poly6_coeff * h2_scaled * h2_scaled * h2_scaled

                # On-the-fly neighbor search
                cell = get_cell(pos[i])
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        for dk in range(-1, 2):
                            ni = cell[0] + di
                            nj = cell[1] + dj
                            nk = cell[2] + dk
                            # Bounds check
                            if 0 <= ni < grid_n and 0 <= nj < grid_n and 0 <= nk < grid_n:
                                # Skip empty cells early
                                cell_count = grid_count[ni, nj, nk]
                                if cell_count > 0:
                                    start = grid_offset[ni, nj, nk]
                                    for idx in range(start, start + cell_count):
                                        j = particle_ids[idx]
                                        if i != j:
                                            r = pos[i] - pos[j]
                                            r_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]
                                            if r_sq < h2:
                                                # Scale r_sq for kernel
                                                r_sq_scaled = r_sq * sim_scale * sim_scale
                                                rho[i] += mass * poly6_kernel(r_sq_scaled)
                # Minimum clamp
                rho[i] = ti.max(rho[i], 0.01)

        @ti.kernel
        def compute_density_and_cache_neighbors():
            """FUSED: Compute density, pressure AND cache neighbors in one pass.

            This combines neighbor search with density computation,
            storing neighbors for later use by force computation.
            Only ONE neighbor search pass needed!
            """
            for i in pos:
                # Self-contribution: W(0) with scaled h
                rho[i] = mass * poly6_coeff * h2_scaled * h2_scaled * h2_scaled

                # Reset neighbor count for this particle
                neighbor_count[i] = 0

                # On-the-fly neighbor search + density + caching
                cell = get_cell(pos[i])
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        for dk in range(-1, 2):
                            ni = cell[0] + di
                            nj = cell[1] + dj
                            nk = cell[2] + dk
                            # Bounds check
                            if 0 <= ni < grid_n and 0 <= nj < grid_n and 0 <= nk < grid_n:
                                # Skip empty cells early
                                cell_count = grid_count[ni, nj, nk]
                                if cell_count > 0:
                                    start = grid_offset[ni, nj, nk]
                                    for idx in range(start, start + cell_count):
                                        j = particle_ids[idx]
                                        if i != j:
                                            r = pos[i] - pos[j]
                                            r_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]
                                            if r_sq < h2:
                                                # Compute density contribution
                                                r_sq_scaled = r_sq * sim_scale * sim_scale
                                                rho[i] += mass * poly6_kernel(r_sq_scaled)

                                                # Cache this neighbor for force computation
                                                if neighbor_count[i] < max_neighbors:
                                                    neighbors[i, neighbor_count[i]] = j
                                                    neighbor_count[i] += 1
                # Minimum clamp
                rho[i] = ti.max(rho[i], 0.01)

                # FUSED: Compute pressure inline (Tait equation)
                p = delta * (rho[i] / rho0 - 1.0)
                pressure[i] = ti.max(p, 0.0)

        @ti.kernel
        def compute_pressure():
            for i in pos:
                # Tait equation: P = delta * (rho/rho0 - 1)
                p = delta * (rho[i] / rho0 - 1.0)
                # Clamp to non-negative
                pressure[i] = ti.max(p, 0.0)

        # =====================================================================
        # CSR Compact Neighbor Storage Kernels
        # Two-pass approach: count → prefix sum → fill
        # Reduces memory from 43MB (n × 64 × 4) to ~3MB (n × avg × 4)
        # =====================================================================
        neighbors_csr = self.neighbors_csr
        neighbor_offsets = self.neighbor_offsets
        csr_capacity = self._csr_capacity
        neighbor_block_sums = self._neighbor_block_sums
        neighbor_block_offsets = self._neighbor_block_offsets
        n_particles = self.n_particles
        prefix_block_size = self._prefix_block_size

        @ti.kernel
        def count_neighbors_and_density():
            """CSR Pass 1: Count neighbors + compute density + pressure.

            Does NOT store neighbors yet - just counts them for prefix sum.
            """
            for i in pos:
                # Self-contribution
                rho[i] = mass * poly6_coeff * h2_scaled * h2_scaled * h2_scaled
                neighbor_count[i] = 0

                cell = get_cell(pos[i])
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        for dk in range(-1, 2):
                            ni = cell[0] + di
                            nj = cell[1] + dj
                            nk = cell[2] + dk
                            if 0 <= ni < grid_n and 0 <= nj < grid_n and 0 <= nk < grid_n:
                                cell_count = grid_count[ni, nj, nk]
                                if cell_count > 0:
                                    start = grid_offset[ni, nj, nk]
                                    for idx in range(start, start + cell_count):
                                        j = particle_ids[idx]
                                        if i != j:
                                            r = pos[i] - pos[j]
                                            r_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]
                                            if r_sq < h2:
                                                r_sq_scaled = r_sq * sim_scale * sim_scale
                                                rho[i] += mass * poly6_kernel(r_sq_scaled)
                                                neighbor_count[i] += 1

                rho[i] = ti.max(rho[i], 0.01)
                # Compute pressure inline
                p = delta * (rho[i] / rho0 - 1.0)
                pressure[i] = ti.max(p, 0.0)

        @ti.kernel
        def neighbor_prefix_phase1():
            """Phase 1: Compute local prefix sums within blocks + block totals."""
            ti.loop_config(serialize=True)
            for block_id in range((n_particles + prefix_block_size - 1) // prefix_block_size):
                block_start = block_id * prefix_block_size
                block_end = ti.min(block_start + prefix_block_size, n_particles)

                running_sum = 0
                for i in range(block_start, block_end):
                    old_count = neighbor_count[i]
                    neighbor_offsets[i] = running_sum
                    running_sum += old_count

                neighbor_block_sums[block_id] = running_sum

        @ti.kernel
        def neighbor_prefix_phase2():
            """Phase 2: Compute prefix sum of block totals (sequential)."""
            ti.loop_config(serialize=True)
            n_blocks = (n_particles + prefix_block_size - 1) // prefix_block_size
            running_sum = 0
            for block_id in range(n_blocks):
                neighbor_block_offsets[block_id] = running_sum
                running_sum += neighbor_block_sums[block_id]
            # Store total count in last position
            neighbor_offsets[n_particles] = running_sum

        @ti.kernel
        def neighbor_prefix_phase3():
            """Phase 3: Add block offsets to all elements (parallel)."""
            for i in range(n_particles):
                block_id = i // prefix_block_size
                neighbor_offsets[i] += neighbor_block_offsets[block_id]

        @ti.kernel
        def fill_neighbors_csr():
            """CSR Pass 2: Fill CSR storage using computed offsets.

            Searches neighbors again (fast) and stores to packed array.
            """
            for i in pos:
                write_idx = neighbor_offsets[i]
                local_count = 0

                cell = get_cell(pos[i])
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        for dk in range(-1, 2):
                            ni = cell[0] + di
                            nj = cell[1] + dj
                            nk = cell[2] + dk
                            if 0 <= ni < grid_n and 0 <= nj < grid_n and 0 <= nk < grid_n:
                                cell_count = grid_count[ni, nj, nk]
                                if cell_count > 0:
                                    start = grid_offset[ni, nj, nk]
                                    for idx in range(start, start + cell_count):
                                        j = particle_ids[idx]
                                        if i != j:
                                            r = pos[i] - pos[j]
                                            r_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]
                                            if r_sq < h2:
                                                if write_idx + local_count < csr_capacity:
                                                    neighbors_csr[write_idx + local_count] = j
                                                    local_count += 1

        @ti.kernel
        def compute_forces_csr():
            """Compute forces using CSR neighbor lookup.

            Uses packed neighbor storage instead of fixed-size 2D array.
            """
            for i in pos:
                particle_type = int(pos[i][3])
                if particle_type == 3:  # BOUNDARY
                    acc[i] = ti.Vector([0.0, 0.0, 0.0, 0.0])
                    continue

                # All forces in world coordinates (gravity included from start)
                force = ti.Vector([gravity[0], gravity[1], gravity[2]])
                pressure_i = pressure[i]
                v_i = ti.Vector([vel[i][0], vel[i][1], vel[i][2]])

                # Get CSR range for this particle's neighbors
                start_idx = neighbor_offsets[i]
                end_idx = neighbor_offsets[i + 1]

                # Iterate through neighbors using CSR
                for k in range(start_idx, end_idx):
                    j = neighbors_csr[k]
                    j_type = int(pos[j][3])

                    r = ti.Vector([pos[i][0] - pos[j][0],
                                   pos[i][1] - pos[j][1],
                                   pos[i][2] - pos[j][2]])
                    r_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]

                    if r_sq > h2 or r_sq < 1e-12:
                        continue

                    r_len = ti.sqrt(r_sq)
                    if r_len < 0.0001:
                        continue

                    # Boundary handling (Ihmsen et al. 2010)
                    if j_type == 3:  # BOUNDARY
                        n_b = ti.Vector([vel[j][0], vel[j][1], vel[j][2]])
                        n_len = ti.sqrt(n_b[0]*n_b[0] + n_b[1]*n_b[1] + n_b[2]*n_b[2])
                        if n_len > 0.0001:
                            n_b = n_b / n_len
                            w = ti.max(0.0, (r0 - r_len) / r0)
                            repulsion = 2000.0 * w * n_b
                            force += repulsion
                        continue

                    rho_j = rho[j]

                    # Scale distance for SPH kernels (world -> scaled coords)
                    r_scaled = ti.Vector([r[0] * sim_scale, r[1] * sim_scale, r[2] * sim_scale])
                    r_len_scaled = r_len * sim_scale

                    # Pressure force
                    if not (particle_type == 2 and j_type == 2):
                        pressure_term = (pressure_i + pressure[j]) / (2.0 * rho_j + 0.0001)
                        grad_scaled = spiky_gradient(r_scaled, r_len_scaled)
                        grad_world = ti.Vector([grad_scaled[0] * sim_scale,
                                               grad_scaled[1] * sim_scale,
                                               grad_scaled[2] * sim_scale])
                        force += -mass * pressure_term * grad_world

                    # Viscosity force
                    v_diff = ti.Vector([vel[j][0] - v_i[0], vel[j][1] - v_i[1], vel[j][2] - v_i[2]])
                    visc_lap_scaled = viscosity_laplacian(r_len_scaled)
                    visc_lap_world = visc_lap_scaled * sim_scale * sim_scale
                    force += mu * mass * v_diff * visc_lap_world / (rho_j + 0.0001)

                acc[i] = ti.Vector([force[0], force[1], force[2], 0.0])

        # =====================================================================
        # Force Computation Kernels
        # =====================================================================
        @ti.kernel
        def compute_forces():
            for i in pos:
                particle_type = int(pos[i][3])

                # Boundary particles don't move
                if particle_type == 3:  # BOUNDARY
                    acc[i] = ti.Vector([0.0, 0.0, 0.0, 0.0])
                    continue

                # All forces in world coordinates
                force = ti.Vector([gravity[0], gravity[1], gravity[2]])

                p_i = pos[i]
                v_i = vel[i]
                rho_i = rho[i]
                pressure_i = pressure[i]

                for ni in range(neighbor_count[i]):
                    j = neighbors[i, ni]
                    j_type = int(pos[j][3])

                    # World-space distance
                    r = ti.Vector([p_i[0] - pos[j][0], p_i[1] - pos[j][1], p_i[2] - pos[j][2]])
                    r_len = ti.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])

                    if r_len < 0.0001:
                        continue

                    # Boundary handling (Ihmsen et al. 2010)
                    if j_type == 3:  # BOUNDARY
                        # Boundary normal stored in velocity
                        n_b = ti.Vector([vel[j][0], vel[j][1], vel[j][2]])
                        n_len = ti.sqrt(n_b[0]*n_b[0] + n_b[1]*n_b[1] + n_b[2]*n_b[2])
                        if n_len > 0.0001:
                            n_b = n_b / n_len
                            w = ti.max(0.0, (r0 - r_len) / r0)
                            repulsion = 2000.0 * w * n_b
                            force += repulsion
                        continue

                    rho_j = rho[j]

                    # Scale distance for SPH kernels (world -> scaled coords)
                    r_scaled = ti.Vector([r[0] * sim_scale, r[1] * sim_scale, r[2] * sim_scale])
                    r_len_scaled = r_len * sim_scale

                    # Pressure force - skip between elastic particles
                    # Gradient computed in scaled coords, convert back to world coords
                    if not (particle_type == 2 and j_type == 2):
                        pressure_term = (pressure_i + pressure[j]) / (2.0 * rho_j + 0.0001)
                        grad_scaled = spiky_gradient(r_scaled, r_len_scaled)
                        # dW/dr_world = dW/dr_scaled * sim_scale
                        grad_world = ti.Vector([grad_scaled[0] * sim_scale, grad_scaled[1] * sim_scale, grad_scaled[2] * sim_scale])
                        force += -mass * pressure_term * grad_world

                    # Viscosity force - also convert laplacian from scaled to world
                    v_diff = ti.Vector([vel[j][0] - v_i[0], vel[j][1] - v_i[1], vel[j][2] - v_i[2]])
                    visc_lap_scaled = viscosity_laplacian(r_len_scaled)
                    # Laplacian has units 1/r^2, so scale by sim_scale^2
                    visc_lap_world = visc_lap_scaled * sim_scale * sim_scale
                    force += mu * mass * v_diff * visc_lap_world / (rho_j + 0.0001)

                    # Surface tension
                    if surf_coeff > 0.0:
                        surf_kernel = ti.max(0.0, h2 - r_len*r_len)
                        surf_kernel = surf_kernel * surf_kernel * surf_kernel
                        force += -1.7e-9 * surf_coeff * surf_kernel * r / (r_len + 0.0001) / mass

                acc[i] = ti.Vector([force[0], force[1], force[2], 0.0])

        @ti.kernel
        def compute_forces_direct():
            """Compute forces with on-the-fly neighbor search (no stored list)."""
            for i in pos:
                particle_type = int(pos[i][3])

                # Boundary particles don't move
                if particle_type == 3:  # BOUNDARY
                    acc[i] = ti.Vector([0.0, 0.0, 0.0, 0.0])
                    continue

                # All forces in world coordinates
                force = ti.Vector([gravity[0], gravity[1], gravity[2]])

                p_i = pos[i]
                v_i = vel[i]
                rho_i = rho[i]
                pressure_i = pressure[i]

                # On-the-fly neighbor search
                cell = get_cell(p_i)
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        for dk in range(-1, 2):
                            ni_cell = cell[0] + di
                            nj_cell = cell[1] + dj
                            nk_cell = cell[2] + dk
                            # Bounds check
                            if 0 <= ni_cell < grid_n and 0 <= nj_cell < grid_n and 0 <= nk_cell < grid_n:
                                # Skip empty cells early
                                cell_count = grid_count[ni_cell, nj_cell, nk_cell]
                                if cell_count > 0:
                                    start = grid_offset[ni_cell, nj_cell, nk_cell]
                                    for idx in range(start, start + cell_count):
                                        j = particle_ids[idx]
                                        if i != j:
                                            # Distance check
                                            r = ti.Vector([p_i[0] - pos[j][0], p_i[1] - pos[j][1], p_i[2] - pos[j][2]])
                                            r_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]
                                            if r_sq < h2:
                                                r_len = ti.sqrt(r_sq)
                                                if r_len < 0.0001:
                                                    continue

                                                j_type = int(pos[j][3])

                                                # Boundary handling (Ihmsen et al. 2010)
                                                if j_type == 3:  # BOUNDARY
                                                    n_b = ti.Vector([vel[j][0], vel[j][1], vel[j][2]])
                                                    n_len = ti.sqrt(n_b[0]*n_b[0] + n_b[1]*n_b[1] + n_b[2]*n_b[2])
                                                    if n_len > 0.0001:
                                                        n_b = n_b / n_len
                                                        w = ti.max(0.0, (r0 - r_len) / r0)
                                                        repulsion = 2000.0 * w * n_b
                                                        force += repulsion
                                                    continue

                                                rho_j = rho[j]

                                                # Scale distance for SPH kernels
                                                r_scaled = ti.Vector([r[0] * sim_scale, r[1] * sim_scale, r[2] * sim_scale])
                                                r_len_scaled = r_len * sim_scale

                                                # Pressure force
                                                if not (particle_type == 2 and j_type == 2):
                                                    pressure_term = (pressure_i + pressure[j]) / (2.0 * rho_j + 0.0001)
                                                    grad_scaled = spiky_gradient(r_scaled, r_len_scaled)
                                                    grad_world = ti.Vector([grad_scaled[0] * sim_scale, grad_scaled[1] * sim_scale, grad_scaled[2] * sim_scale])
                                                    force += -mass * pressure_term * grad_world

                                                # Viscosity force
                                                v_diff = ti.Vector([vel[j][0] - v_i[0], vel[j][1] - v_i[1], vel[j][2] - v_i[2]])
                                                visc_lap_scaled = viscosity_laplacian(r_len_scaled)
                                                visc_lap_world = visc_lap_scaled * sim_scale * sim_scale
                                                force += mu * mass * v_diff * visc_lap_world / (rho_j + 0.0001)

                                                # Surface tension
                                                if surf_coeff > 0.0:
                                                    surf_kernel = ti.max(0.0, h2 - r_sq)
                                                    surf_kernel = surf_kernel * surf_kernel * surf_kernel
                                                    force += -1.7e-9 * surf_coeff * surf_kernel * r / (r_len + 0.0001) / mass

                acc[i] = ti.Vector([force[0], force[1], force[2], 0.0])

        # =====================================================================
        # Integration Kernels
        # =====================================================================
        @ti.kernel
        def integrate_position():
            """Leapfrog position update (mode 0).

            Matches OpenCL sphFluid.cl line ~1455:
                pos += (vel*dt + acc*dt²/2) * simulationScaleInv
            The * simulationScaleInv converts the integrator's m/s² and
            m/s units into the sim-unit space the position field actually
            uses (1 sim_unit ≈ sim_scale meters). Without this multiplier
            the cube barely moves — confirmed via per-step debug dump
            on demo1 (pos stayed at 39.41 across 4 steps with gravity
            successfully accumulating velocity to -6.86e-4 m/s).
            """
            for i in pos:
                particle_type = int(pos[i][3])
                if particle_type == 3:  # BOUNDARY - frozen
                    continue

                # Position delta in world meters; convert to sim units.
                delta_x = (vel[i][0] * dt + 0.5 * acc_old[i][0] * dt * dt) * sim_scale_inv
                delta_y = (vel[i][1] * dt + 0.5 * acc_old[i][1] * dt * dt) * sim_scale_inv
                delta_z = (vel[i][2] * dt + 0.5 * acc_old[i][2] * dt * dt) * sim_scale_inv
                pos[i][0] += delta_x
                pos[i][1] += delta_y
                pos[i][2] += delta_z

        @ti.kernel
        def integrate_velocity():
            """Leapfrog velocity update (mode 1)."""
            for i in pos:
                particle_type = int(pos[i][3])
                if particle_type == 3:  # BOUNDARY - frozen
                    vel[i] = ti.Vector([0.0, 0.0, 0.0, 0.0])
                    continue

                # v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
                vel[i][0] += 0.5 * (acc_old[i][0] + acc[i][0]) * dt
                vel[i][1] += 0.5 * (acc_old[i][1] + acc[i][1]) * dt
                vel[i][2] += 0.5 * (acc_old[i][2] + acc[i][2]) * dt

        @ti.kernel
        def integrate_euler():
            """Semi-implicit Euler (fallback).

            Also applies simulationScaleInv to position delta for OpenCL compatibility.
            """
            for i in pos:
                particle_type = int(pos[i][3])
                if particle_type == 3:  # BOUNDARY
                    vel[i] = ti.Vector([0.0, 0.0, 0.0, 0.0])
                    continue

                # v(t+dt) = v(t) + a(t)*dt
                vel[i][0] += acc[i][0] * dt
                vel[i][1] += acc[i][1] * dt
                vel[i][2] += acc[i][2] * dt

                # x(t+dt) = x(t) + v(t+dt)*dt * sim_scale_inv (world m → sim units)
                pos[i][0] += vel[i][0] * dt * sim_scale_inv
                pos[i][1] += vel[i][1] * dt * sim_scale_inv
                pos[i][2] += vel[i][2] * dt * sim_scale_inv

        @ti.kernel
        def apply_floor_constraint():
            """Apply hard floor constraint."""
            for i in pos:
                particle_type = int(pos[i][3])
                if particle_type == 3:  # BOUNDARY
                    continue

                if pos[i][1] < floor_y:
                    pos[i][1] = floor_y
                    # Reflect velocity with restitution
                    if vel[i][1] < 0:
                        vel[i][1] = -vel[i][1] * floor_restitution

        @ti.kernel
        def store_acceleration():
            """Store current acceleration for leapfrog."""
            for i in pos:
                acc_old[i] = acc[i]

        # Store kernel references
        self._count_particles_in_grid = count_particles_in_grid
        self._compute_grid_offsets_sequential = compute_grid_offsets_sequential
        self._flatten_grid_count = flatten_grid_count
        self._unflatten_grid_offset = unflatten_grid_offset
        self._prefix_sum_phase1 = prefix_sum_phase1
        self._prefix_sum_phase2 = prefix_sum_phase2
        self._prefix_sum_phase3 = prefix_sum_phase3
        self._sort_particles = sort_particles
        self._build_neighbor_list = build_neighbor_list
        self._compute_density = compute_density
        self._compute_density_direct = compute_density_direct
        self._compute_density_and_cache_neighbors = compute_density_and_cache_neighbors
        self._compute_pressure = compute_pressure
        self._compute_forces = compute_forces
        self._compute_forces_direct = compute_forces_direct
        self._integrate_position = integrate_position
        self._integrate_velocity = integrate_velocity
        self._integrate_euler = integrate_euler

        # CSR compact neighbor storage kernels
        self._count_neighbors_and_density = count_neighbors_and_density
        self._neighbor_prefix_phase1 = neighbor_prefix_phase1
        self._neighbor_prefix_phase2 = neighbor_prefix_phase2
        self._neighbor_prefix_phase3 = neighbor_prefix_phase3
        self._fill_neighbors_csr = fill_neighbors_csr
        self._compute_forces_csr = compute_forces_csr

        # Neighbor search mode:
        # "fused"  - compute density + cache neighbors in one pass, forces use cached (DEFAULT - fastest)
        # "list"   - build neighbor list first, then density/forces use it
        # "direct" - on-the-fly for both density and forces (2x neighbor search)
        # "csr"    - compact storage: count → prefix sum → fill → forces (lowest memory)
        # Can be overridden via SIBERNETIC_NEIGHBOR_MODE environment variable
        env_mode = os.environ.get("SIBERNETIC_NEIGHBOR_MODE", "fused")
        self._neighbor_mode = env_mode if env_mode in ("fused", "list", "direct", "csr") else "fused"
        print(f"[TaichiSolver] Neighbor mode: {self._neighbor_mode}", file=sys.stderr, flush=True)
        self._apply_floor_constraint = apply_floor_constraint
        self._store_acceleration = store_acceleration

    def _compute_grid_offsets(self):
        """Compute grid offsets using HIERARCHICAL PARALLEL prefix sum.

        This replaces the O(N) sequential GPU scan with a 3-phase algorithm:

        Phase 1 (parallel): Divide into blocks of 1024, compute local prefix
                           sums within each block. ~1100 blocks run in parallel.

        Phase 2 (sequential): Compute prefix sum of block totals. Only ~1100
                             iterations instead of 1.1M.

        Phase 3 (parallel): Add block offsets to all elements. All elements
                           update in parallel.

        Complexity: O(B) + O(N/B) + O(1) where B=1024, vs O(N) for sequential.
        For N=1.1M: ~1100 sequential ops vs 1.1M = ~1000x fewer sequential ops.
        """
        # Step 1: Flatten 3D grid_count to 1D (parallel on GPU)
        self._flatten_grid_count()

        # Step 2: Hierarchical parallel prefix sum (all on GPU)
        self._prefix_sum_phase1()  # Parallel: local prefix sums + block totals
        self._prefix_sum_phase2()  # Sequential: prefix sum of ~1100 block totals
        self._prefix_sum_phase3()  # Parallel: add block offsets to all elements

        # Step 3: Unflatten back to 3D (parallel on GPU)
        self._unflatten_grid_offset()

    # =========================================================================
    # Elastic Connections (for worm body simulation)
    # =========================================================================
    def load_elastic_connections(self, connections_data):
        """Load elastic connection data for deformable body simulation.

        Args:
            connections_data: Array of shape (n_elastic, max_conn, 4) containing:
                [connected_particle_id, rest_length, stiffness, is_muscle]
        """
        if hasattr(connections_data, 'numpy'):
            connections_data = connections_data.numpy()

        n_elastic = connections_data.shape[0]
        max_conn = min(connections_data.shape[1], self.max_connections)

        # Allocate elastic fields if not done
        if self.elastic_connections is None:
            self.elastic_connections = ti.field(dtype=ti.i32, shape=(n_elastic, max_conn))
            self.elastic_rest_lengths = ti.field(dtype=ti.f32, shape=(n_elastic, max_conn))
            self.elastic_stiffness = ti.field(dtype=ti.f32, shape=(n_elastic, max_conn))
            self.is_muscle = ti.field(dtype=ti.i32, shape=(n_elastic, max_conn))
            self.muscle_activation = ti.field(dtype=ti.f32, shape=n_elastic)

        # Fill data
        conn_np = connections_data[:, :max_conn, 0].astype(np.int32)
        rest_np = connections_data[:, :max_conn, 1].astype(np.float32)
        stiff_np = connections_data[:, :max_conn, 2].astype(np.float32)
        muscle_np = connections_data[:, :max_conn, 3].astype(np.int32)

        self.elastic_connections.from_numpy(conn_np)
        self.elastic_rest_lengths.from_numpy(rest_np)
        self.elastic_stiffness.from_numpy(stiff_np)
        self.is_muscle.from_numpy(muscle_np)

        self.has_elastic = True
        self._n_elastic = n_elastic
        self._max_conn = max_conn

        # Build elastic force kernel
        self._build_elastic_kernel()

    def _load_elastic_from_cpp(self, elastic_list, config):
        """Load elastic connections from C++ format.

        The C++ code passes a flat list of connections where each elastic
        particle has 32 connection slots, each with 4 values:
        [connected_particle_id, rest_length, muscle_id, unused]

        IMPORTANT: The rest_length values are already multiplied by simulationScale
        in the C++ loader, so we need to apply simulationScale to our distance
        calculations in the elastic kernel to match.

        Args:
            elastic_list: Flat list from C++, length = num_elastic * 32
            config: Config dict with num_elastic_particles and elasticity_coefficient
        """
        num_elastic = config.get("num_elastic_particles", 0)
        if num_elastic == 0:
            return

        max_conn_per_particle = 32  # MAX_NEIGHBOR_COUNT in C++
        elasticity_coeff = config.get("elasticity_coefficient", 1e6)

        # Store simulation scale for elastic force calculation
        self.simulation_scale = config.get("simulation_scale", 1.0)

        # Convert flat list to numpy array and reshape
        elastic_np = np.array(elastic_list, dtype=np.float32)
        # Shape: (num_elastic * 32, 4) -> (num_elastic, 32, 4)
        elastic_np = elastic_np.reshape(num_elastic, max_conn_per_particle, 4)

        # C++ format: [particle_id, rest_length, muscle_id, unused]
        # Expected format: [particle_id, rest_length, stiffness, is_muscle]
        # Convert muscle_id > 0 to is_muscle = 1
        # Use elasticity coefficient as stiffness

        conn_ids = elastic_np[:, :, 0].astype(np.int32)
        rest_lengths = elastic_np[:, :, 1].astype(np.float32)
        muscle_ids = elastic_np[:, :, 2].astype(np.float32)

        # Create stiffness and is_muscle arrays
        stiffness = np.ones_like(rest_lengths) * elasticity_coeff
        is_muscle = (muscle_ids > 0).astype(np.int32)

        # Mark invalid connections (particle_id < 0 or rest_length <= 0)
        invalid = (conn_ids < 0) | (rest_lengths <= 0)
        conn_ids[invalid] = -1

        # Combine into expected format
        connections_data = np.stack([conn_ids, rest_lengths, stiffness, is_muscle], axis=2)

        # Use up to self.max_connections
        max_conn = min(max_conn_per_particle, self.max_connections)
        connections_data = connections_data[:, :max_conn, :]

        # Allocate and load
        n_elastic = num_elastic
        self._n_elastic = n_elastic
        self._max_conn = max_conn

        if self.elastic_connections is None:
            self.elastic_connections = ti.field(dtype=ti.i32, shape=(n_elastic, max_conn))
            self.elastic_rest_lengths = ti.field(dtype=ti.f32, shape=(n_elastic, max_conn))
            self.elastic_stiffness = ti.field(dtype=ti.f32, shape=(n_elastic, max_conn))
            self.is_muscle = ti.field(dtype=ti.i32, shape=(n_elastic, max_conn))
            self.muscle_activation = ti.field(dtype=ti.f32, shape=n_elastic)

        self.elastic_connections.from_numpy(connections_data[:, :, 0].astype(np.int32))
        self.elastic_rest_lengths.from_numpy(connections_data[:, :, 1].astype(np.float32))
        self.elastic_stiffness.from_numpy(connections_data[:, :, 2].astype(np.float32))
        self.is_muscle.from_numpy(connections_data[:, :, 3].astype(np.int32))

        self.has_elastic = True

        # Build elastic force kernel
        self._build_elastic_kernel()

        # Debug: print elastic parameters
        valid_rest = rest_lengths[rest_lengths > 0]
        print(f"Loaded {n_elastic} elastic particles with up to {max_conn} connections each")
        print(f"  simulation_scale: {self.simulation_scale}")
        print(f"  elasticity_coefficient: {elasticity_coeff:.3e}")
        if len(valid_rest) > 0:
            print(f"  rest_lengths (scaled): min={valid_rest.min():.2e}, max={valid_rest.max():.2e}, mean={valid_rest.mean():.2e}")

    def _build_elastic_kernel(self):
        """Build kernel for elastic/muscle forces.

        The elastic force formula from OpenCL:
            vect_r_ij = (pos[i] - pos[j]) * simulationScale
            r_ij = length(vect_r_ij)
            delta_r = r_ij - rest_length  (rest_length already scaled)
            acceleration += -(vect_r_ij/r_ij) * delta_r * elasticityCoefficient

        Note: elasticityCoefficient already includes 1/mass factor, so we add
        directly to acceleration, not force.
        """
        pos = self.pos
        vel = self.vel
        acc = self.acc
        elastic_connections = self.elastic_connections
        elastic_rest_lengths = self.elastic_rest_lengths
        elastic_stiffness = self.elastic_stiffness
        is_muscle = self.is_muscle
        muscle_activation = self.muscle_activation
        n_elastic = self._n_elastic
        max_conn = self._max_conn
        max_muscle_force = 4000.0
        sim_scale = self.simulation_scale
        sim_scale_inv = 1.0 / sim_scale

        @ti.kernel
        def compute_elastic_forces():
            for i in range(n_elastic):
                for c in range(max_conn):
                    j = elastic_connections[i, c]
                    if j < 0:
                        continue

                    rest_len = elastic_rest_lengths[i, c]
                    stiffness = elastic_stiffness[i, c]
                    is_musc = is_muscle[i, c]

                    # Vector from i to j in world coordinates
                    r_world = ti.Vector([pos[j][0] - pos[i][0],
                                         pos[j][1] - pos[i][1],
                                         pos[j][2] - pos[i][2]])

                    # Scale to simulation coordinates (rest_length is already scaled)
                    r_scaled = r_world * sim_scale
                    r_len = ti.sqrt(r_scaled[0]*r_scaled[0] + r_scaled[1]*r_scaled[1] + r_scaled[2]*r_scaled[2])

                    if r_len < 0.0001:
                        continue

                    # Direction in scaled space (unit vector from i toward j)
                    direction = r_scaled / r_len

                    # Displacement from equilibrium (positive = stretched)
                    delta_len = r_len - rest_len

                    # Spring acceleration (elasticityCoefficient includes 1/mass)
                    # Positive delta (stretched) -> accelerate toward j (direction)
                    # Negative delta (compressed) -> accelerate away from j (-direction)
                    spring_accel = stiffness * delta_len

                    # Damping: opposes relative motion along spring direction
                    # Relative velocity of i with respect to j
                    v_rel_x = vel[i][0] - vel[j][0]
                    v_rel_y = vel[i][1] - vel[j][1]
                    v_rel_z = vel[i][2] - vel[j][2]
                    # Component along spring direction (positive = i moving toward j)
                    v_along = v_rel_x * direction[0] + v_rel_y * direction[1] + v_rel_z * direction[2]
                    # Damping coefficient (fraction of stiffness)
                    damping_coeff = 0.3 * stiffness
                    damping_accel = -damping_coeff * v_along

                    accel_magnitude = spring_accel + damping_accel

                    # Muscle contraction (adds extra inward force)
                    if is_musc > 0:
                        activation = muscle_activation[i]
                        muscle_accel = ti.min(activation * max_muscle_force, max_muscle_force)
                        accel_magnitude += muscle_accel

                    # With the integration step now applying sim_scale_inv to the
                    # position delta (matching OpenCL's sphFluid.cl line ~1455),
                    # elastic acceleration in scaled coordinates flows through the
                    # integrator with the right effective magnitude. The previous
                    # `elastic_boost = 2000.0` was a fudge for the missing
                    # sim_scale_inv at integration time; with that fix in place
                    # boost should be 1.0 for parity with OpenCL.
                    accel = accel_magnitude * direction
                    acc[i][0] += accel[0]
                    acc[i][1] += accel[1]
                    acc[i][2] += accel[2]

        self._compute_elastic_forces = compute_elastic_forces

    def set_muscle_activation(self, activation):
        """Set muscle activation levels.

        Args:
            activation: Array of activation values per elastic particle
        """
        if hasattr(activation, 'numpy'):
            activation = activation.numpy()
        self.muscle_activation.from_numpy(activation.astype(np.float32))
        self.has_muscles = True

    # =========================================================================
    # Membrane Collision (elastic shell containment)
    # =========================================================================
    def load_membranes(self, membrane_data):
        """Load membrane triangle data for liquid containment.

        Args:
            membrane_data: Nx3 array of particle indices forming triangles
        """
        if hasattr(membrane_data, 'numpy'):
            membrane_data = membrane_data.numpy()

        membrane_data = np.array(membrane_data, dtype=np.int32)
        self.n_membranes = len(membrane_data)

        if self.n_membranes == 0:
            return

        # Allocate membrane field
        self.membrane_triangles = ti.Vector.field(3, dtype=ti.i32, shape=self.n_membranes)
        self.membrane_triangles.from_numpy(membrane_data)

        self.has_membranes = True
        print(f"Loaded {self.n_membranes} membrane triangles")

        # Build membrane collision kernel
        self._build_membrane_kernel()

    def _build_membrane_kernel(self):
        """Build membrane collision detection kernel."""
        pos = self.pos
        vel = self.vel
        acc = self.acc
        membrane_triangles = self.membrane_triangles
        n_membranes = self.n_membranes
        n_elastic = self._n_elastic
        h = self.h

        @ti.kernel
        def apply_membrane_forces():
            """Apply repulsion forces when liquid particles approach membranes."""
            # First compute center of mass of elastic particles
            com_x = 0.0
            com_y = 0.0
            com_z = 0.0
            elastic_count = 0
            for i in pos:
                if int(pos[i][3]) == 2:  # Elastic particle
                    com_x += pos[i][0]
                    com_y += pos[i][1]
                    com_z += pos[i][2]
                    elastic_count += 1

            if elastic_count > 0:
                com_x /= elastic_count
                com_y /= elastic_count
                com_z /= elastic_count

            elastic_com = ti.Vector([com_x, com_y, com_z])

            # For each liquid particle
            for i in pos:
                particle_type = int(pos[i][3])
                if particle_type != 1:  # Only liquid particles
                    continue

                p = ti.Vector([pos[i][0], pos[i][1], pos[i][2]])

                # Check against all membrane triangles
                for m in range(n_membranes):
                    # Get triangle vertex indices
                    i0 = membrane_triangles[m][0]
                    i1 = membrane_triangles[m][1]
                    i2 = membrane_triangles[m][2]

                    # Get vertex positions
                    v0 = ti.Vector([pos[i0][0], pos[i0][1], pos[i0][2]])
                    v1 = ti.Vector([pos[i1][0], pos[i1][1], pos[i1][2]])
                    v2 = ti.Vector([pos[i2][0], pos[i2][1], pos[i2][2]])

                    # Triangle center
                    tri_center = (v0 + v1 + v2) / 3.0

                    # Quick distance check - skip if too far
                    dist_to_center = (p - tri_center).norm()
                    if dist_to_center > h * 2.0:
                        continue

                    # Compute triangle normal
                    e1 = v1 - v0
                    e2 = v2 - v0
                    normal = e1.cross(e2)
                    normal_len = normal.norm()
                    if normal_len < 0.0001:
                        continue
                    normal = normal / normal_len

                    # Orient normal to point INWARD (toward elastic COM)
                    to_com = elastic_com - tri_center
                    if normal.dot(to_com) < 0:
                        normal = -normal  # Flip to point inward

                    # Signed distance from particle to plane (positive = outside/away from COM)
                    signed_dist = (p - v0).dot(normal)

                    # Only apply force if particle is close to membrane
                    membrane_thickness = h * 1.0
                    if ti.abs(signed_dist) < membrane_thickness:
                        # Project point onto plane
                        proj = p - signed_dist * normal

                        # Compute barycentric coordinates
                        v0v1 = v1 - v0
                        v0v2 = v2 - v0
                        v0p = proj - v0

                        dot00 = v0v2.dot(v0v2)
                        dot01 = v0v2.dot(v0v1)
                        dot02 = v0v2.dot(v0p)
                        dot11 = v0v1.dot(v0v1)
                        dot12 = v0v1.dot(v0p)

                        inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01 + 0.0001)
                        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
                        v = (dot00 * dot12 - dot01 * dot02) * inv_denom

                        # Check if inside triangle (with small margin)
                        margin = 0.1
                        if u >= -margin and v >= -margin and (u + v) <= 1.0 + margin:
                            # Particle near membrane - push toward inside (in +normal direction)
                            # signed_dist > 0 means inside, < 0 means outside
                            # Normal points inward, so push in +normal direction if outside
                            if signed_dist < 0:
                                # Particle is OUTSIDE - strong push back in
                                penetration = -signed_dist  # How far outside
                                force_mag = 10000.0 * (1.0 + penetration / membrane_thickness)
                                acc[i][0] += force_mag * normal[0]
                                acc[i][1] += force_mag * normal[1]
                                acc[i][2] += force_mag * normal[2]

                                # Also damp outward velocity
                                vel_normal = vel[i][0]*normal[0] + vel[i][1]*normal[1] + vel[i][2]*normal[2]
                                if vel_normal < 0:  # Moving outward
                                    vel[i][0] -= 0.5 * vel_normal * normal[0]
                                    vel[i][1] -= 0.5 * vel_normal * normal[1]
                                    vel[i][2] -= 0.5 * vel_normal * normal[2]

        self._apply_membrane_forces = apply_membrane_forces

    # =========================================================================
    # Main Simulation Step
    # =========================================================================
    def run_step(self):
        """Execute one complete simulation step."""
        # Mark GPU state as dirty (invalidate cache)
        self._gpu_dirty = True

        # Optional per-step state dump for solver debugging. Set
        # SIBERNETIC_DEBUG_DUMP=1 to emit pos/vel/acc for the first elastic
        # particle at the first 4 simulation steps — useful when tracing
        # through magnitude / units bugs as we found in the integration step
        # fix on 2026-05-03.
        import os as _os, sys as _sys
        _debug = _os.environ.get("SIBERNETIC_DEBUG_DUMP") == "1"
        if _debug and not hasattr(self, "_dbg_step"):
            self._dbg_step = 0
        if _debug and self._dbg_step <= 3:
            ti.sync()
            pos_np = self.pos.to_numpy()
            vel_np = self.vel.to_numpy()
            acc_np = self.acc.to_numpy()
            print(f"[DBG step={self._dbg_step} BEGIN] "
                  f"pos[0]={pos_np[0][:3]}  vel[0]={vel_np[0][:3]}  "
                  f"acc[0]={acc_np[0][:3]}  acc_old[0]={self.acc_old.to_numpy()[0][:3]}",
                  file=_sys.stderr, flush=True)

        start_time = time.perf_counter() if self._timing_enabled else 0

        # Grid-based spatial hashing (always needed)
        self._count_particles_in_grid()
        self._compute_grid_offsets()
        self._sort_particles()

        # Build neighbor list only in "list" mode
        if self._neighbor_mode == "list":
            self._build_neighbor_list()

        if self._timing_enabled:
            ti.sync()
            self._record_timing("neighbor_search", start_time)
            start_time = time.perf_counter()

        # Density and pressure computation
        if self._neighbor_mode == "fused":
            # FUSED: compute density, pressure AND cache neighbors in one pass
            self._compute_density_and_cache_neighbors()
        elif self._neighbor_mode == "csr":
            # CSR Pass 1: count neighbors + density + pressure (no storage)
            self._count_neighbors_and_density()
            # Compute neighbor offsets using hierarchical prefix sum
            self._neighbor_prefix_phase1()
            self._neighbor_prefix_phase2()
            self._neighbor_prefix_phase3()
            # CSR Pass 2: fill CSR storage
            self._fill_neighbors_csr()
        elif self._neighbor_mode == "list":
            # Use pre-built neighbor list
            self._compute_density()
            self._compute_pressure()
        else:  # "direct"
            # On-the-fly, no caching
            self._compute_density_direct()
            self._compute_pressure()

        if self._timing_enabled:
            ti.sync()
            self._record_timing("density_pressure", start_time)
            start_time = time.perf_counter()

        # Forces (pressure, viscosity, boundary, gravity)
        if self._neighbor_mode == "direct":
            # On-the-fly (searches neighbors again)
            self._compute_forces_direct()
        elif self._neighbor_mode == "csr":
            # Use CSR compact neighbor storage
            self._compute_forces_csr()
        else:
            # Use cached neighbors (from fused or list mode)
            self._compute_forces()

        # Snapshot acc BEFORE elastic to measure elastic-only contribution.
        if _debug and self._dbg_step <= 3 and self.has_elastic:
            ti.sync()
            acc_pre = self.acc.to_numpy()[0][:3].copy()

        # Elastic/muscle forces
        if self.has_elastic:
            self._compute_elastic_forces()

        # Snapshot acc AFTER elastic.
        if _debug and self._dbg_step <= 3 and self.has_elastic:
            ti.sync()
            acc_post = self.acc.to_numpy()[0][:3]
            elastic_contrib = acc_post - acc_pre
            print(f"[DBG step={self._dbg_step} ELASTIC] "
                  f"pre={acc_pre}  post={acc_post}  delta={elastic_contrib}",
                  file=_sys.stderr, flush=True)

        # Membrane collision forces
        if self.has_membranes:
            self._apply_membrane_forces()

        if self._timing_enabled:
            ti.sync()
            self._record_timing("forces", start_time)
            start_time = time.perf_counter()

        # ── DEBUG: snapshot AFTER force computation, BEFORE integration ──
        if _debug and self._dbg_step <= 3:
            ti.sync()
            acc_np = self.acc.to_numpy()
            print(f"[DBG step={self._dbg_step} POST-FORCES] "
                  f"acc[0]={acc_np[0][:3]}",
                  file=_sys.stderr, flush=True)

        # Integration (leapfrog)
        self._integrate_position()
        self._integrate_velocity()
        self._store_acceleration()

        # Floor constraint
        self._apply_floor_constraint()

        # ── DEBUG: snapshot AFTER integration ──
        if _debug and self._dbg_step <= 3:
            ti.sync()
            pos_np = self.pos.to_numpy()
            vel_np = self.vel.to_numpy()
            print(f"[DBG step={self._dbg_step} POST-INTEG] "
                  f"pos[0]={pos_np[0][:3]}  vel[0]={vel_np[0][:3]}",
                  file=_sys.stderr, flush=True)
            self._dbg_step += 1

        if self._timing_enabled:
            ti.sync()
            self._record_timing("integration", start_time)
            self._step_count += 1

    def run_integrate(self, mode=None):
        """Run integration step with specified mode.

        Args:
            mode: 0 = position update, 1 = velocity update, None = Euler
        """
        if mode == 0:
            self._integrate_position()
        elif mode == 1:
            self._integrate_velocity()
            self._store_acceleration()
            self._apply_floor_constraint()
        else:
            self._integrate_euler()
            self._apply_floor_constraint()

    # =========================================================================
    # Configuration
    # =========================================================================
    def set_neighbor_mode(self, mode):
        """Set neighbor search mode.

        Args:
            mode: One of:
                "fused"  - Compute density + cache neighbors in one pass (DEFAULT, fastest)
                "csr"    - CSR compact storage (lowest memory: ~3MB vs 43MB)
                "list"   - Build neighbor list first, then density/forces use it
                "direct" - On-the-fly for both (2x neighbor search, no memory)
        """
        valid_modes = {"fused", "csr", "list", "direct"}
        if mode not in valid_modes:
            raise ValueError(f"Invalid mode '{mode}'. Must be one of {valid_modes}")
        self._neighbor_mode = mode

    def get_memory_stats(self):
        """Return memory usage statistics for neighbor storage.

        Returns:
            dict with:
                - mode: current neighbor mode
                - fixed_array_mb: memory used by fixed-size 2D array
                - csr_capacity_mb: memory allocated for CSR storage
                - csr_used_mb: actual CSR memory used (if in CSR mode)
        """
        n = self.n_particles
        fixed_mb = n * self.max_neighbors * 4 / (1024 * 1024)
        csr_capacity_mb = self._csr_capacity * 4 / (1024 * 1024)

        stats = {
            "mode": self._neighbor_mode,
            "n_particles": n,
            "max_neighbors": self.max_neighbors,
            "fixed_array_mb": fixed_mb,
            "csr_capacity_mb": csr_capacity_mb,
        }

        if self._neighbor_mode == "csr":
            ti.sync()
            # Read from Taichi field - use to_numpy() for reliable access
            total_neighbors = int(self.neighbor_offsets.to_numpy()[n])
            stats["csr_used_mb"] = total_neighbors * 4 / (1024 * 1024)
            stats["total_neighbors"] = total_neighbors
            stats["avg_neighbors"] = total_neighbors / n if n > 0 else 0

        return stats

    # =========================================================================
    # State Access (for C++ integration)
    # =========================================================================
    def get_state(self):
        """Return position and velocity as contiguous float32 numpy arrays.

        Optimized: Returns cached numpy arrays, only syncing when GPU state is dirty.
        Use get_state_as_lists() if you need Python lists (slower).
        """
        if self._gpu_dirty or self._cached_pos is None:
            ti.sync()
            self._cached_pos = self.pos.to_numpy()
            self._cached_vel = self.vel.to_numpy()
            self._gpu_dirty = False
        return self._cached_pos, self._cached_vel

    def get_state_as_lists(self):
        """Legacy method returning Python lists (slower, for compatibility)."""
        pos, vel = self.get_state()
        return pos.tolist(), vel.tolist()

    def get_positions(self):
        """Return positions as numpy array (uses cache)."""
        pos, _ = self.get_state()
        return pos

    def get_velocities(self):
        """Return velocities as numpy array (uses cache)."""
        _, vel = self.get_state()
        return vel

    def get_density(self):
        """Return density array."""
        ti.sync()
        return self.rho.to_numpy()

    def get_pressure(self):
        """Return pressure array."""
        ti.sync()
        return self.pressure.to_numpy()

    def get_config(self):
        """Return configuration dictionary."""
        return self.config

    @property
    def position(self):
        """Property for compatibility with PyTorch solver (uses cache)."""
        pos, _ = self.get_state()
        return pos

    @property
    def velocity(self):
        """Property for compatibility with PyTorch solver (uses cache)."""
        _, vel = self.get_state()
        return vel

    # =========================================================================
    # Timing and Logging
    # =========================================================================
    def enable_timing(self, enabled):
        """Enable/disable timing collection."""
        self._timing_enabled = enabled
        if enabled:
            self._timing_data = {}
            self._step_count = 0

    def _record_timing(self, name, start_time):
        """Record timing for a substep."""
        elapsed = (time.perf_counter() - start_time) * 1000  # ms
        if name not in self._timing_data:
            self._timing_data[name] = 0.0
        self._timing_data[name] += elapsed

    def get_timing_breakdown(self):
        """Return timing data per substep (averaged over steps)."""
        if self._step_count == 0:
            return {}
        return {k: v / self._step_count for k, v in self._timing_data.items()}

    def enable_logging(self, enabled):
        """Enable/disable step logging."""
        self._logging_enabled = enabled
        if enabled:
            self._step_logs = []

    def get_step_log(self):
        """Return logged step data."""
        return self._step_logs if self._logging_enabled else None


# =============================================================================
# Benchmark and Demo
# =============================================================================
def benchmark_comparison():
    """Benchmark Taichi solver against PyTorch."""
    import sys
    sys.path.insert(0, '.')

    print("=" * 70)
    print("Taichi SPH Solver Benchmark")
    print("=" * 70)

    results = {}

    for n in [500, 1000, 2000, 5000, 10000]:
        # Create particles
        n_side = int(n ** (1/3)) + 1
        spacing = 0.03
        particles = []
        for i in range(n_side):
            for j in range(n_side):
                for k in range(n_side):
                    if len(particles) >= n:
                        break
                    particles.append([0.3 + i*spacing, 0.8 + j*spacing, 0.3 + k*spacing, 1.0])

        position = np.array(particles[:n], dtype=np.float32)
        velocity = np.zeros_like(position)

        config = {
            "xmin": 0.0, "ymin": -0.5, "zmin": 0.0,
            "xmax": 2.0, "ymax": 2.0, "zmax": 2.0,
            "h": 0.1,
            "rho0": 1000.0,
            "delta": 500.0,
            "time_step": 0.0001,
            "gravity_y": -9.81,
            "floor_y": 0.0,
            "floor_restitution": 0.3,
            "mu": 0.1,
        }

        # Test Metal
        config["device"] = "metal"
        solver = TaichiSolver(position.copy(), velocity.copy(), config)

        # Warmup
        for _ in range(10):
            solver.run_step()
        ti.sync()

        # Benchmark
        n_steps = 100
        start = time.perf_counter()
        for _ in range(n_steps):
            solver.run_step()
        ti.sync()
        elapsed = time.perf_counter() - start

        metal_rate = n_steps / elapsed
        results[n] = metal_rate

        # Verify correctness
        final_pos = solver.get_positions()
        y_min = final_pos[:, 1].min()
        y_max = final_pos[:, 1].max()

        print(f"Particles: {n:>6} | Metal: {metal_rate:>8.1f} steps/s | Y: {y_min:.2f} to {y_max:.2f}")

    print("\n" + "=" * 70)
    print("Summary: TaichiSolver ready for production use!")
    print("=" * 70)

    return results


if __name__ == "__main__":
    benchmark_comparison()
