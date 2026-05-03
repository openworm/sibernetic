#include <metal_stdlib>
using namespace metal;

// active × static pairwise squared distance.
// dist[i*n_static + j] = ||active[i] - static_p[j]||^2
//
// Why active×static (not all×all): demo1 has ~343 dynamic + ~17K static
// (boundary) particles. Static-static interactions never fire — the
// boundary doesn't move. Computing only active×static drops the matrix
// from 17K×17K (1.16 GB) to 343×17K (23 MB).
kernel void dist_active_static(
    device const packed_float3 *active   [[buffer(0)]],
    device const packed_float3 *static_p [[buffer(1)]],
    device float                *dist    [[buffer(2)]],
    constant uint              &n_active [[buffer(3)]],
    constant uint              &n_static [[buffer(4)]],
    uint2 gid                            [[thread_position_in_grid]])
{
    uint i = gid.x;  // active particle index
    uint j = gid.y;  // static particle index
    if (i >= n_active || j >= n_static) return;
    float3 a = float3(active[i]);
    float3 s = float3(static_p[j]);
    float3 d = a - s;
    dist[i * n_static + j] = d.x*d.x + d.y*d.y + d.z*d.z;
}

// M6.3 — active×active squared distance.
// Same pattern as dist_active_static but a single position buffer used
// as both rows and cols. Diagonal entries are 0.0 (self-distance);
// downstream consumers (e.g., density_constraint_grad) skip i==j.
//
// For demo1: 343×343 = 117K entries = 470 KB — tiny.
kernel void dist_active_active(
    device const packed_float3 *active   [[buffer(0)]],
    device float                *dist    [[buffer(1)]],
    constant uint              &n_active [[buffer(2)]],
    uint2 gid                            [[thread_position_in_grid]])
{
    uint i = gid.x;
    uint j = gid.y;
    if (i >= n_active || j >= n_active) return;
    float3 a = float3(active[i]);
    float3 b = float3(active[j]);
    float3 d = a - b;
    dist[i * n_active + j] = d.x*d.x + d.y*d.y + d.z*d.z;
}

// M6.1 — Apply Müller 2003 Wpoly6 SPH kernel elementwise on a squared
// distance matrix. In-place: input is r², output is W(r²) overwriting
// the same buffer.
//
//     W_poly6(r²) = poly6_const · (h² - r²)³   if r < h else 0
//     poly6_const = 315 / (64π h⁹)             (computed host-side)
//
// One thread per matrix element. Trivially parallel, no atomics.
kernel void wpoly6_inplace(
    device float       *r2_or_W   [[buffer(0)]],
    constant float     &h2         [[buffer(1)]],
    constant float     &poly6_const [[buffer(2)]],
    constant uint      &n_total    [[buffer(3)]],
    uint gid                        [[thread_position_in_grid]])
{
    if (gid >= n_total) return;
    float r2 = r2_or_W[gid];
    if (r2 >= h2) {
        r2_or_W[gid] = 0.0;
        return;
    }
    float diff = h2 - r2;
    r2_or_W[gid] = poly6_const * diff * diff * diff;
}

// M6.2 — Density via row-sum reduction.
//
//     density[i] = mass · Σ_j W[i, j]
//
// Threadgroup-level reduction: each row is one threadgroup of TG=256
// threads. Each thread strides through the row, accumulates a partial
// sum, then a tree reduction in threadgroup memory combines partials.
// For demo1 (n_rows=343) this dispatches 87,808 threads — saturates
// the Apple GPU. The naive "one thread per row" version was the
// bottleneck of the M6 SPH pipeline (3.65 ms vs ~0.5 ms expected).
//
// Dispatch: grid (256, n_rows, 1), threadgroup (256, 1, 1).
kernel void rowsum_density(
    device const float *W        [[buffer(0)]],
    device float       *density  [[buffer(1)]],
    constant float     &mass     [[buffer(2)]],
    constant uint      &n_cols   [[buffer(3)]],
    constant uint      &n_rows   [[buffer(4)]],
    uint2 gid                     [[thread_position_in_grid]],
    uint2 lid                     [[thread_position_in_threadgroup]],
    uint2 tg_size                 [[threads_per_threadgroup]])
{
    uint i = gid.y;
    if (i >= n_rows) return;

    threadgroup float partials[256];
    uint t = lid.x;
    uint T = tg_size.x;

    float sum = 0.0;
    for (uint j = t; j < n_cols; j += T) {
        sum += W[i * n_cols + j];
    }
    partials[t] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint stride = T / 2; stride > 0; stride >>= 1) {
        if (t < stride) {
            partials[t] += partials[t + stride];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (t == 0) {
        density[i] = mass * partials[0];
    }
}

// PERF M6.4-fused — density + grad_C + denom_helper in ONE kernel.
// The standard XPBD pipeline computes density, grad_C, and denom_helper
// from the SAME per-pair iterator over (active+static) neighbors. This
// kernel does all three reductions in one threadgroup pass over r².
//
// Replaces (per inner XPBD iter):
//   wpoly6_rowsum_density_fused × 2 (aa + as)
//   density_constraint_grad
// → single density_grad_combined kernel.
//
// Outputs:
//   density[i]      = mass · Σ_j W_poly6(r_ij)
//   grad_C[i]       = (mass/ρ_rest) · Σ_j ∇W_spiky(p_i - p_j)
//   denom_helper[i] = Σ_j |∇W_spiky(p_i - p_j)|²
//
// Three threadgroup reduction buffers (float for density, float3 for
// grad, float for denom). Apple Silicon has plenty of threadgroup mem.
//
// Dispatch: grid (256, n_active, 1), threadgroup (256, 1, 1).
kernel void density_grad_combined(
    device const packed_float3  *active        [[buffer(0)]],
    device const packed_float3  *static_p      [[buffer(1)]],
    device const float          *r2_aa         [[buffer(2)]],
    device const float          *r2_as         [[buffer(3)]],
    device float                *density       [[buffer(4)]],   // out
    device packed_float3        *grad_C        [[buffer(5)]],   // out
    device float                *denom_helper  [[buffer(6)]],   // out
    constant float              &h             [[buffer(7)]],
    constant float              &poly6_const   [[buffer(8)]],
    constant float              &spiky_const   [[buffer(9)]],
    constant float              &mass          [[buffer(10)]],
    constant float              &rho_rest      [[buffer(11)]],
    constant uint               &n_active      [[buffer(12)]],
    constant uint               &n_static      [[buffer(13)]],
    uint2 gid                                   [[thread_position_in_grid]],
    uint2 lid                                   [[thread_position_in_threadgroup]],
    uint2 tg_size                               [[threads_per_threadgroup]])
{
    uint i = gid.y;
    if (i >= n_active) return;

    threadgroup float  partials_dens[256];
    threadgroup float3 partials_grad[256];
    threadgroup float  partials_denom[256];
    uint t = lid.x;
    uint T = tg_size.x;

    float h2 = h * h;
    float3 p_i = float3(active[i]);
    float  partial_dens  = 0.0;
    float3 partial_grad  = float3(0.0);
    float  partial_denom = 0.0;
    uint n_total = n_active + n_static;

    for (uint k = t; k < n_total; k += T) {
        float r2;
        float3 dir;
        if (k < n_active) {
            uint j = k;
            // Self contribution to density (W(0)) IS included in the
            // standard SPH formulation but NOT for grad/denom (avoid
            // singularity at r=0).
            r2 = r2_aa[i * n_active + j];
            if (r2 >= h2) continue;
            // Wpoly6 contribution to density (always, including self).
            float diff = h2 - r2;
            partial_dens += poly6_const * diff * diff * diff;
            if (j == i) continue;  // skip self for grad/denom only
            dir = p_i - float3(active[j]);
        } else {
            uint j = k - n_active;
            r2 = r2_as[i * n_static + j];
            if (r2 >= h2) continue;
            float diff = h2 - r2;
            partial_dens += poly6_const * diff * diff * diff;
            dir = p_i - float3(static_p[j]);
        }
        float r = sqrt(r2);
        float h_minus_r = h - r;
        float coef = spiky_const * h_minus_r * h_minus_r / (r + 1e-7);
        float3 grad_W = coef * dir;
        partial_grad  += mass * grad_W;
        partial_denom += grad_W.x*grad_W.x + grad_W.y*grad_W.y + grad_W.z*grad_W.z;
    }
    partials_dens[t]  = partial_dens;
    partials_grad[t]  = partial_grad;
    partials_denom[t] = partial_denom;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint stride = T / 2; stride > 0; stride >>= 1) {
        if (t < stride) {
            partials_dens[t]  += partials_dens[t + stride];
            partials_grad[t]  += partials_grad[t + stride];
            partials_denom[t] += partials_denom[t + stride];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (t == 0) {
        density[i] = mass * partials_dens[0];
        grad_C[i] = packed_float3(partials_grad[0] / rho_rest);
        denom_helper[i] = partials_denom[0];
    }
}

// PERF MEGA-GRID — same as mega_fused but uses a precomputed spatial
// hash grid for static neighbors. Each active particle looks up its
// 3×3×3 cell neighborhood (27 cells) instead of iterating all 17498
// static particles. At demo1's particle density, ~99.8% of pair
// checks are wasted; grid lookup recovers that work.
//
// Grid buffers (built host-side, one-time, written by load_config.py):
//   sorted_static[i_sorted]       — static particle positions, reordered
//                                    so particles in the same cell are
//                                    contiguous in memory
//   cell_start[c]                  — first index into sorted_static for
//                                    cell c. Range [cell_start[c],
//                                    cell_start[c+1]) is the cell's
//                                    particles. Length = n_cells + 1.
// Grid params (constants): grid_dim (int3), grid_origin (float3), h.
//
// Active-active still uses dense iteration (n_active = 343 — N² is
// ~117K pair checks, fine without a grid).
kernel void density_grad_mega_grid(
    device const packed_float3 *active        [[buffer(0)]],
    device const packed_float3 *sorted_static [[buffer(1)]],
    device const int           *cell_start    [[buffer(2)]],
    device float                *density      [[buffer(3)]],
    device packed_float3        *grad_C       [[buffer(4)]],
    device float                *denom_helper [[buffer(5)]],
    constant float             &h             [[buffer(6)]],
    constant float             &h2            [[buffer(7)]],
    constant float             &poly6_const   [[buffer(8)]],
    constant float             &spiky_const   [[buffer(9)]],
    constant float             &mass          [[buffer(10)]],
    constant float             &rho_rest      [[buffer(11)]],
    constant uint              &n_active      [[buffer(12)]],
    constant int3              &grid_dim      [[buffer(13)]],
    constant packed_float3     &grid_origin   [[buffer(14)]],
    uint2 gid                                  [[thread_position_in_grid]],
    uint2 lid                                  [[thread_position_in_threadgroup]],
    uint2 tg_size                              [[threads_per_threadgroup]])
{
    uint i = gid.y;
    if (i >= n_active) return;

    uint t = lid.x;
    uint T = tg_size.x;

    float3 p_i = float3(active[i]);
    float  partial_dens  = 0.0;
    float3 partial_grad  = float3(0.0);
    float  partial_denom = 0.0;

    // ── active-active (dense) ──
    for (uint k = t; k < n_active; k += T) {
        float3 p_j = float3(active[k]);
        float3 dir = p_i - p_j;
        float r2 = dot(dir, dir);
        if (r2 >= h2) continue;
        float diff = h2 - r2;
        partial_dens += poly6_const * diff * diff * diff;
        if (k == i) continue;
        float r = sqrt(r2);
        float h_minus_r = h - r;
        float coef = spiky_const * h_minus_r * h_minus_r / (r + 1e-7);
        float3 grad_W = coef * dir;
        partial_grad  += mass * grad_W;
        partial_denom += grad_W.x*grad_W.x + grad_W.y*grad_W.y + grad_W.z*grad_W.z;
    }

    // ── static via 27-cell grid lookup ──
    int3 my_cell = int3(floor((p_i - float3(grid_origin)) / h));
    int n_cells_xy = grid_dim.x * grid_dim.y;
    for (int dz = -1; dz <= 1; dz++) {
        int cz = my_cell.z + dz;
        if (cz < 0 || cz >= grid_dim.z) continue;
        for (int dy = -1; dy <= 1; dy++) {
            int cy = my_cell.y + dy;
            if (cy < 0 || cy >= grid_dim.y) continue;
            for (int dx = -1; dx <= 1; dx++) {
                int cx = my_cell.x + dx;
                if (cx < 0 || cx >= grid_dim.x) continue;
                int c_id = cx + cy * grid_dim.x + cz * n_cells_xy;
                int start = cell_start[c_id];
                int end   = cell_start[c_id + 1];
                // 256 threads cooperatively iterate over [start, end).
                for (int j = start + (int)t; j < end; j += (int)T) {
                    float3 p_j = float3(sorted_static[j]);
                    float3 dir = p_i - p_j;
                    float r2 = dot(dir, dir);
                    if (r2 >= h2) continue;
                    float diff = h2 - r2;
                    partial_dens += poly6_const * diff * diff * diff;
                    float r = sqrt(r2);
                    float h_minus_r = h - r;
                    float coef = spiky_const * h_minus_r * h_minus_r / (r + 1e-7);
                    float3 grad_W = coef * dir;
                    partial_grad  += mass * grad_W;
                    partial_denom += grad_W.x*grad_W.x
                                   + grad_W.y*grad_W.y
                                   + grad_W.z*grad_W.z;
                }
            }
        }
    }

    // simdgroup reduction. Number of simdgroups = T/32. We size the
    // arrays for max 8 (i.e. 256-thread TG) and only iterate the
    // populated entries — keeps the kernel threadgroup-size-agnostic.
    float  s_dens  = simd_sum(partial_dens);
    float3 s_grad  = float3(simd_sum(partial_grad.x),
                            simd_sum(partial_grad.y),
                            simd_sum(partial_grad.z));
    float  s_denom = simd_sum(partial_denom);

    threadgroup float  simd_dens_arr[8];
    threadgroup float3 simd_grad_arr[8];
    threadgroup float  simd_denom_arr[8];
    uint n_simds = (T + 31) / 32;
    uint simd_id = t / 32;
    uint lane = t % 32;
    if (lane == 0 && simd_id < 8) {
        simd_dens_arr[simd_id]  = s_dens;
        simd_grad_arr[simd_id]  = s_grad;
        simd_denom_arr[simd_id] = s_denom;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (t == 0) {
        float total_d = 0.0;
        float3 total_g = float3(0.0);
        float total_h = 0.0;
        for (uint s = 0; s < n_simds; s++) {
            total_d += simd_dens_arr[s];
            total_g += simd_grad_arr[s];
            total_h += simd_denom_arr[s];
        }
        density[i] = mass * total_d;
        grad_C[i] = packed_float3(total_g / rho_rest);
        denom_helper[i] = total_h;
    }
}

// ──────────────────────────────────────────────────────────────────────
// pair_forces_grid — viscosity + surface tension as a single fused
// SPH pair-sum kernel. Mirrors Sibernetic's `pcisph_computeForces…`
// (src/sphFluid.cl:519) but specialized to demo1's particle types
// (cube = elastic 2.x, liquid = 1.x, walls = boundary 3.x).
//
// Two physics components per neighbor pair within r < h:
//
//   visc_force_i = (1.5 · m · ∇²W_visc / ρ_i)
//                  · Σ_j coef_pair · (v_j - v_i) · (h_s - r_s) / 1000
//   surf_force_i = (-1.7e-9 · surfTensCoeff / m)
//                  · Σ_j (h_s² - r_s²)³ · (p_i - p_j)_units
//
// where subscript _s means scaled to meters (× sim_scale). The
// constant 1000 is Sibernetic's literal — replicating their formula.
// `coef_pair = 1e-4` for our cube case (elastic↔elastic,
// elastic↔liquid, elastic↔boundary).
//
// Boundary particles have v_j = 0 (they don't move) — their visc
// contribution is `coef_pair · (-v_i) · (h_s - r_s) / 1000`, which
// is the drag-on-impact force that splat-stretches OpenCL's cube.
//
// Output is per-active-particle external acceleration in m/s² —
// add it to gravity in the predict step.
//
// Differentiability note: this kernel reads (pos, vel) and writes
// ext_accel. A matching backward kernel would chain ∂L/∂ext_accel
// back to ∂L/∂pos and ∂L/∂vel; not needed for the trajectory-match
// task but listed in M-plan for the trainable-with-viscosity path.
// ──────────────────────────────────────────────────────────────────────
kernel void pair_forces_grid(
    device const packed_float3 *active_pos      [[buffer(0)]],
    device const packed_float3 *active_vel      [[buffer(1)]],
    device const packed_float3 *sorted_static   [[buffer(2)]],
    device const int           *cell_start      [[buffer(3)]],
    device const float         *density         [[buffer(4)]],
    device packed_float3       *ext_accel       [[buffer(5)]],
    constant float             &h               [[buffer(6)]],
    constant float             &h2              [[buffer(7)]],
    constant float             &mass            [[buffer(8)]],
    constant float             &sim_scale       [[buffer(9)]],
    constant float             &visc_pair_coef  [[buffer(10)]],
    constant float             &visc_amp        [[buffer(11)]],
    constant float             &surf_amp        [[buffer(12)]],
    constant uint              &n_active        [[buffer(13)]],
    constant int3              &grid_dim        [[buffer(14)]],
    constant packed_float3     &grid_origin     [[buffer(15)]],
    uint2 gid                                    [[thread_position_in_grid]],
    uint2 lid                                    [[thread_position_in_threadgroup]],
    uint2 tg_size                                [[threads_per_threadgroup]])
{
    uint i = gid.y;
    if (i >= n_active) return;

    uint t = lid.x;
    uint T = tg_size.x;

    float3 p_i = float3(active_pos[i]);
    float3 v_i = float3(active_vel[i]);
    float h_scaled  = h * sim_scale;
    float h2_scaled = h_scaled * h_scaled;
    float ss2       = sim_scale * sim_scale;

    float3 visc_partial = float3(0.0);
    float3 surf_partial = float3(0.0);

    // ── active-active (dense) ──
    for (uint k = t; k < n_active; k += T) {
        if (k == i) continue;
        float3 p_j = float3(active_pos[k]);
        float3 v_j = float3(active_vel[k]);
        float3 dir_unit = p_i - p_j;
        float r2_unit = dot(dir_unit, dir_unit);
        if (r2_unit >= h2) continue;
        float r_unit = sqrt(r2_unit);
        float r_scaled = r_unit * sim_scale;
        float h_minus_r = h_scaled - r_scaled;
        float r2_scaled = r2_unit * ss2;
        float surf_kern = (h2_scaled - r2_scaled);
        surf_kern = surf_kern * surf_kern * surf_kern;

        visc_partial += visc_pair_coef * (v_j - v_i) * h_minus_r * (1.0/1000.0);
        surf_partial += surf_kern * dir_unit;
    }

    // ── static (boundary) via 27-cell grid lookup. v_j ≡ 0. ──
    int3 my_cell = int3(floor((p_i - float3(grid_origin)) / h));
    int n_cells_xy = grid_dim.x * grid_dim.y;
    for (int dz = -1; dz <= 1; dz++) {
        int cz = my_cell.z + dz;
        if (cz < 0 || cz >= grid_dim.z) continue;
        for (int dy = -1; dy <= 1; dy++) {
            int cy = my_cell.y + dy;
            if (cy < 0 || cy >= grid_dim.y) continue;
            for (int dx = -1; dx <= 1; dx++) {
                int cx = my_cell.x + dx;
                if (cx < 0 || cx >= grid_dim.x) continue;
                int c_id = cx + cy * grid_dim.x + cz * n_cells_xy;
                int start = cell_start[c_id];
                int end   = cell_start[c_id + 1];
                for (int j = start + (int)t; j < end; j += (int)T) {
                    float3 p_j = float3(sorted_static[j]);
                    float3 dir_unit = p_i - p_j;
                    float r2_unit = dot(dir_unit, dir_unit);
                    if (r2_unit >= h2) continue;
                    float r_unit = sqrt(r2_unit);
                    float r_scaled = r_unit * sim_scale;
                    float h_minus_r = h_scaled - r_scaled;
                    float r2_scaled = r2_unit * ss2;
                    float surf_kern = (h2_scaled - r2_scaled);
                    surf_kern = surf_kern * surf_kern * surf_kern;

                    // Boundary velocity is zero → (v_j - v_i) = -v_i.
                    visc_partial += visc_pair_coef * (-v_i) * h_minus_r * (1.0/1000.0);
                    surf_partial += surf_kern * dir_unit;
                }
            }
        }
    }

    // simdgroup reduction over T threads (32, 64, 128, or 256).
    float3 s_visc = float3(simd_sum(visc_partial.x),
                           simd_sum(visc_partial.y),
                           simd_sum(visc_partial.z));
    float3 s_surf = float3(simd_sum(surf_partial.x),
                           simd_sum(surf_partial.y),
                           simd_sum(surf_partial.z));

    threadgroup float3 simd_visc_arr[8];
    threadgroup float3 simd_surf_arr[8];
    uint n_simds = (T + 31) / 32;
    uint simd_id = t / 32;
    uint lane = t % 32;
    if (lane == 0 && simd_id < 8) {
        simd_visc_arr[simd_id] = s_visc;
        simd_surf_arr[simd_id] = s_surf;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (t == 0) {
        float3 v_total = float3(0.0);
        float3 s_total = float3(0.0);
        for (uint s = 0; s < n_simds; s++) {
            v_total += simd_visc_arr[s];
            s_total += simd_surf_arr[s];
        }
        // Final scaling. For density-zero rows (shouldn't happen at
        // settled state but might at step 0), clamp to avoid 1/0.
        float rho_i = max(density[i], 1e-30);
        float3 a_visc = (visc_amp / rho_i) * v_total;
        float3 a_surf = surf_amp * s_total;
        ext_accel[i] = packed_float3(a_visc + a_surf);
    }
}

// apply_ext_accel — adds external acceleration to velocity:
//   v += dt · a_ext   (separate from gravity; gravity is in predict)
// Used between pair_forces_grid (computes ext_accel) and predict_positions.
kernel void apply_ext_accel(
    device packed_float3       *vel        [[buffer(0)]],
    device const packed_float3 *ext_accel  [[buffer(1)]],
    constant float             &dt         [[buffer(2)]],
    constant uint              &n          [[buffer(3)]],
    uint gid                                [[thread_position_in_grid]])
{
    if (gid >= n) return;
    vel[gid] = packed_float3(float3(vel[gid]) + float3(ext_accel[gid]) * dt);
}

// ──────────────────────────────────────────────────────────────────────
// pair_forces_grid_backward — gradients of (visc + surf) ext_accel
// w.r.t. positions and velocities.
//
// Setup: forward computes ext_accel_i = visc_amp/ρ_i · visc_partial_i
//                                    + surf_amp · surf_partial_i
// where the partials are pair sums over r < h neighbors.
//
// Symmetric backward trick: each thread (for particle i) reads both
//   - its OWN upstream ga_i = ∂L/∂(ext_accel_i), and
//   - the upstream of every neighbor ga_j (active only; boundary frozen)
// then accumulates contributions to ONLY ∂L/∂p_i and ∂L/∂v_i. The
// cross-particle gradients ∂L/∂p_j are picked up when thread j runs
// its own loop and sees i as a neighbor. No atomics needed.
//
// Per-pair (i, j) gradient contributions to ∂L/∂p_i:
//   visc:  K_v / r · <G_v_j − G_v_i, v_diff> · dir
//   surf:  −6·ss²·(h_s²−r_s²)² · dir · <dir, G_s_i − G_s_j>
//        + surf_kern · (G_s_i − G_s_j)
// where dir = p_i − p_j, v_diff = v_j − v_i,
//   K_v = visc_pair_coef · sim_scale / 1000,
//   G_v_a = (visc_amp / ρ_a) · ga_a,
//   G_s_a = surf_amp · ga_a.
// For boundary j: G_v_j = 0, G_s_j = 0, v_j = 0 (so v_diff = −v_i).
// Sign flips fall out cleanly.
//
// Per-pair contributions to ∂L/∂v_i:
//   visc:  K_v_v · (G_v_j − G_v_i)
// where K_v_v = visc_pair_coef · (h_s − r_s) / 1000. (No surf v-dep.)
//
// Note: density ρ_i is treated as a stop-gradient (we don't backprop
// through density-compute via this path). For tasks where you train
// fluid params through trajectory loss, this is sufficient.
// ──────────────────────────────────────────────────────────────────────
kernel void pair_forces_grid_backward(
    device const packed_float3 *active_pos       [[buffer(0)]],
    device const packed_float3 *active_vel       [[buffer(1)]],
    device const packed_float3 *sorted_static    [[buffer(2)]],
    device const int           *cell_start       [[buffer(3)]],
    device const float         *density          [[buffer(4)]],
    device const packed_float3 *grad_ext_accel   [[buffer(5)]],
    device packed_float3       *grad_pos         [[buffer(6)]],
    device packed_float3       *grad_vel         [[buffer(7)]],
    constant float             &h                [[buffer(8)]],
    constant float             &h2               [[buffer(9)]],
    constant float             &mass             [[buffer(10)]],
    constant float             &sim_scale        [[buffer(11)]],
    constant float             &visc_pair_coef   [[buffer(12)]],
    constant float             &visc_amp         [[buffer(13)]],
    constant float             &surf_amp         [[buffer(14)]],
    constant uint              &n_active         [[buffer(15)]],
    constant int3              &grid_dim         [[buffer(16)]],
    constant packed_float3     &grid_origin      [[buffer(17)]],
    uint gid                                      [[thread_position_in_grid]])
{
    if (gid >= n_active) return;
    uint i = gid;

    float3 p_i = float3(active_pos[i]);
    float3 v_i = float3(active_vel[i]);
    float3 ga_i = float3(grad_ext_accel[i]);
    float rho_i = max(density[i], 1e-30);

    float h_scaled  = h * sim_scale;
    float h2_scaled = h_scaled * h_scaled;
    float ss2       = sim_scale * sim_scale;

    float3 G_v_i = (visc_amp / rho_i) * ga_i;
    float3 G_s_i = surf_amp * ga_i;

    float3 dp_i = float3(0.0);
    float3 dv_i = float3(0.0);

    // ── active-active loop ──
    for (uint k = 0; k < n_active; k++) {
        if (k == i) continue;
        float3 p_j = float3(active_pos[k]);
        float3 v_j = float3(active_vel[k]);
        float3 dir = p_i - p_j;
        float r2_unit = dot(dir, dir);
        if (r2_unit >= h2) continue;
        float r_unit = sqrt(r2_unit);
        if (r_unit < 1e-7) continue;

        float r_scaled = r_unit * sim_scale;
        float h_minus_r = h_scaled - r_scaled;
        float r2_scaled = r2_unit * ss2;
        float h2r2 = h2_scaled - r2_scaled;
        float surf_kern = h2r2 * h2r2 * h2r2;
        float surf_kern_deriv = h2r2 * h2r2;
        float3 v_diff = v_j - v_i;

        float rho_j = max(density[k], 1e-30);
        float3 ga_j = float3(grad_ext_accel[k]);
        float3 G_v_j = (visc_amp / rho_j) * ga_j;
        float3 G_s_j = surf_amp * ga_j;

        float3 G_v_diff = G_v_j - G_v_i;
        float3 G_s_diff = G_s_i - G_s_j;

        // viscosity ∂L/∂p_i
        dp_i += (visc_pair_coef * sim_scale / (1000.0 * r_unit))
              * dot(G_v_diff, v_diff) * dir;
        // surface tension ∂L/∂p_i
        dp_i += -6.0 * ss2 * surf_kern_deriv * dir * dot(dir, G_s_diff)
              +  surf_kern * G_s_diff;
        // viscosity ∂L/∂v_i
        dv_i += (visc_pair_coef * h_minus_r / 1000.0) * G_v_diff;
    }

    // ── active-static loop (boundary j, v_j = 0, ga_j = 0) ──
    int3 my_cell = int3(floor((p_i - float3(grid_origin)) / h));
    int n_cells_xy = grid_dim.x * grid_dim.y;
    for (int dz = -1; dz <= 1; dz++) {
        int cz = my_cell.z + dz;
        if (cz < 0 || cz >= grid_dim.z) continue;
        for (int dy = -1; dy <= 1; dy++) {
            int cy = my_cell.y + dy;
            if (cy < 0 || cy >= grid_dim.y) continue;
            for (int dx = -1; dx <= 1; dx++) {
                int cx = my_cell.x + dx;
                if (cx < 0 || cx >= grid_dim.x) continue;
                int c_id = cx + cy * grid_dim.x + cz * n_cells_xy;
                int start = cell_start[c_id];
                int end   = cell_start[c_id + 1];
                for (int j = start; j < end; j++) {
                    float3 p_j = float3(sorted_static[j]);
                    float3 dir = p_i - p_j;
                    float r2_unit = dot(dir, dir);
                    if (r2_unit >= h2) continue;
                    float r_unit = sqrt(r2_unit);
                    if (r_unit < 1e-7) continue;

                    float r_scaled = r_unit * sim_scale;
                    float h_minus_r = h_scaled - r_scaled;
                    float r2_scaled = r2_unit * ss2;
                    float h2r2 = h2_scaled - r2_scaled;
                    float surf_kern = h2r2 * h2r2 * h2r2;
                    float surf_kern_deriv = h2r2 * h2r2;
                    float3 v_diff = -v_i;  // v_j = 0

                    // Boundary: G_v_j = 0, G_s_j = 0
                    float3 G_v_diff = -G_v_i;
                    float3 G_s_diff = G_s_i;

                    dp_i += (visc_pair_coef * sim_scale / (1000.0 * r_unit))
                          * dot(G_v_diff, v_diff) * dir;
                    dp_i += -6.0 * ss2 * surf_kern_deriv * dir * dot(dir, G_s_diff)
                          +  surf_kern * G_s_diff;
                    dv_i += (visc_pair_coef * h_minus_r / 1000.0) * G_v_diff;
                }
            }
        }
    }

    grad_pos[i] = packed_float3(float3(grad_pos[i]) + dp_i);
    grad_vel[i] = packed_float3(float3(grad_vel[i]) + dv_i);
}

// apply_ext_accel_backward — straight-through for v, scale by dt for ext_accel.
// Forward: v_new = v_old + dt · ext_accel
//   dL/d(v_old)     += dL/d(v_new)
//   dL/d(ext_accel) += dt · dL/d(v_new)
kernel void apply_ext_accel_backward(
    device const packed_float3 *grad_v_new        [[buffer(0)]],
    device packed_float3       *grad_v_old        [[buffer(1)]],
    device packed_float3       *grad_ext_accel    [[buffer(2)]],
    constant float             &dt                [[buffer(3)]],
    constant uint              &n                 [[buffer(4)]],
    uint gid                                       [[thread_position_in_grid]])
{
    if (gid >= n) return;
    float3 g_new = float3(grad_v_new[gid]);
    grad_v_old[gid]     = packed_float3(float3(grad_v_old[gid]) + g_new);
    grad_ext_accel[gid] = packed_float3(float3(grad_ext_accel[gid]) + dt * g_new);
}

// PERF MEGA-FUSED — distance + density + grad_C + denom_helper, all
// inline. NO r² matrix materialization. Each inner XPBD iter reads
// positions fresh and computes everything; for our demo1 scale this
// drops memory bandwidth from 96 MB/step to ~750 KB/step, even though
// per-pair compute grows ~30%.
//
// As a bonus: each inner iter uses CURRENT positions for distance
// (fresh from constraint corrections), so this is strictly more
// accurate than distance-reuse across iters.
//
// Tradeoff: r² is no longer available for backward kernels. Use
// density_grad_mega_fused only in forward-only paths (xpbd_step).
// xpbd_full_fwd (which saves state for backward) keeps using the
// non-mega path so r² is recomputable in backward.
kernel void density_grad_mega_fused(
    device const packed_float3 *active        [[buffer(0)]],
    device const packed_float3 *static_p      [[buffer(1)]],
    device float                *density      [[buffer(2)]],
    device packed_float3        *grad_C       [[buffer(3)]],
    device float                *denom_helper [[buffer(4)]],
    constant float             &h             [[buffer(5)]],
    constant float             &h2            [[buffer(6)]],
    constant float             &poly6_const   [[buffer(7)]],
    constant float             &spiky_const   [[buffer(8)]],
    constant float             &mass          [[buffer(9)]],
    constant float             &rho_rest      [[buffer(10)]],
    constant uint              &n_active      [[buffer(11)]],
    constant uint              &n_static      [[buffer(12)]],
    uint2 gid                                  [[thread_position_in_grid]],
    uint2 lid                                  [[thread_position_in_threadgroup]],
    uint2 tg_size                              [[threads_per_threadgroup]])
{
    uint i = gid.y;
    if (i >= n_active) return;

    uint t = lid.x;
    uint T = tg_size.x;

    float3 p_i = float3(active[i]);
    float  partial_dens  = 0.0;
    float3 partial_grad  = float3(0.0);
    float  partial_denom = 0.0;
    uint n_total = n_active + n_static;

    for (uint k = t; k < n_total; k += T) {
        float3 p_j;
        bool is_active_neighbor = (k < n_active);
        if (is_active_neighbor) {
            p_j = float3(active[k]);
        } else {
            p_j = float3(static_p[k - n_active]);
        }
        // Inline distance.
        float3 dir = p_i - p_j;
        float r2 = dot(dir, dir);
        if (r2 >= h2) continue;
        // Wpoly6 → density (always, including self at r²=0).
        float diff = h2 - r2;
        partial_dens += poly6_const * diff * diff * diff;
        // Skip self for grad/denom (avoid 0/0 in r̂).
        if (is_active_neighbor && k == i) continue;
        float r = sqrt(r2);
        float h_minus_r = h - r;
        float coef = spiky_const * h_minus_r * h_minus_r / (r + 1e-7);
        float3 grad_W = coef * dir;
        partial_grad  += mass * grad_W;
        partial_denom += grad_W.x*grad_W.x + grad_W.y*grad_W.y + grad_W.z*grad_W.z;
    }
    // PERF — simdgroup-aware reduction. Apple Silicon GPU has 32-lane
    // SIMD groups; simd_sum is a single hardware instruction. We reduce
    // within each simdgroup (32 threads → 1 result) then combine the 8
    // simdgroup results in a quick threadgroup pass. This replaces the
    // 8-iteration tree reduction.
    float  s_dens  = simd_sum(partial_dens);
    float3 s_grad  = float3(simd_sum(partial_grad.x),
                            simd_sum(partial_grad.y),
                            simd_sum(partial_grad.z));
    float  s_denom = simd_sum(partial_denom);

    threadgroup float  simd_dens_arr[8];
    threadgroup float3 simd_grad_arr[8];
    threadgroup float  simd_denom_arr[8];
    uint simd_id = t / 32;
    uint lane = t % 32;
    if (lane == 0) {
        simd_dens_arr[simd_id]  = s_dens;
        simd_grad_arr[simd_id]  = s_grad;
        simd_denom_arr[simd_id] = s_denom;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (t == 0) {
        float total_d = 0.0;
        float3 total_g = float3(0.0);
        float total_h = 0.0;
        for (uint s = 0; s < 8; s++) {
            total_d += simd_dens_arr[s];
            total_g += simd_grad_arr[s];
            total_h += simd_denom_arr[s];
        }
        density[i] = mass * total_d;
        grad_C[i] = packed_float3(total_g / rho_rest);
        denom_helper[i] = total_h;
    }
}

// PERF M6.1+M6.2 fused — apply Wpoly6 to r² AND row-reduce in one pass.
// Skips materializing the W matrix (saves ~24 MB of memory bandwidth at
// demo1's 343×17498 = 6M-entry W matrix). Same threadgroup-reduction
// pattern as rowsum_density, but the Wpoly6 evaluation happens inline
// inside the per-thread accumulator loop.
//
// Forward chain replaced:
//   wpoly6_inplace + rowsum_density   →   wpoly6_rowsum_density_fused
// Back-compat: the un-fused kernels stay so backward (which separates
// the two stages) still works unchanged.
//
// Dispatch: grid (256, n_rows, 1), threadgroup (256, 1, 1).
kernel void wpoly6_rowsum_density_fused(
    device const float *r2          [[buffer(0)]],   // input
    device float       *density     [[buffer(1)]],   // output
    constant float     &h2          [[buffer(2)]],
    constant float     &poly6_const [[buffer(3)]],
    constant float     &mass        [[buffer(4)]],
    constant uint      &n_cols      [[buffer(5)]],
    constant uint      &n_rows      [[buffer(6)]],
    uint2 gid                        [[thread_position_in_grid]],
    uint2 lid                        [[thread_position_in_threadgroup]],
    uint2 tg_size                    [[threads_per_threadgroup]])
{
    uint i = gid.y;
    if (i >= n_rows) return;

    threadgroup float partials[256];
    uint t = lid.x;
    uint T = tg_size.x;

    float partial = 0.0;
    for (uint j = t; j < n_cols; j += T) {
        float r2_val = r2[i * n_cols + j];
        if (r2_val < h2) {
            float diff = h2 - r2_val;
            partial += poly6_const * diff * diff * diff;
        }
    }
    partials[t] = partial;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint stride = T / 2; stride > 0; stride >>= 1) {
        if (t < stride) {
            partials[t] += partials[t + stride];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (t == 0) {
        density[i] = mass * partials[0];
    }
}

// M6.4 — Density constraint gradient (fused).
//
// XPBD's density constraint is C_i = ρ_i / ρ_rest - 1. We need
// ∇_p_i C_i for the constraint projection step:
//
//     ∇_p_i C_i = (1/ρ_rest) Σ_j m_j ∇W(r_ij)
//
// Following Macklin 2013 PBD-Fluids: density uses Wpoly6, but the
// gradient uses Wspiky to avoid the well-known pressure-clustering
// pathology when both are Wpoly6 (Wpoly6's gradient is zero at r=0).
//
//     ∇W_spiky(r) = -(45/(πh⁶)) · (h - |r|)² · (r̂)
//     where r = p_i - p_j, r̂ = r/|r|
//
// Kernel structure: one threadgroup per active particle row, 256
// threads per row. Each thread strides through (active + static)
// neighbors, accumulating ∇W contributions into a float3 partial.
// Tree reduction in threadgroup memory combines partials, then
// scaled by mass/ρ_rest and written.
//
// Dispatch: grid (256, n_active, 1), threadgroup (256, 1, 1).
//
// Note: spiky_const is computed host-side and passed as a constant
// to avoid recomputing M_PI_F·h⁶ per thread.
// Outputs TWO things in one pass:
//   grad_C[i]      = ∇_p_i C_i  (3-vector)
//   denom_helper[i] = Σ_j |∇W(r_ij)|²  (scalar — needed for the
//                    XPBD denominator's per-neighbor term, which
//                    matters when ∇C_i is small for symmetric particles)
//
// The proper XPBD denominator for the density constraint is:
//   Σ_k |∂C_i/∂p_k|²/m_k = |∇C_i|²/m + (m/ρ₀²)·Σ_j |∇W(r_ij)|²
// The first term is what we'd compute from grad_C alone; the second
// term (per-neighbor contributions) is the missing piece that the
// constraint solver needs.
kernel void density_constraint_grad(
    device const packed_float3  *active        [[buffer(0)]],
    device const packed_float3  *static_p      [[buffer(1)]],
    device const float          *r2_aa         [[buffer(2)]],
    device const float          *r2_as         [[buffer(3)]],
    device packed_float3        *grad_C        [[buffer(4)]],
    device float                *denom_helper  [[buffer(5)]],
    constant float              &h             [[buffer(6)]],
    constant float              &spiky_const   [[buffer(7)]],
    constant float              &mass          [[buffer(8)]],
    constant float              &rho_rest      [[buffer(9)]],
    constant uint               &n_active      [[buffer(10)]],
    constant uint               &n_static      [[buffer(11)]],
    uint2 gid                                   [[thread_position_in_grid]],
    uint2 lid                                   [[thread_position_in_threadgroup]],
    uint2 tg_size                               [[threads_per_threadgroup]])
{
    uint i = gid.y;
    if (i >= n_active) return;

    threadgroup float3 partials_grad[256];
    threadgroup float  partials_denom[256];
    uint t = lid.x;
    uint T = tg_size.x;

    float h2 = h * h;
    float3 p_i = float3(active[i]);
    float3 partial_grad = float3(0.0);
    float partial_denom = 0.0;
    uint n_total = n_active + n_static;

    for (uint k = t; k < n_total; k += T) {
        float r2;
        float3 dir;
        if (k < n_active) {
            uint j = k;
            if (j == i) continue;
            r2 = r2_aa[i * n_active + j];
            if (r2 >= h2) continue;
            dir = p_i - float3(active[j]);
        } else {
            uint j = k - n_active;
            r2 = r2_as[i * n_static + j];
            if (r2 >= h2) continue;
            dir = p_i - float3(static_p[j]);
        }
        float r = sqrt(r2);
        float h_minus_r = h - r;
        // ∇W_spiky vector at this pair: spiky_const · (h-r)² · r̂ ; r̂ = dir/r
        float coef = spiky_const * h_minus_r * h_minus_r / (r + 1e-7);
        float3 grad_W = coef * dir;
        partial_grad  += mass * grad_W;
        partial_denom += grad_W.x*grad_W.x + grad_W.y*grad_W.y + grad_W.z*grad_W.z;
    }
    partials_grad[t]  = partial_grad;
    partials_denom[t] = partial_denom;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint stride = T / 2; stride > 0; stride >>= 1) {
        if (t < stride) {
            partials_grad[t]  += partials_grad[t + stride];
            partials_denom[t] += partials_denom[t + stride];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (t == 0) {
        grad_C[i] = packed_float3(partials_grad[0] / rho_rest);
        denom_helper[i] = partials_denom[0];
    }
}

// ────────────────────────────────────────────────────────────────────
// M9: density-chain backward kernels (gradient chain Option 3).
//   M6.0_bwd dist_active_static_backward
//   M6.1_bwd wpoly6_inplace_backward (in-place: ∂L/∂W → ∂L/∂r²)
//   M6.2_bwd rowsum_density_backward (broadcast)
//   M6.3_bwd dist_active_active_backward
// Each pairs with its forward kernel above.
// ────────────────────────────────────────────────────────────────────

// M6.0_bwd — Backward of dist_active_static.
// Forward: r²[i, j] = ||active[i] - static_p[j]||²
// Static positions are FROZEN, so we don't write gradients for them.
// Active gradient: ∂L/∂active[i] += Σ_j ∂L/∂r²[i,j] · 2·(active[i] - static_p[j]).
//
// One thread per active particle row — straightforward inner sum over
// static neighbors.
kernel void dist_active_static_backward(
    device const packed_float3 *active     [[buffer(0)]],
    device const packed_float3 *static_p   [[buffer(1)]],
    device const float          *grad_r2   [[buffer(2)]],
    device packed_float3        *grad_active [[buffer(3)]],   // accumulate
    constant uint              &n_active   [[buffer(4)]],
    constant uint              &n_static   [[buffer(5)]],
    uint i                                  [[thread_position_in_grid]])
{
    if (i >= n_active) return;
    float3 grad = float3(0.0);
    float3 ai = float3(active[i]);
    for (uint j = 0; j < n_static; j++) {
        float3 diff = ai - float3(static_p[j]);
        grad += 2.0 * grad_r2[i * n_static + j] * diff;
    }
    grad_active[i] = packed_float3(float3(grad_active[i]) + grad);
}

// M6.3_bwd — Backward of dist_active_active.
// Both i and j are active, so both get gradient updates. Per particle i:
//   ∂L/∂active[i] += Σ_j (grad_r2[i,j] + grad_r2[j,i]) · 2·(active[i] - active[j])
// (because r²[j,i] also depends on active[i] via the symmetric formula).
kernel void dist_active_active_backward(
    device const packed_float3 *active     [[buffer(0)]],
    device const float          *grad_r2   [[buffer(1)]],
    device packed_float3        *grad_active [[buffer(2)]],
    constant uint              &n_active   [[buffer(3)]],
    uint i                                  [[thread_position_in_grid]])
{
    if (i >= n_active) return;
    float3 grad = float3(0.0);
    float3 ai = float3(active[i]);
    for (uint j = 0; j < n_active; j++) {
        if (j == i) continue;
        float3 diff = ai - float3(active[j]);
        grad += 2.0 * (grad_r2[i * n_active + j]
                       + grad_r2[j * n_active + i]) * diff;
    }
    grad_active[i] = packed_float3(float3(grad_active[i]) + grad);
}

// M6.1_bwd — In-place backward of wpoly6_inplace.
// Forward: W = poly6_const · (h² - r²)³  if r² < h² else 0
//   dW/dr² = -3 · poly6_const · (h² - r²)²
// In-place: input buffer holds ∂L/∂W; on return holds ∂L/∂r².
// Requires saved r² from forward (passed as separate buffer).
kernel void wpoly6_inplace_backward(
    device const float *r2          [[buffer(0)]],   // saved r² from forward
    device float       *grad_W_or_r2 [[buffer(1)]],  // in: ∂L/∂W; out: ∂L/∂r²
    constant float     &h2           [[buffer(2)]],
    constant float     &poly6_const  [[buffer(3)]],
    constant uint      &n_total      [[buffer(4)]],
    uint gid                          [[thread_position_in_grid]])
{
    if (gid >= n_total) return;
    float r2_val = r2[gid];
    if (r2_val >= h2) {
        grad_W_or_r2[gid] = 0.0;
        return;
    }
    float diff = h2 - r2_val;
    float dW_dr2 = -3.0 * poly6_const * diff * diff;
    grad_W_or_r2[gid] = grad_W_or_r2[gid] * dW_dr2;
}

// M6.2_bwd — Backward of rowsum_density.
// Forward: density[i] = mass · Σ_j W[i,j]
//   ∂L/∂W[i,j] = mass · ∂L/∂density[i]   (constant across j)
// Trivial broadcast: dispatch (n_rows × n_cols), each thread writes one
// element of grad_W. (We could instead pre-multiply ∂L/∂density by mass
// host-side and skip the kernel entirely — but having an explicit
// backward kernel keeps the chain symmetric with the forward.)
kernel void rowsum_density_backward(
    device const float *grad_density [[buffer(0)]],   // [n_rows]
    device float       *grad_W       [[buffer(1)]],   // [n_rows × n_cols]
    constant float     &mass         [[buffer(2)]],
    constant uint      &n_cols       [[buffer(3)]],
    uint2 gid                         [[thread_position_in_grid]])
{
    uint i = gid.y;  // row
    uint j = gid.x;  // col
    grad_W[i * n_cols + j] = mass * grad_density[i];
}

// M9.C — density_constraint_grad_backward.
//
// Forward (M6.4) produces, per active i:
//   grad_C[i]      = (mass/ρ_rest) · Σ_k ∇W_spiky(p_i - p_k)     (3-vec)
//   denom_helper[i] = Σ_k |∇W_spiky(p_i - p_k)|²                  (scalar)
// where k ranges over active neighbors (k ≠ i) AND static neighbors.
//
// Backward: given ω_i = ∂L/∂grad_C[i] and ψ_i = ∂L/∂denom_helper[i],
// compute ∂L/∂active[i]. Static positions are frozen (no gradient).
//
// For each pair (i, k) the gradient on the per-pair grad_W is:
//   u = (mass/ρ_rest) · ω_row + 2 · ψ_row · grad_W
// where "row" is whichever particle owns the row this pair belongs to.
// The chain through ∂grad_W/∂p_i is:
//   J = c·(h-r) · [(h-r)/r · (I - ĝĝᵀ) - 2·ĝĝᵀ]  with c = spiky_const
//   u^T · J = c·(h-r) · [(h-r)/r · (u - (u·ĝ)·ĝ) - 2·(u·ĝ)·ĝ]
// And ∂grad_W/∂p_j = -J (anti-symmetric in v = p_i - p_j).
//
// Per active particle i, gradient accumulates from:
//   1. "self" rows (this i is the row): for each neighbor k (active or
//      static), contribution = +u_self^T · J(p_i, p_k)
//   2. "as-neighbor" rows (active j has i as a neighbor): for each
//      active j ≠ i with r_ij < h, contribution = -u_j^T · J(p_i, p_j)
//      (negative because ∂grad_W(p_j, p_i)/∂p_i = -J(p_i, p_j))
//
// Single thread per row i — all writes are to grad_active[i] only,
// reads are from omega/psi/positions buffers (read-only). No races.
kernel void density_constraint_grad_backward(
    device const packed_float3  *active         [[buffer(0)]],
    device const packed_float3  *static_p       [[buffer(1)]],
    device const float          *r2_aa          [[buffer(2)]],
    device const float          *r2_as          [[buffer(3)]],
    device const packed_float3  *grad_grad_C    [[buffer(4)]],   // ω
    device const float          *grad_denom_h   [[buffer(5)]],   // ψ
    device packed_float3        *grad_active    [[buffer(6)]],   // out, accumulate
    constant float              &h              [[buffer(7)]],
    constant float              &spiky_const    [[buffer(8)]],
    constant float              &mass           [[buffer(9)]],
    constant float              &rho_rest       [[buffer(10)]],
    constant uint               &n_active       [[buffer(11)]],
    constant uint               &n_static       [[buffer(12)]],
    uint i                                       [[thread_position_in_grid]])
{
    if (i >= n_active) return;
    float h2 = h * h;
    float scale = mass / rho_rest;
    float3 p_i = float3(active[i]);
    float3 omega_i = float3(grad_grad_C[i]);
    float  psi_i   = grad_denom_h[i];
    float3 grad = float3(0.0);

    // Active neighbors: handle both "self" (row i) and "as-neighbor" (row j) contributions.
    for (uint j = 0; j < n_active; j++) {
        if (j == i) continue;
        float r2 = r2_aa[i * n_active + j];
        if (r2 >= h2) continue;
        float r = sqrt(r2);
        float h_minus_r = h - r;
        float3 dir = p_i - float3(active[j]);    // v
        float3 g_hat = dir / r;
        float coef = spiky_const * h_minus_r * h_minus_r / (r + 1e-7);
        float3 grad_W = coef * dir;

        // "self" (row i):
        float3 u_self = scale * omega_i + 2.0 * psi_i * grad_W;
        float ug_s = dot(u_self, g_hat);
        float3 u_perp_s = u_self - ug_s * g_hat;
        float3 J_T_u_s = spiky_const * h_minus_r *
            (h_minus_r / r * u_perp_s - 2.0 * ug_s * g_hat);
        grad += J_T_u_s;

        // "as-neighbor" (row j, where this pair is (j, i)):
        // grad_W_for_j = grad_W(p_j, p_i) = -grad_W
        float3 omega_j = float3(grad_grad_C[j]);
        float  psi_j   = grad_denom_h[j];
        float3 u_neigh = scale * omega_j + 2.0 * psi_j * (-grad_W);
        float ug_n = dot(u_neigh, g_hat);
        float3 u_perp_n = u_neigh - ug_n * g_hat;
        float3 J_T_u_n = spiky_const * h_minus_r *
            (h_minus_r / r * u_perp_n - 2.0 * ug_n * g_hat);
        // ∂grad_W(p_j, p_i)/∂p_i = -J(p_i, p_j) → flip sign
        grad -= J_T_u_n;
    }

    // Static neighbors: only "self" contribution (no gradient flows to frozen static).
    for (uint k = 0; k < n_static; k++) {
        float r2 = r2_as[i * n_static + k];
        if (r2 >= h2) continue;
        float r = sqrt(r2);
        float h_minus_r = h - r;
        float3 dir = p_i - float3(static_p[k]);
        float3 g_hat = dir / r;
        float coef = spiky_const * h_minus_r * h_minus_r / (r + 1e-7);
        float3 grad_W = coef * dir;

        float3 u_self = scale * omega_i + 2.0 * psi_i * grad_W;
        float ug_s = dot(u_self, g_hat);
        float3 u_perp_s = u_self - ug_s * g_hat;
        float3 J_T_u_s = spiky_const * h_minus_r *
            (h_minus_r / r * u_perp_s - 2.0 * ug_s * g_hat);
        grad += J_T_u_s;
    }

    grad_active[i] = packed_float3(float3(grad_active[i]) + grad);
}

// ────────────────────────────────────────────────────────────────────
// M7: XPBD orchestration kernels
// ────────────────────────────────────────────────────────────────────

// M7.0 — Predict positions under external forces (gravity).
// Semi-implicit Euler:  v_pred = v + dt·g;  x_pred = x + dt·v_pred·s
//
// `s` is `sim_scale_inv` — the unit-system bridge that lets velocity
// stay in physical SI (m/s) while positions stay in Sibernetic's
// "particle units" (where 1 unit = simulationScale meters ≈ 7.4e-6 m).
// For toy tests where positions and velocity share a unit system,
// pass s=1.0 and behavior is unchanged from the original Euler step.
//
// pos_old and vel are read; pos_pred is written (separate buffer so
// pos_old is preserved for the velocity-recovery step at end of step).
kernel void predict_positions(
    device const packed_float3 *pos_old        [[buffer(0)]],
    device const packed_float3 *vel            [[buffer(1)]],
    device packed_float3       *pos_pred       [[buffer(2)]],
    constant float             &dt             [[buffer(3)]],
    constant float             &gravity_y      [[buffer(4)]],
    constant uint              &n              [[buffer(5)]],
    constant float             &sim_scale_inv  [[buffer(6)]],
    uint gid                                    [[thread_position_in_grid]])
{
    if (gid >= n) return;
    float3 x = float3(pos_old[gid]);
    float3 v = float3(vel[gid]);
    float3 v_pred = v + float3(0.0, gravity_y * dt, 0.0);
    pos_pred[gid] = packed_float3(x + v_pred * dt * sim_scale_inv);
}

// M7.1 — Density constraint projection (per active particle).
//
//   C_i = max(0, ρ_i/ρ_rest - 1)        (one-sided: only correct overcompression)
//   Δλ  = -(C + α/dt² · λ) / (|∇C|²/m + α/dt²)
//   Δx  = ∇C · Δλ / m
//   λ  += Δλ
//
// One-sided projection (skip when C ≤ 0) is the Macklin 2013 PBD-Fluids
// convention — fluid can "expand" freely but compression triggers
// projection. This avoids the well-known SPH pressure-clustering issue.
kernel void solve_density_constraint(
    device packed_float3       *pos_pred       [[buffer(0)]],
    device float               *lambda         [[buffer(1)]],
    device const float         *density        [[buffer(2)]],
    device const packed_float3 *grad_C         [[buffer(3)]],
    device const float         *denom_helper   [[buffer(4)]],
    constant float             &rho_rest       [[buffer(5)]],
    constant float             &mass           [[buffer(6)]],
    constant float             &alpha_inv_dt2  [[buffer(7)]],
    constant uint              &n_active       [[buffer(8)]],
    uint gid                                    [[thread_position_in_grid]])
{
    if (gid >= n_active) return;
    float C = density[gid] / rho_rest - 1.0;
    if (C <= 0.0) return;
    float3 g = float3(grad_C[gid]);
    float g2 = g.x*g.x + g.y*g.y + g.z*g.z;
    // Full XPBD denominator: |∇C_i|²/m + (m/ρ₀²)·Σ|∇W|² + α/dt²
    float helper = denom_helper[gid];
    float denom = g2 / mass + (mass / (rho_rest * rho_rest)) * helper + alpha_inv_dt2;
    if (denom < 1e-9) return;
    float dlambda = -(C + alpha_inv_dt2 * lambda[gid]) / denom;
    float3 dx = g * dlambda / mass;
    pos_pred[gid] = packed_float3(float3(pos_pred[gid]) + dx);
    lambda[gid] += dlambda;
}

// M7.B — Distance constraints (elastic bonds).
//
// For each bond (i, j) with rest length L:
//   C(p_i, p_j) = ||p_i - p_j|| - L     (signed: + stretched, - compressed)
//   ∇C_i =  (p_i - p_j)/||p_i - p_j||
//   ∇C_j = -(p_i - p_j)/||p_i - p_j||
//   For uniform mass m:  Σ|∇C|²/m = 2/m
//   Δλ = -(C + α/dt²·λ) / (2/m + α/dt²)
//   Δp_i = +∇C_i · Δλ / m
//   Δp_j = -∇C_i · Δλ / m
//
// Sequential Gauss-Seidel: one thread iterates all bonds, each bond's
// projection sees the latest particle positions. Trivially correct,
// no race conditions. For 64-particle cube with 144 bonds, this is
// ~150 ns/bond × 144 = 22 µs per iter — negligible. Parallelize via
// graph coloring later if bond count grows.
//
// bond_ij: int2 per bond (i, j), packed as int32×2.
kernel void solve_distance_constraints_seq(
    device packed_float3 *pos_pred       [[buffer(0)]],
    device float         *lambdas        [[buffer(1)]],
    device const int2    *bond_ij        [[buffer(2)]],
    device const float   *rest_len       [[buffer(3)]],
    constant float       &alpha_inv_dt2  [[buffer(4)]],
    constant float       &mass_inv       [[buffer(5)]],
    constant uint        &n_bonds        [[buffer(6)]],
    uint gid                              [[thread_position_in_grid]])
{
    if (gid != 0) return;
    for (uint b = 0; b < n_bonds; b++) {
        int i = bond_ij[b].x;
        int j = bond_ij[b].y;
        float3 pi = float3(pos_pred[i]);
        float3 pj = float3(pos_pred[j]);
        float3 dij = pi - pj;
        float d = length(dij);
        if (d < 1e-7) continue;
        float C = d - rest_len[b];
        float3 grad = dij / d;
        float denom = 2.0 * mass_inv + alpha_inv_dt2;
        float dlambda = -(C + alpha_inv_dt2 * lambdas[b]) / denom;
        pos_pred[i] = packed_float3(pi + grad * dlambda * mass_inv);
        pos_pred[j] = packed_float3(pj - grad * dlambda * mass_inv);
        lambdas[b] += dlambda;
    }
}

// M8 — solve_distance_constraints_seq_with_save: same as
// solve_distance_constraints_seq but ALSO writes per-bond pre-state
// (pi, pj, λ) to a state buffer for the backward pass.
//
// state[7*b + 0..2] = pi_pre (xyz)
// state[7*b + 3..5] = pj_pre (xyz)
// state[7*b + 6]    = λ_pre
kernel void solve_distance_constraints_seq_with_save(
    device packed_float3 *pos_pred       [[buffer(0)]],
    device float         *lambdas        [[buffer(1)]],
    device const int2    *bond_ij        [[buffer(2)]],
    device const float   *rest_len       [[buffer(3)]],
    device float         *state          [[buffer(4)]],
    constant float       &alpha_inv_dt2  [[buffer(5)]],
    constant float       &mass_inv       [[buffer(6)]],
    constant uint        &n_bonds        [[buffer(7)]],
    uint gid                              [[thread_position_in_grid]])
{
    if (gid != 0) return;
    for (uint b = 0; b < n_bonds; b++) {
        int i = bond_ij[b].x;
        int j = bond_ij[b].y;
        float3 pi = float3(pos_pred[i]);
        float3 pj = float3(pos_pred[j]);
        float lambda_pre = lambdas[b];

        // Save pre-state for the backward pass.
        uint base = b * 7;
        state[base + 0] = pi.x; state[base + 1] = pi.y; state[base + 2] = pi.z;
        state[base + 3] = pj.x; state[base + 4] = pj.y; state[base + 5] = pj.z;
        state[base + 6] = lambda_pre;

        float3 dij = pi - pj;
        float d = length(dij);
        if (d < 1e-7) continue;
        float C = d - rest_len[b];
        float3 g = dij / d;
        float D = 2.0 * mass_inv + alpha_inv_dt2;
        float dlambda = -(C + alpha_inv_dt2 * lambda_pre) / D;
        pos_pred[i] = packed_float3(pi + g * dlambda * mass_inv);
        pos_pred[j] = packed_float3(pj - g * dlambda * mass_inv);
        lambdas[b] = lambda_pre + dlambda;
    }
}

// M8 — solve_distance_constraints_seq_backward: hand-derived backward
// of one bond projection, applied sequentially in REVERSE order.
//
// For each bond (in reverse), the chain rule is:
//   δ = ω - φ                                (ω=∂L/∂pi_post, φ=∂L/∂pj_post)
//   δg = δ·g                                 (scalar)
//   δJ = (Δλ/(m·d))·(δ - δg·g) - (δg/(m·D))·g
//
//   ∂L/∂pi_pre = ω + δJ - ψ·g/D
//   ∂L/∂pj_pre = φ - δJ + ψ·g/D
//   ∂L/∂λ_pre  = -δg·α/(m·dt²·D) + ψ·(2/(m·D))
//   ∂L/∂α    += [δg/m + ψ] · (C - 2λ_pre/m)/(dt²·D²)     (per bond)
//
// where Δλ, C, g are computed from saved pre-state.
//
// pos_grad and lambda_grad are read-modify-write — sequential reverse
// order means each bond's backward sees the gradient as updated by
// LATER bonds' backward calls (consistent with forward Gauss-Seidel
// being adjoint to backward reverse-Gauss-Seidel).
kernel void solve_distance_constraints_seq_backward(
    device packed_float3 *pos_grad      [[buffer(0)]],
    device float         *lambda_grad   [[buffer(1)]],
    device float         *alpha_grad    [[buffer(2)]],   // single-element accumulator
    device const int2    *bond_ij       [[buffer(3)]],
    device const float   *rest_len      [[buffer(4)]],
    device const float   *state         [[buffer(5)]],
    constant float       &alpha_inv_dt2 [[buffer(6)]],
    constant float       &alpha_param   [[buffer(7)]],
    constant float       &dt2           [[buffer(8)]],
    constant float       &mass_inv      [[buffer(9)]],
    constant uint        &n_bonds       [[buffer(10)]],
    uint gid                             [[thread_position_in_grid]])
{
    if (gid != 0) return;

    float local_alpha_grad = 0.0;
    for (int bi = (int)n_bonds - 1; bi >= 0; bi--) {
        uint b = (uint)bi;
        int i = bond_ij[b].x;
        int j = bond_ij[b].y;
        uint base = b * 7;

        float3 pi_pre = float3(state[base+0], state[base+1], state[base+2]);
        float3 pj_pre = float3(state[base+3], state[base+4], state[base+5]);
        float lambda_pre = state[base+6];

        float3 v = pi_pre - pj_pre;
        float d = length(v);
        if (d < 1e-7) continue;
        float3 g = v / d;
        float L = rest_len[b];
        float C = d - L;
        float D = 2.0 * mass_inv + alpha_inv_dt2;
        float dlambda = -(C + alpha_inv_dt2 * lambda_pre) / D;

        float3 omega = float3(pos_grad[i]);
        float3 phi   = float3(pos_grad[j]);
        float psi    = lambda_grad[b];

        float3 delta = omega - phi;
        float delta_g = dot(delta, g);

        float coef1 = dlambda * mass_inv / d;
        float coef2 = -delta_g * mass_inv / D;
        float3 delta_J = coef1 * (delta - delta_g * g) + coef2 * g;

        pos_grad[i] = packed_float3(omega + delta_J - psi * g / D);
        pos_grad[j] = packed_float3(phi - delta_J + psi * g / D);

        lambda_grad[b] = -delta_g * alpha_param * mass_inv / (dt2 * D)
                         + psi * 2.0 * mass_inv / D;

        // Parameter gradient accumulation (single scalar α).
        float dlambda_dalpha = (C - 2.0 * lambda_pre * mass_inv) / (dt2 * D * D);
        float chain_param = delta_g * mass_inv + psi;
        local_alpha_grad += chain_param * dlambda_dalpha;
    }
    alpha_grad[0] += local_alpha_grad;
}

// M9.D — solve_density_constraint_backward.
//
// Forward:
//   if C ≤ 0: skip (one-sided projection, identity backward)
//   else:
//     C = density/ρ_rest - 1
//     g = grad_C
//     D = |g|²/m + (m/ρ²)·helper + α/dt²
//     Δλ = -(C + α/dt² · λ_pre) / D
//     pos_post = pos_pre + g·Δλ/m
//     λ_post = λ_pre + Δλ
//
// Backward (per active i, derived in DEVELOPMENT_LOG M9.D section):
//   chain = (ω·g)/m + ψ
//   ∂L/∂pos_pre   = ω                                       (identity)
//   ∂L/∂λ_pre    = -(ω·g)·α/dt²/(m·D) + ψ·(1 - α/dt²/D)
//   ∂L/∂density   = -chain/(D·ρ_rest)
//   ∂L/∂grad_C    = (Δλ/m)·ω - 2·Δλ/(m·D)·chain·g
//   ∂L/∂helper    = -Δλ·m/(ρ²·D)·chain
//   ∂L/∂ρ_rest   += chain · [density/(ρ²·D) + 2·Δλ·m·helper/(ρ³·D)]
// (∂L/∂ρ_rest is per-particle; host sums to a scalar after the kernel.)
kernel void solve_density_constraint_backward(
    device const float          *density          [[buffer(0)]],
    device const packed_float3  *grad_C           [[buffer(1)]],
    device const float          *denom_helper     [[buffer(2)]],
    device const float          *lambda_pre       [[buffer(3)]],

    device const packed_float3  *grad_pos_post    [[buffer(4)]],
    device const float          *grad_lambda_post [[buffer(5)]],

    device packed_float3        *grad_pos_pre     [[buffer(6)]],   // accumulate
    device float                *grad_lambda_pre  [[buffer(7)]],   // write
    device float                *grad_density     [[buffer(8)]],   // write
    device packed_float3        *grad_grad_C      [[buffer(9)]],   // write
    device float                *grad_denom_h     [[buffer(10)]],  // write
    device float                *grad_rho_rest    [[buffer(11)]],  // per-particle (host sums)

    constant float              &rho_rest         [[buffer(12)]],
    constant float              &mass             [[buffer(13)]],
    constant float              &alpha_inv_dt2    [[buffer(14)]],
    constant uint               &n_active         [[buffer(15)]],
    uint i                                          [[thread_position_in_grid]])
{
    if (i >= n_active) return;

    float C = density[i] / rho_rest - 1.0;
    float3 omega = float3(grad_pos_post[i]);
    float psi = grad_lambda_post[i];

    if (C <= 0.0) {
        // Forward skipped → identity backward, downstream gradients zero.
        grad_pos_pre[i] = packed_float3(float3(grad_pos_pre[i]) + omega);
        grad_lambda_pre[i] = psi;
        grad_density[i] = 0.0;
        grad_grad_C[i] = packed_float3(0.0);
        grad_denom_h[i] = 0.0;
        grad_rho_rest[i] = 0.0;
        return;
    }

    float3 g = float3(grad_C[i]);
    float g2 = g.x*g.x + g.y*g.y + g.z*g.z;
    float helper = denom_helper[i];
    float lam = lambda_pre[i];
    float D = g2 / mass + (mass / (rho_rest * rho_rest)) * helper + alpha_inv_dt2;
    float dlambda = -(C + alpha_inv_dt2 * lam) / D;

    float omega_dot_g = dot(omega, g);
    float chain = omega_dot_g / mass + psi;

    grad_pos_pre[i] = packed_float3(float3(grad_pos_pre[i]) + omega);

    grad_lambda_pre[i] = -omega_dot_g * alpha_inv_dt2 / (mass * D)
                        + psi * (1.0 - alpha_inv_dt2 / D);

    grad_density[i] = -chain / (D * rho_rest);

    float3 grad_g = (dlambda / mass) * omega
                  - (2.0 * dlambda / (mass * D)) * chain * g;
    grad_grad_C[i] = packed_float3(grad_g);

    grad_denom_h[i] = -dlambda * mass / (rho_rest * rho_rest * D) * chain;

    float dlambda_dr = density[i] / (rho_rest * rho_rest * D)
                     + 2.0 * dlambda * mass * helper
                       / (rho_rest * rho_rest * rho_rest * D);
    grad_rho_rest[i] = chain * dlambda_dr;
}

// M7.3 — Floor constraint (hard inequality).
//
//   if x.y < floor_y:  x.y = floor_y
//
// Trivial per-particle clamp. No Lagrange multiplier needed since this
// is an unilateral hard constraint with infinite stiffness.
kernel void solve_floor_constraint(
    device packed_float3 *pos_pred  [[buffer(0)]],
    constant float       &floor_y   [[buffer(1)]],
    constant uint        &n         [[buffer(2)]],
    uint gid                        [[thread_position_in_grid]])
{
    if (gid >= n) return;
    float3 x = float3(pos_pred[gid]);
    if (x.y < floor_y) {
        x.y = floor_y;
        pos_pred[gid] = packed_float3(x);
    }
}

// M7.4 — Velocity update from position change.
//
//   v_new = (x_pred_after_projection - x_old) / dt
//
// This recovers velocity from the constraint corrections, which is the
// XPBD/PBD pattern that makes constraints "physical" — particles
// gain velocity proportional to how far they were projected.
kernel void update_velocities(
    device packed_float3       *vel           [[buffer(0)]],
    device const packed_float3 *pos_old       [[buffer(1)]],
    device const packed_float3 *pos_pred      [[buffer(2)]],
    constant float             &dt            [[buffer(3)]],
    constant uint              &n             [[buffer(4)]],
    constant float             &sim_scale     [[buffer(5)]],
    uint gid                                   [[thread_position_in_grid]])
{
    if (gid >= n) return;
    // Inverse of predict's `pos += dt·v·sim_scale_inv`. dx is in particle
    // units; multiply by sim_scale to recover meters, divide by dt to
    // recover m/s. With sim_scale=1.0 this is the original v = dx/dt.
    float3 dx = float3(pos_pred[gid]) - float3(pos_old[gid]);
    vel[gid] = packed_float3(dx * sim_scale / dt);
}

// Utility — Element-wise add: dst += src.
// Used in the XPBD step to combine active-active and active-static
// densities into a single density vector before constraint projection.
kernel void add_inplace(
    device float       *dst   [[buffer(0)]],
    device const float *src   [[buffer(1)]],
    constant uint      &n     [[buffer(2)]],
    uint gid                  [[thread_position_in_grid]])
{
    if (gid >= n) return;
    dst[gid] += src[gid];
}

// ────────────────────────────────────────────────────────────────────
// M7.C: Backward kernels for the gradient chain (Option 3).
// Each pairs with its forward kernel above. Hand-derived gradients,
// no general autograd. Output gradients accumulate (use += not =) so
// multiple kernels can contribute to the same input gradient.
// ────────────────────────────────────────────────────────────────────

// M7.C-predict — Backward of predict_positions.
//
// Forward:
//   v_pred = v + dt·g     (g is a 3-vec but only g.y nonzero in our use)
//   x_pred = x + dt·v_pred = x + dt·v + dt²·g
// Inputs gradient flow:
//   dL/dx     = dL/dx_pred                  (identity)
//   dL/dv     = dt · dL/dx_pred             (chain through v_pred → x_pred)
//   dL/dg.y  += dt² · sum_i dL/dx_pred[i].y (accumulated across particles)
//
// Note: dL/dg is a scalar (just y component). For demo we accumulate it
// host-side after this kernel writes per-particle dL/dg.y contributions.
kernel void predict_positions_backward(
    device const packed_float3 *grad_x_pred  [[buffer(0)]],   // ∂L/∂x_pred
    device packed_float3       *grad_x_old   [[buffer(1)]],   // ∂L/∂x_old (output, accumulate)
    device packed_float3       *grad_vel     [[buffer(2)]],   // ∂L/∂v (output, accumulate)
    device float               *grad_g_y     [[buffer(3)]],   // per-particle dL/dg.y (output)
    constant float             &dt           [[buffer(4)]],
    constant uint              &n            [[buffer(5)]],
    uint gid                                  [[thread_position_in_grid]])
{
    if (gid >= n) return;
    float3 g_xp = float3(grad_x_pred[gid]);
    grad_x_old[gid] = packed_float3(float3(grad_x_old[gid]) + g_xp);
    grad_vel[gid]   = packed_float3(float3(grad_vel[gid])   + g_xp * dt);
    grad_g_y[gid]   = g_xp.y * dt * dt;
}

// M7.C-floor-fwd-mask — Forward floor constraint with mask emission.
// Same effect as solve_floor_constraint, but also writes clamped[i] = 1
// when the particle was clamped, 0 otherwise. The mask is consumed by
// solve_floor_constraint_backward.
kernel void solve_floor_constraint_with_mask(
    device packed_float3 *pos_pred  [[buffer(0)]],
    device int           *clamped   [[buffer(1)]],
    constant float       &floor_y   [[buffer(2)]],
    constant uint        &n         [[buffer(3)]],
    uint gid                        [[thread_position_in_grid]])
{
    if (gid >= n) return;
    float3 x = float3(pos_pred[gid]);
    if (x.y < floor_y) {
        x.y = floor_y;
        pos_pred[gid] = packed_float3(x);
        clamped[gid] = 1;
    } else {
        clamped[gid] = 0;
    }
}

// M7.C-floor-bwd — Backward of solve_floor_constraint_with_mask.
//
// Forward semantics:
//   if clamped[i]: pos_post[i] = (pos_pre.x, floor_y, pos_pre.z)
//                  pos_post.y/dpos_pre.y = 0 (clamping kills gradient)
//                  pos_post.y/dfloor_y   = 1
//   else:          pos_post = pos_pre
//
// We accumulate ∂L/∂pos_pre (caller may already have partial grads).
// ∂L/∂floor_y is per-particle and host-summed into the scalar gradient.
kernel void solve_floor_constraint_backward(
    device const packed_float3 *grad_pos_post  [[buffer(0)]],
    device packed_float3       *grad_pos_pre   [[buffer(1)]],
    device float               *grad_floor     [[buffer(2)]],
    device const int           *clamped        [[buffer(3)]],
    constant uint              &n              [[buffer(4)]],
    uint gid                                    [[thread_position_in_grid]])
{
    if (gid >= n) return;
    float3 g_post = float3(grad_pos_post[gid]);
    if (clamped[gid] != 0) {
        grad_pos_pre[gid] = packed_float3(float3(grad_pos_pre[gid]) +
                                          float3(g_post.x, 0.0, g_post.z));
        grad_floor[gid] = g_post.y;
    } else {
        grad_pos_pre[gid] = packed_float3(float3(grad_pos_pre[gid]) + g_post);
        grad_floor[gid] = 0.0;
    }
}

// M7.C-update_vel — Backward of update_velocities.
//
// Forward:
//   v_new = (x_pred - x_old) / dt
// Backward:
//   dL/dx_pred  +=  (1/dt) · dL/dv_new
//   dL/dx_old   += -(1/dt) · dL/dv_new
kernel void update_velocities_backward(
    device const packed_float3 *grad_v_new   [[buffer(0)]],   // ∂L/∂v_new
    device packed_float3       *grad_x_pred  [[buffer(1)]],   // ∂L/∂x_pred (output, accumulate)
    device packed_float3       *grad_x_old   [[buffer(2)]],   // ∂L/∂x_old (output, accumulate)
    constant float             &dt           [[buffer(3)]],
    constant uint              &n            [[buffer(4)]],
    uint gid                                  [[thread_position_in_grid]])
{
    if (gid >= n) return;
    float3 g_v = float3(grad_v_new[gid]) / dt;
    grad_x_pred[gid] = packed_float3(float3(grad_x_pred[gid]) + g_v);
    grad_x_old[gid]  = packed_float3(float3(grad_x_old[gid])  - g_v);
}
