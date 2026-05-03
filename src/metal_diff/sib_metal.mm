// sib_metal — hand-written Metal compute kernels for differentiable
// Sibernetic. M6 substrate: pairwise active×static squared distance
// matrix as the foundation for matrix-formulated SPH.
//
// Build: ./build.sh
// Run:   ./sib_metal <op> <args>
//
// Currently supported ops:
//   dist_active_static     <n_active> <n_static> <active.bin> <static.bin> <out.bin> [iters]
//   dist_active_active     <n_active> <active.bin> <out.bin> [iters]
//   wpoly6_inplace         <n_total> <h> <inout.bin> [iters]
//   rowsum_density         <n_rows> <n_cols> <mass> <W.bin> <density.bin> [iters]
//   density_constraint_grad <n_active> <n_static> <h> <mass> <rho_rest> \
//                            <active.bin> <static.bin> <r2_aa.bin> <r2_as.bin> \
//                            <gradC.bin> [iters]
//
// `iters` (optional, default 1): re-run the compute kernel `iters` times
// against the same buffers and print steady-state per-iter wall time on
// stderr. Useful for amortizing Metal startup cost (device init + shader
// compile, ~700 ms) when measuring per-step kernel performance.
//
// Binary file formats (all little-endian, fp32):
//   active.bin: [n_active * 3] floats (xyz interleaved)
//   static.bin: [n_static * 3] floats
//   out.bin:    [n_active * n_static] floats (row-major)
//
// Embedded shader source — kept inline so the binary is single-file
// and the Python driver doesn't need to track a .metallib alongside.

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *kShaderSrc = R"METAL(
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
// Semi-implicit Euler:  v_pred = v + dt·g;  x_pred = x + dt·v_pred
//
// pos_old and vel are read; pos_pred is written (separate buffer so
// pos_old is preserved for the velocity-recovery step at end of step).
kernel void predict_positions(
    device const packed_float3 *pos_old   [[buffer(0)]],
    device const packed_float3 *vel       [[buffer(1)]],
    device packed_float3       *pos_pred  [[buffer(2)]],
    constant float             &dt        [[buffer(3)]],
    constant float             &gravity_y [[buffer(4)]],
    constant uint              &n         [[buffer(5)]],
    uint gid                              [[thread_position_in_grid]])
{
    if (gid >= n) return;
    float3 x = float3(pos_old[gid]);
    float3 v = float3(vel[gid]);
    float3 v_pred = v + float3(0.0, gravity_y * dt, 0.0);
    pos_pred[gid] = packed_float3(x + v_pred * dt);
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
    device packed_float3       *vel       [[buffer(0)]],
    device const packed_float3 *pos_old   [[buffer(1)]],
    device const packed_float3 *pos_pred  [[buffer(2)]],
    constant float             &dt        [[buffer(3)]],
    constant uint              &n         [[buffer(4)]],
    uint gid                              [[thread_position_in_grid]])
{
    if (gid >= n) return;
    float3 dx = float3(pos_pred[gid]) - float3(pos_old[gid]);
    vel[gid] = packed_float3(dx / dt);
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
)METAL";

static int run_dist_active_static(int argc, char **argv) {
    if (argc != 7 && argc != 8) {
        fprintf(stderr,
            "usage: sib_metal dist_active_static "
            "<n_active> <n_static> <active.bin> <static.bin> <out.bin> [iters]\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    const char *path_active = argv[4];
    const char *path_static = argv[5];
    const char *path_out    = argv[6];
    int iters = (argc == 8) ? atoi(argv[7]) : 1;
    if (iters < 1) iters = 1;

    FILE *fa = fopen(path_active, "rb");
    FILE *fs = fopen(path_static, "rb");
    if (!fa || !fs) {
        fprintf(stderr, "cannot open input files\n");
        return 1;
    }
    float *active = (float *)malloc((size_t)n_active * 3 * sizeof(float));
    float *static_p = (float *)malloc((size_t)n_static * 3 * sizeof(float));
    if (fread(active, sizeof(float), n_active * 3, fa) != n_active * 3 ||
        fread(static_p, sizeof(float), n_static * 3, fs) != n_static * 3) {
        fprintf(stderr, "short read on input file\n");
        return 1;
    }
    fclose(fa);
    fclose(fs);

    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            fprintf(stderr, "no Metal device\n");
            return 1;
        }
        NSError *err = nil;
        id<MTLLibrary> lib = [device
            newLibraryWithSource:[NSString stringWithUTF8String:kShaderSrc]
                         options:nil
                           error:&err];
        if (!lib) {
            fprintf(stderr, "shader compile: %s\n",
                [[err localizedDescription] UTF8String]);
            return 1;
        }
        id<MTLFunction> fn = [lib newFunctionWithName:@"dist_active_static"];
        id<MTLComputePipelineState> pso =
            [device newComputePipelineStateWithFunction:fn error:&err];
        if (!pso) {
            fprintf(stderr, "pipeline: %s\n",
                [[err localizedDescription] UTF8String]);
            return 1;
        }
        id<MTLCommandQueue> queue = [device newCommandQueue];

        id<MTLBuffer> bufActive = [device
            newBufferWithBytes:active
                        length:(size_t)n_active * 3 * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufStatic = [device
            newBufferWithBytes:static_p
                        length:(size_t)n_static * 3 * sizeof(float)
                       options:MTLResourceStorageModeShared];
        size_t out_bytes = (size_t)n_active * n_static * sizeof(float);
        id<MTLBuffer> bufDist = [device
            newBufferWithLength:out_bytes
                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNa = [device
            newBufferWithBytes:&n_active
                        length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNs = [device
            newBufferWithBytes:&n_static
                        length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];

        // 16x16 threadgroup: a reasonable default for 2D grids on Apple GPUs.
        // Total threads = 256 per group, fits Apple's 1024 max comfortably.
        MTLSize threadgroup = MTLSizeMake(16, 16, 1);
        MTLSize grid = MTLSizeMake(n_active, n_static, 1);

        // Steady-state benchmark: run `iters` times, time only the
        // dispatch+wait, skip device/library/buffer setup.
        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int it = 0; it < iters; it++) {
            id<MTLCommandBuffer> cmd = [queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufActive offset:0 atIndex:0];
            [enc setBuffer:bufStatic offset:0 atIndex:1];
            [enc setBuffer:bufDist  offset:0 atIndex:2];
            [enc setBuffer:bufNa    offset:0 atIndex:3];
            [enc setBuffer:bufNs    offset:0 atIndex:4];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_iter_ms = (t1 - t0) * 1000.0 / iters;
        fprintf(stderr, "kernel: %d iters, %.3f ms/iter (steady state)\n",
                iters, per_iter_ms);

        FILE *fo = fopen(path_out, "wb");
        if (!fo) { fprintf(stderr, "cannot open output\n"); return 1; }
        fwrite([bufDist contents], 1, out_bytes, fo);
        fclose(fo);
    }
    free(active);
    free(static_p);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// Shared Metal context — created once, reused by all ops in a process.
// Keeps the per-process Metal startup (device init + shader compile,
// ~700ms) amortized across multiple op invocations if a future driver
// calls multiple ops in one binary launch.
// ──────────────────────────────────────────────────────────────────────
typedef struct {
    id<MTLDevice>       device;
    id<MTLLibrary>      lib;
    id<MTLCommandQueue> queue;
} MetalCtx;

static MetalCtx make_ctx(void) {
    MetalCtx ctx = {};
    ctx.device = MTLCreateSystemDefaultDevice();
    if (!ctx.device) { fprintf(stderr, "no Metal device\n"); exit(1); }
    NSError *err = nil;
    ctx.lib = [ctx.device
        newLibraryWithSource:[NSString stringWithUTF8String:kShaderSrc]
                     options:nil
                       error:&err];
    if (!ctx.lib) {
        fprintf(stderr, "shader compile: %s\n",
            [[err localizedDescription] UTF8String]);
        exit(1);
    }
    ctx.queue = [ctx.device newCommandQueue];
    return ctx;
}

static id<MTLComputePipelineState>
make_pso(MetalCtx ctx, const char *fn_name) {
    NSError *err = nil;
    id<MTLFunction> fn = [ctx.lib newFunctionWithName:
        [NSString stringWithUTF8String:fn_name]];
    id<MTLComputePipelineState> pso =
        [ctx.device newComputePipelineStateWithFunction:fn error:&err];
    if (!pso) {
        fprintf(stderr, "pipeline %s: %s\n", fn_name,
            [[err localizedDescription] UTF8String]);
        exit(1);
    }
    return pso;
}

// ──────────────────────────────────────────────────────────────────────
// M6.1 — wpoly6_inplace
// ──────────────────────────────────────────────────────────────────────
static int run_wpoly6_inplace(int argc, char **argv) {
    if (argc != 5 && argc != 6) {
        fprintf(stderr,
            "usage: sib_metal wpoly6_inplace "
            "<n_total> <h> <inout.bin> [iters]\n");
        return 1;
    }
    uint32_t n_total = (uint32_t)atoi(argv[2]);
    float h = (float)atof(argv[3]);
    const char *path_inout = argv[4];
    int iters = (argc == 6) ? atoi(argv[5]) : 1;
    if (iters < 1) iters = 1;

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));

    FILE *f = fopen(path_inout, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path_inout); return 1; }
    float *data = (float *)malloc((size_t)n_total * sizeof(float));
    if (fread(data, sizeof(float), n_total, f) != n_total) {
        fprintf(stderr, "short read on %s\n", path_inout); return 1;
    }
    fclose(f);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso = make_pso(ctx, "wpoly6_inplace");

        id<MTLBuffer> bufData = [ctx.device
            newBufferWithBytes:data
                        length:(size_t)n_total * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufH2 = [ctx.device
            newBufferWithBytes:&h2 length:sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufC = [ctx.device
            newBufferWithBytes:&poly6_const length:sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufN = [ctx.device
            newBufferWithBytes:&n_total length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];

        MTLSize threadgroup = MTLSizeMake(256, 1, 1);
        MTLSize grid = MTLSizeMake(n_total, 1, 1);

        // For steady-state timing of an in-place kernel, the per-iter
        // memcpy needed to restore the input would inflate the number.
        // We measure ONLY kernel dispatch+wait. After iter 1 the buffer
        // holds the output (W values, mostly zeros) — running the
        // kernel again gives wrong outputs but the same compute cost,
        // which is what we want to measure. Final on-disk output is
        // restored to the iter-1 result before write-back.
        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int it = 0; it < iters; it++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufData offset:0 atIndex:0];
            [enc setBuffer:bufH2   offset:0 atIndex:1];
            [enc setBuffer:bufC    offset:0 atIndex:2];
            [enc setBuffer:bufN    offset:0 atIndex:3];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_iter_ms = (t1 - t0) * 1000.0 / iters;
        fprintf(stderr, "kernel: %d iters, %.3f ms/iter (steady state)\n",
                iters, per_iter_ms);

        // Restore the iter-1 result for write-back: re-upload input,
        // run kernel once more.
        if (iters > 1) {
            memcpy([bufData contents], data, (size_t)n_total * sizeof(float));
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufData offset:0 atIndex:0];
            [enc setBuffer:bufH2   offset:0 atIndex:1];
            [enc setBuffer:bufC    offset:0 atIndex:2];
            [enc setBuffer:bufN    offset:0 atIndex:3];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }

        FILE *fo = fopen(path_inout, "wb");
        if (!fo) { fprintf(stderr, "cannot open output\n"); return 1; }
        fwrite([bufData contents], 1, (size_t)n_total * sizeof(float), fo);
        fclose(fo);
    }
    free(data);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M6.2 — rowsum_density
// ──────────────────────────────────────────────────────────────────────
static int run_rowsum_density(int argc, char **argv) {
    if (argc != 7 && argc != 8) {
        fprintf(stderr,
            "usage: sib_metal rowsum_density "
            "<n_rows> <n_cols> <mass> <W.bin> <density.bin> [iters]\n");
        return 1;
    }
    uint32_t n_rows = (uint32_t)atoi(argv[2]);
    uint32_t n_cols = (uint32_t)atoi(argv[3]);
    float    mass   = (float)atof(argv[4]);
    const char *path_W       = argv[5];
    const char *path_density = argv[6];
    int iters = (argc == 8) ? atoi(argv[7]) : 1;
    if (iters < 1) iters = 1;

    size_t w_count = (size_t)n_rows * n_cols;
    FILE *f = fopen(path_W, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path_W); return 1; }
    float *W = (float *)malloc(w_count * sizeof(float));
    if (fread(W, sizeof(float), w_count, f) != w_count) {
        fprintf(stderr, "short read on %s\n", path_W); return 1;
    }
    fclose(f);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso = make_pso(ctx, "rowsum_density");

        id<MTLBuffer> bufW = [ctx.device
            newBufferWithBytes:W
                        length:w_count * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufRho = [ctx.device
            newBufferWithLength:(size_t)n_rows * sizeof(float)
                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufMass = [ctx.device
            newBufferWithBytes:&mass length:sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNc = [ctx.device
            newBufferWithBytes:&n_cols length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNr = [ctx.device
            newBufferWithBytes:&n_rows length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];

        // 2D dispatch: X axis = threads within threadgroup (256),
        // Y axis = row index. Threadgroup size 256x1 means each row
        // gets its own threadgroup of 256 threads doing tree reduction.
        MTLSize threadgroup = MTLSizeMake(256, 1, 1);
        MTLSize grid = MTLSizeMake(256, n_rows, 1);

        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int it = 0; it < iters; it++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufW    offset:0 atIndex:0];
            [enc setBuffer:bufRho  offset:0 atIndex:1];
            [enc setBuffer:bufMass offset:0 atIndex:2];
            [enc setBuffer:bufNc   offset:0 atIndex:3];
            [enc setBuffer:bufNr   offset:0 atIndex:4];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_iter_ms = (t1 - t0) * 1000.0 / iters;
        fprintf(stderr, "kernel: %d iters, %.3f ms/iter (steady state)\n",
                iters, per_iter_ms);

        FILE *fo = fopen(path_density, "wb");
        if (!fo) { fprintf(stderr, "cannot open output\n"); return 1; }
        fwrite([bufRho contents], 1, (size_t)n_rows * sizeof(float), fo);
        fclose(fo);
    }
    free(W);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M6.3 — dist_active_active
// ──────────────────────────────────────────────────────────────────────
static int run_dist_active_active(int argc, char **argv) {
    if (argc != 5 && argc != 6) {
        fprintf(stderr,
            "usage: sib_metal dist_active_active "
            "<n_active> <active.bin> <out.bin> [iters]\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    const char *path_active = argv[3];
    const char *path_out    = argv[4];
    int iters = (argc == 6) ? atoi(argv[5]) : 1;
    if (iters < 1) iters = 1;

    FILE *f = fopen(path_active, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path_active); return 1; }
    float *active = (float *)malloc((size_t)n_active * 3 * sizeof(float));
    if (fread(active, sizeof(float), n_active * 3, f) != n_active * 3) {
        fprintf(stderr, "short read on %s\n", path_active); return 1;
    }
    fclose(f);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso =
            make_pso(ctx, "dist_active_active");

        id<MTLBuffer> bufActive = [ctx.device
            newBufferWithBytes:active
                        length:(size_t)n_active * 3 * sizeof(float)
                       options:MTLResourceStorageModeShared];
        size_t out_bytes = (size_t)n_active * n_active * sizeof(float);
        id<MTLBuffer> bufDist = [ctx.device
            newBufferWithLength:out_bytes
                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufN = [ctx.device
            newBufferWithBytes:&n_active length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];

        MTLSize threadgroup = MTLSizeMake(16, 16, 1);
        MTLSize grid = MTLSizeMake(n_active, n_active, 1);

        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int it = 0; it < iters; it++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufActive offset:0 atIndex:0];
            [enc setBuffer:bufDist   offset:0 atIndex:1];
            [enc setBuffer:bufN      offset:0 atIndex:2];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_iter_ms = (t1 - t0) * 1000.0 / iters;
        fprintf(stderr, "kernel: %d iters, %.3f ms/iter (steady state)\n",
                iters, per_iter_ms);

        FILE *fo = fopen(path_out, "wb");
        if (!fo) { fprintf(stderr, "cannot open output\n"); return 1; }
        fwrite([bufDist contents], 1, out_bytes, fo);
        fclose(fo);
    }
    free(active);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M6.4 — density_constraint_grad (also outputs denom_helper)
// ──────────────────────────────────────────────────────────────────────
static int run_density_constraint_grad(int argc, char **argv) {
    if (argc != 13 && argc != 14) {
        fprintf(stderr,
            "usage: sib_metal density_constraint_grad "
            "<n_active> <n_static> <h> <mass> <rho_rest> "
            "<active.bin> <static.bin> <r2_aa.bin> <r2_as.bin> "
            "<gradC.bin> <denom_helper.bin> [iters]\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h          = (float)atof(argv[4]);
    float mass       = (float)atof(argv[5]);
    float rho_rest   = (float)atof(argv[6]);
    const char *path_active = argv[7];
    const char *path_static = argv[8];
    const char *path_r2_aa  = argv[9];
    const char *path_r2_as  = argv[10];
    const char *path_gradC  = argv[11];
    const char *path_denomH = argv[12];
    int iters = (argc == 14) ? atoi(argv[13]) : 1;
    if (iters < 1) iters = 1;

    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));

    // Load all four input buffers.
    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        if (fread(buf, sizeof(float), n_floats, f) != n_floats) {
            fprintf(stderr, "short read on %s\n", path); exit(1);
        }
        fclose(f);
        return buf;
    };

    float *active   = read_floats(path_active, (size_t)n_active * 3);
    float *static_p = read_floats(path_static, (size_t)n_static * 3);
    float *r2_aa    = read_floats(path_r2_aa,  (size_t)n_active * n_active);
    float *r2_as    = read_floats(path_r2_as,  (size_t)n_active * n_static);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso =
            make_pso(ctx, "density_constraint_grad");

        id<MTLBuffer> bufActive = [ctx.device
            newBufferWithBytes:active
                        length:(size_t)n_active * 3 * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufStatic = [ctx.device
            newBufferWithBytes:static_p
                        length:(size_t)n_static * 3 * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufR2aa = [ctx.device
            newBufferWithBytes:r2_aa
                        length:(size_t)n_active * n_active * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufR2as = [ctx.device
            newBufferWithBytes:r2_as
                        length:(size_t)n_active * n_static * sizeof(float)
                       options:MTLResourceStorageModeShared];
        // packed_float3 = 12 bytes per element.
        size_t out_bytes = (size_t)n_active * 3 * sizeof(float);
        size_t denom_bytes = (size_t)n_active * sizeof(float);
        id<MTLBuffer> bufGradC = [ctx.device
            newBufferWithLength:out_bytes
                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufDenomH = [ctx.device
            newBufferWithLength:denom_bytes
                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufH        = [ctx.device newBufferWithBytes:&h
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufSpiky    = [ctx.device newBufferWithBytes:&spiky_const
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufMass     = [ctx.device newBufferWithBytes:&mass
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufRho      = [ctx.device newBufferWithBytes:&rho_rest
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNa       = [ctx.device newBufferWithBytes:&n_active
            length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNs       = [ctx.device newBufferWithBytes:&n_static
            length:sizeof(uint32_t) options:MTLResourceStorageModeShared];

        MTLSize threadgroup = MTLSizeMake(256, 1, 1);
        MTLSize grid = MTLSizeMake(256, n_active, 1);

        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int it = 0; it < iters; it++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufActive  offset:0 atIndex:0];
            [enc setBuffer:bufStatic  offset:0 atIndex:1];
            [enc setBuffer:bufR2aa    offset:0 atIndex:2];
            [enc setBuffer:bufR2as    offset:0 atIndex:3];
            [enc setBuffer:bufGradC   offset:0 atIndex:4];
            [enc setBuffer:bufDenomH  offset:0 atIndex:5];
            [enc setBuffer:bufH       offset:0 atIndex:6];
            [enc setBuffer:bufSpiky   offset:0 atIndex:7];
            [enc setBuffer:bufMass    offset:0 atIndex:8];
            [enc setBuffer:bufRho     offset:0 atIndex:9];
            [enc setBuffer:bufNa      offset:0 atIndex:10];
            [enc setBuffer:bufNs      offset:0 atIndex:11];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_iter_ms = (t1 - t0) * 1000.0 / iters;
        fprintf(stderr, "kernel: %d iters, %.3f ms/iter (steady state)\n",
                iters, per_iter_ms);

        FILE *fo = fopen(path_gradC, "wb");
        if (!fo) { fprintf(stderr, "cannot open output\n"); return 1; }
        fwrite([bufGradC contents], 1, out_bytes, fo);
        fclose(fo);
        FILE *fd = fopen(path_denomH, "wb");
        if (!fd) { fprintf(stderr, "cannot open denom output\n"); return 1; }
        fwrite([bufDenomH contents], 1, denom_bytes, fd);
        fclose(fd);
    }
    free(active); free(static_p); free(r2_aa); free(r2_as);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M7 — xpbd_step: one full XPBD timestep.
//
// Pipeline per step:
//   1. predict_positions (apply gravity to vel, integrate to x_pred)
//   2. for iter in 0..n_iters:
//        a. dist_active_active   (recompute distances after each
//        b. dist_active_static    projection — particles moved)
//        c. wpoly6_inplace × 2    (kernel-evaluate r²_aa and r²_as)
//        d. rowsum_density × 2   (compute ρ_aa and ρ_as)
//        e. add_inplace          (ρ = ρ_aa + ρ_as)
//        f. density_constraint_grad
//        g. solve_density_constraint
//        h. solve_floor_constraint
//   3. update_velocities (recover v from position change)
//
// Inputs (binary fp32 files, little-endian):
//   pos_active.bin: [n_active*3] starting positions of dynamic particles
//   vel_active.bin: [n_active*3] starting velocities
//   pos_static.bin: [n_static*3] frozen boundary positions
// Outputs:
//   pos_active_out.bin, vel_active_out.bin
//
// Lagrange multipliers (λ) are reset to zero at the start of each
// step — XPBD's "warm restart" form rather than persistent λ. This
// matches Macklin 2013 PBD-Fluids.
// ──────────────────────────────────────────────────────────────────────
static int run_xpbd_step(int argc, char **argv) {
    // Required: op + 14 base + 3 bonds = 18 args minimum
    // Optional: + bench_steps = 19
    // Pass n_bonds=0 to skip distance constraints; in that case
    // bonds.bin and alpha_dist are read but ignored.
    if (argc != 18 && argc != 19) {
        fprintf(stderr,
            "usage: sib_metal xpbd_step "
            "<n_active> <n_static> <h> <mass> <rho_rest> <dt> <gravity_y> "
            "<floor_y> <alpha_density> <n_iters> "
            "<pos_active.bin> <vel_active.bin> <pos_static.bin> "
            "<n_bonds> <bonds.bin> <alpha_dist> [bench_steps]\n"
            "       (outputs written to /tmp/xpbd_{pos,vel}_out.bin)\n"
            "       bonds.bin format: per bond, [int32 i, int32 j, "
            "float32 rest_len, float32 _pad] (16 bytes each)\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h           = (float)atof(argv[4]);
    float mass        = (float)atof(argv[5]);
    float rho_rest    = (float)atof(argv[6]);
    float dt          = (float)atof(argv[7]);
    float gravity_y   = (float)atof(argv[8]);
    float floor_y     = (float)atof(argv[9]);
    float alpha_dens  = (float)atof(argv[10]);
    uint32_t n_iters  = (uint32_t)atoi(argv[11]);
    const char *path_pos_active  = argv[12];
    const char *path_vel_active  = argv[13];
    const char *path_pos_static  = argv[14];
    uint32_t n_bonds   = (uint32_t)atoi(argv[15]);
    const char *path_bonds       = argv[16];
    float alpha_dist  = (float)atof(argv[17]);
    int bench_steps = (argc == 19) ? atoi(argv[18]) : 1;
    if (bench_steps < 1) bench_steps = 1;

    const char *out_pos_path = "/tmp/xpbd_pos_out.bin";
    const char *out_vel_path = "/tmp/xpbd_vel_out.bin";

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));
    float alpha_inv_dt2 = alpha_dens / (dt * dt);
    float alpha_dist_inv_dt2 = alpha_dist / (dt * dt);
    float mass_inv = 1.0f / mass;

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        if (fread(buf, sizeof(float), n_floats, f) != n_floats) {
            fprintf(stderr, "short read on %s\n", path); exit(1);
        }
        fclose(f);
        return buf;
    };

    float *pos_active_init = read_floats(path_pos_active, (size_t)n_active * 3);
    float *vel_active_init = read_floats(path_vel_active, (size_t)n_active * 3);
    float *pos_static      = read_floats(path_pos_static, (size_t)n_static * 3);

    // Load bonds. Format per bond: [i:int32, j:int32, rest_len:float32, pad:float32]
    // Total 16 bytes per bond. Read raw and unpack.
    void *bonds_raw = NULL;
    if (n_bonds > 0) {
        FILE *f = fopen(path_bonds, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path_bonds); exit(1); }
        bonds_raw = malloc((size_t)n_bonds * 16);
        if (fread(bonds_raw, 16, n_bonds, f) != n_bonds) {
            fprintf(stderr, "short read on %s\n", path_bonds); exit(1);
        }
        fclose(f);
    }

    @autoreleasepool {
        MetalCtx ctx = make_ctx();

        // Compile all the PSOs we need (this is the expensive Metal
        // setup; per-step only command-buffer creation/dispatch).
        id<MTLComputePipelineState> pso_dist_aa  = make_pso(ctx, "dist_active_active");
        id<MTLComputePipelineState> pso_dist_as  = make_pso(ctx, "dist_active_static");
        id<MTLComputePipelineState> pso_wpoly6   = make_pso(ctx, "wpoly6_inplace");
        id<MTLComputePipelineState> pso_density  = make_pso(ctx, "rowsum_density");
        id<MTLComputePipelineState> pso_addin    = make_pso(ctx, "add_inplace");
        id<MTLComputePipelineState> pso_grad_C   = make_pso(ctx, "density_constraint_grad");
        id<MTLComputePipelineState> pso_predict  = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_solve_d  = make_pso(ctx, "solve_density_constraint");
        id<MTLComputePipelineState> pso_solve_f  = make_pso(ctx, "solve_floor_constraint");
        id<MTLComputePipelineState> pso_updvel   = make_pso(ctx, "update_velocities");
        id<MTLComputePipelineState> pso_solve_b  = (n_bonds > 0)
            ? make_pso(ctx, "solve_distance_constraints_seq") : nil;

        size_t pos_a_bytes = (size_t)n_active * 3 * sizeof(float);
        size_t pos_s_bytes = (size_t)n_static * 3 * sizeof(float);
        size_t r2_aa_bytes = (size_t)n_active * n_active * sizeof(float);
        size_t r2_as_bytes = (size_t)n_active * n_static * sizeof(float);
        size_t dens_bytes  = (size_t)n_active * sizeof(float);

        // Persistent buffers across all steps + iters.
        id<MTLBuffer> bufPosOld   = [ctx.device newBufferWithBytes:pos_active_init
            length:pos_a_bytes options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufPosPred  = [ctx.device newBufferWithLength:pos_a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufVel      = [ctx.device newBufferWithBytes:vel_active_init
            length:pos_a_bytes options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufPosStat  = [ctx.device newBufferWithBytes:pos_static
            length:pos_s_bytes options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufR2aa     = [ctx.device newBufferWithLength:r2_aa_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufR2as     = [ctx.device newBufferWithLength:r2_as_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufWaa      = [ctx.device newBufferWithLength:r2_aa_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufWas      = [ctx.device newBufferWithLength:r2_as_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufDensAa   = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufDens     = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufGradC    = [ctx.device newBufferWithLength:pos_a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufDenomH   = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufLambda   = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];

        // Bond buffers: layout in memory is [int32 i, int32 j, float32 rest, float32 pad]
        // We split into bond_ij (int2 per bond) and rest_len (float per bond)
        // for clean kernel access. Done by repacking bonds_raw on host.
        id<MTLBuffer> bufBondIJ   = nil;
        id<MTLBuffer> bufBondRest = nil;
        id<MTLBuffer> bufBondLam  = nil;
        id<MTLBuffer> bufNbonds   = [ctx.device newBufferWithBytes:&n_bonds
            length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufAdistDt2 = [ctx.device newBufferWithBytes:&alpha_dist_inv_dt2
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufMassInv  = [ctx.device newBufferWithBytes:&mass_inv
            length:sizeof(float) options:MTLResourceStorageModeShared];
        if (n_bonds > 0) {
            int32_t *bond_ij   = (int32_t *)malloc((size_t)n_bonds * 2 * sizeof(int32_t));
            float   *bond_rest = (float *)malloc((size_t)n_bonds * sizeof(float));
            uint8_t *raw       = (uint8_t *)bonds_raw;
            for (uint32_t b = 0; b < n_bonds; b++) {
                memcpy(&bond_ij[b * 2], raw + b * 16, 8);  // i, j
                memcpy(&bond_rest[b],   raw + b * 16 + 8, 4); // rest_len
                // raw + b*16 + 12 is padding, ignored
            }
            bufBondIJ   = [ctx.device newBufferWithBytes:bond_ij
                length:(size_t)n_bonds * 2 * sizeof(int32_t)
                options:MTLResourceStorageModeShared];
            bufBondRest = [ctx.device newBufferWithBytes:bond_rest
                length:(size_t)n_bonds * sizeof(float)
                options:MTLResourceStorageModeShared];
            bufBondLam  = [ctx.device newBufferWithLength:(size_t)n_bonds * sizeof(float)
                options:MTLResourceStorageModeShared];
            free(bond_ij); free(bond_rest);
        }

        // Constants — wrapped in shared buffers so kernels can read them.
        id<MTLBuffer> bufNa  = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNs  = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufH   = [ctx.device newBufferWithBytes:&h       length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufH2  = [ctx.device newBufferWithBytes:&h2      length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufPoly6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufSpiky = [ctx.device newBufferWithBytes:&spiky_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufMass  = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufRho   = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufDt    = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufG     = [ctx.device newBufferWithBytes:&gravity_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufFloor = [ctx.device newBufferWithBytes:&floor_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufAdt2  = [ctx.device newBufferWithBytes:&alpha_inv_dt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        // n_total for wpoly6_inplace (different per call).
        uint32_t n_aa_total = n_active * n_active;
        uint32_t n_as_total = n_active * n_static;
        id<MTLBuffer> bufNaaTot = [ctx.device newBufferWithBytes:&n_aa_total length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNasTot = [ctx.device newBufferWithBytes:&n_as_total length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int step = 0; step < bench_steps; step++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];

            // ── 1. predict_positions ──
            {
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_predict];
                [enc setBuffer:bufPosOld  offset:0 atIndex:0];
                [enc setBuffer:bufVel     offset:0 atIndex:1];
                [enc setBuffer:bufPosPred offset:0 atIndex:2];
                [enc setBuffer:bufDt      offset:0 atIndex:3];
                [enc setBuffer:bufG       offset:0 atIndex:4];
                [enc setBuffer:bufNa      offset:0 atIndex:5];
                [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [enc endEncoding];
            }

            // Reset Lagrange multipliers to zero at start of step.
            // (The memset is safe here because the previous step's
            // [cmd waitUntilCompleted] already finished.)
            memset([bufLambda contents], 0, dens_bytes);
            if (n_bonds > 0) {
                memset([bufBondLam contents], 0,
                       (size_t)n_bonds * sizeof(float));
            }

            // ── 2. inner XPBD iterations ──
            for (uint32_t it = 0; it < n_iters; it++) {
                // a. dist_active_active
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_dist_aa];
                    [enc setBuffer:bufPosPred offset:0 atIndex:0];
                    [enc setBuffer:bufR2aa    offset:0 atIndex:1];
                    [enc setBuffer:bufNa      offset:0 atIndex:2];
                    [enc dispatchThreads:MTLSizeMake(n_active, n_active, 1)
                        threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
                    [enc endEncoding];
                }
                // b. dist_active_static
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_dist_as];
                    [enc setBuffer:bufPosPred offset:0 atIndex:0];
                    [enc setBuffer:bufPosStat offset:0 atIndex:1];
                    [enc setBuffer:bufR2as    offset:0 atIndex:2];
                    [enc setBuffer:bufNa      offset:0 atIndex:3];
                    [enc setBuffer:bufNs      offset:0 atIndex:4];
                    [enc dispatchThreads:MTLSizeMake(n_active, n_static, 1)
                        threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
                    [enc endEncoding];
                }
                // c. Wpoly6 on r2_aa → W_aa
                //    GPU-side blit copy first (NOT a CPU memcpy — the
                //    previous encoder hasn't executed yet, so bufR2aa
                //    on CPU is still stale. Blit encoder runs on GPU
                //    timeline, properly ordered after dist_active_active.)
                {
                    id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
                    [blit copyFromBuffer:bufR2aa sourceOffset:0
                                toBuffer:bufWaa  destinationOffset:0
                                    size:r2_aa_bytes];
                    [blit endEncoding];
                }
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_wpoly6];
                    [enc setBuffer:bufWaa    offset:0 atIndex:0];
                    [enc setBuffer:bufH2     offset:0 atIndex:1];
                    [enc setBuffer:bufPoly6  offset:0 atIndex:2];
                    [enc setBuffer:bufNaaTot offset:0 atIndex:3];
                    [enc dispatchThreads:MTLSizeMake(n_aa_total, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
                    [enc endEncoding];
                }
                // d. Wpoly6 on r2_as → W_as (same blit pattern)
                {
                    id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
                    [blit copyFromBuffer:bufR2as sourceOffset:0
                                toBuffer:bufWas  destinationOffset:0
                                    size:r2_as_bytes];
                    [blit endEncoding];
                }
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_wpoly6];
                    [enc setBuffer:bufWas    offset:0 atIndex:0];
                    [enc setBuffer:bufH2     offset:0 atIndex:1];
                    [enc setBuffer:bufPoly6  offset:0 atIndex:2];
                    [enc setBuffer:bufNasTot offset:0 atIndex:3];
                    [enc dispatchThreads:MTLSizeMake(n_as_total, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
                    [enc endEncoding];
                }
                // e. rowsum_density(W_aa, mass) → density_aa
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_density];
                    [enc setBuffer:bufWaa    offset:0 atIndex:0];
                    [enc setBuffer:bufDensAa offset:0 atIndex:1];
                    [enc setBuffer:bufMass   offset:0 atIndex:2];
                    [enc setBuffer:bufNa     offset:0 atIndex:3];  // n_cols = n_active
                    [enc setBuffer:bufNa     offset:0 atIndex:4];  // n_rows = n_active
                    [enc dispatchThreads:MTLSizeMake(256, n_active, 1)
                        threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
                    [enc endEncoding];
                }
                // f. rowsum_density(W_as, mass) → density (overall accumulator)
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_density];
                    [enc setBuffer:bufWas    offset:0 atIndex:0];
                    [enc setBuffer:bufDens   offset:0 atIndex:1];
                    [enc setBuffer:bufMass   offset:0 atIndex:2];
                    [enc setBuffer:bufNs     offset:0 atIndex:3];  // n_cols = n_static
                    [enc setBuffer:bufNa     offset:0 atIndex:4];  // n_rows = n_active
                    [enc dispatchThreads:MTLSizeMake(256, n_active, 1)
                        threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
                    [enc endEncoding];
                }
                // g. density += density_aa  (combine into one density vector)
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_addin];
                    [enc setBuffer:bufDens   offset:0 atIndex:0];
                    [enc setBuffer:bufDensAa offset:0 atIndex:1];
                    [enc setBuffer:bufNa     offset:0 atIndex:2];
                    [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [enc endEncoding];
                }
                // h. density_constraint_grad (also outputs denom_helper)
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_grad_C];
                    [enc setBuffer:bufPosPred offset:0 atIndex:0];
                    [enc setBuffer:bufPosStat offset:0 atIndex:1];
                    [enc setBuffer:bufR2aa    offset:0 atIndex:2];
                    [enc setBuffer:bufR2as    offset:0 atIndex:3];
                    [enc setBuffer:bufGradC   offset:0 atIndex:4];
                    [enc setBuffer:bufDenomH  offset:0 atIndex:5];
                    [enc setBuffer:bufH       offset:0 atIndex:6];
                    [enc setBuffer:bufSpiky   offset:0 atIndex:7];
                    [enc setBuffer:bufMass    offset:0 atIndex:8];
                    [enc setBuffer:bufRho     offset:0 atIndex:9];
                    [enc setBuffer:bufNa      offset:0 atIndex:10];
                    [enc setBuffer:bufNs      offset:0 atIndex:11];
                    [enc dispatchThreads:MTLSizeMake(256, n_active, 1)
                        threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
                    [enc endEncoding];
                }
                // i. solve_density_constraint (uses denom_helper)
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_solve_d];
                    [enc setBuffer:bufPosPred offset:0 atIndex:0];
                    [enc setBuffer:bufLambda  offset:0 atIndex:1];
                    [enc setBuffer:bufDens    offset:0 atIndex:2];
                    [enc setBuffer:bufGradC   offset:0 atIndex:3];
                    [enc setBuffer:bufDenomH  offset:0 atIndex:4];
                    [enc setBuffer:bufRho     offset:0 atIndex:5];
                    [enc setBuffer:bufMass    offset:0 atIndex:6];
                    [enc setBuffer:bufAdt2    offset:0 atIndex:7];
                    [enc setBuffer:bufNa      offset:0 atIndex:8];
                    [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [enc endEncoding];
                }
                // j. solve_distance_constraints (elastic bonds, sequential)
                if (n_bonds > 0) {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_solve_b];
                    [enc setBuffer:bufPosPred offset:0 atIndex:0];
                    [enc setBuffer:bufBondLam offset:0 atIndex:1];
                    [enc setBuffer:bufBondIJ  offset:0 atIndex:2];
                    [enc setBuffer:bufBondRest offset:0 atIndex:3];
                    [enc setBuffer:bufAdistDt2 offset:0 atIndex:4];
                    [enc setBuffer:bufMassInv  offset:0 atIndex:5];
                    [enc setBuffer:bufNbonds   offset:0 atIndex:6];
                    [enc dispatchThreads:MTLSizeMake(1, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(1, 1, 1)];
                    [enc endEncoding];
                }
                // k. solve_floor_constraint
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_solve_f];
                    [enc setBuffer:bufPosPred offset:0 atIndex:0];
                    [enc setBuffer:bufFloor   offset:0 atIndex:1];
                    [enc setBuffer:bufNa      offset:0 atIndex:2];
                    [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [enc endEncoding];
                }
            }

            // ── 3. update_velocities ──
            {
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_updvel];
                [enc setBuffer:bufVel     offset:0 atIndex:0];
                [enc setBuffer:bufPosOld  offset:0 atIndex:1];
                [enc setBuffer:bufPosPred offset:0 atIndex:2];
                [enc setBuffer:bufDt      offset:0 atIndex:3];
                [enc setBuffer:bufNa      offset:0 atIndex:4];
                [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [enc endEncoding];
            }

            [cmd commit];
            [cmd waitUntilCompleted];

            // For multi-step bench: x_new becomes x_old for next step.
            memcpy([bufPosOld contents], [bufPosPred contents], pos_a_bytes);
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_step_ms = (t1 - t0) * 1000.0 / bench_steps;
        fprintf(stderr, "xpbd_step: %d steps × %u inner iters, "
                        "%.3f ms/step (steady state)\n",
                bench_steps, n_iters, per_step_ms);

        // Write outputs (final state).
        FILE *fp = fopen(out_pos_path, "wb");
        if (!fp) { fprintf(stderr, "cannot open %s\n", out_pos_path); return 1; }
        fwrite([bufPosPred contents], 1, pos_a_bytes, fp);
        fclose(fp);
        FILE *fv = fopen(out_vel_path, "wb");
        if (!fv) { fprintf(stderr, "cannot open %s\n", out_vel_path); return 1; }
        fwrite([bufVel contents], 1, pos_a_bytes, fv);
        fclose(fv);
    }
    free(pos_active_init); free(vel_active_init); free(pos_static);
    if (bonds_raw) free(bonds_raw);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// ────────────────────────────────────────────────────────────────────
// M7.C floor — backward of solve_floor_constraint, with mask-emitting
// forward variant so the backward knows which particles were clamped.
// ────────────────────────────────────────────────────────────────────
// Forward variant of solve_floor_constraint that ALSO writes a per-
// particle "was-clamped" mask buffer. The original solve_floor_constraint
// stays unchanged so xpbd_step (M7.A+B forward-only path) doesn't need
// to allocate the mask buffer.
//
// (Defined inside the embedded shader source — see KERNELS block above.)
//
// Host runner that uses the with-mask forward + backward lives further
// down (run_step_floor_fwd / run_step_floor_bwd).
// ────────────────────────────────────────────────────────────────────

// M7.C — step_simple_fwd: ONE step of predict + update_velocities only.
// No constraints. Inputs: x_old, v_old. Outputs: x_pred, v_new.
//
// Used for validating the M7.C backward kernels via numerical-gradient
// comparison. The simplified pipeline has known closed-form derivatives
// so we can check the hand-derived backward matches finite-difference.
// ──────────────────────────────────────────────────────────────────────
static int run_step_simple_fwd(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr,
            "usage: sib_metal step_simple_fwd "
            "<n> <dt> <gravity_y> <pos_old.bin> <vel_old.bin> "
            "<pos_pred_out.bin> <vel_new_out.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    float dt = (float)atof(argv[3]);
    float g_y = (float)atof(argv[4]);

    FILE *fp = fopen(argv[5], "rb");
    FILE *fv = fopen(argv[6], "rb");
    if (!fp || !fv) { fprintf(stderr, "cannot open input\n"); return 1; }
    float *x_old = (float *)malloc((size_t)n * 3 * sizeof(float));
    float *v_old = (float *)malloc((size_t)n * 3 * sizeof(float));
    fread(x_old, sizeof(float), n * 3, fp);
    fread(v_old, sizeof(float), n * 3, fv);
    fclose(fp); fclose(fv);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_predict = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_updvel  = make_pso(ctx, "update_velocities");

        size_t nb = (size_t)n * 3 * sizeof(float);
        id<MTLBuffer> bX = [ctx.device newBufferWithBytes:x_old length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bV = [ctx.device newBufferWithBytes:v_old length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bXp = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bVn = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDt  = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bG   = [ctx.device newBufferWithBytes:&g_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN   = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // predict
        {
            id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
            [e setComputePipelineState:pso_predict];
            [e setBuffer:bX  offset:0 atIndex:0];
            [e setBuffer:bV  offset:0 atIndex:1];
            [e setBuffer:bXp offset:0 atIndex:2];
            [e setBuffer:bDt offset:0 atIndex:3];
            [e setBuffer:bG  offset:0 atIndex:4];
            [e setBuffer:bN  offset:0 atIndex:5];
            [e dispatchThreads:MTLSizeMake(n, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
            [e endEncoding];
        }
        // update_vel
        {
            id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
            [e setComputePipelineState:pso_updvel];
            [e setBuffer:bVn offset:0 atIndex:0];
            [e setBuffer:bX  offset:0 atIndex:1];
            [e setBuffer:bXp offset:0 atIndex:2];
            [e setBuffer:bDt offset:0 atIndex:3];
            [e setBuffer:bN  offset:0 atIndex:4];
            [e dispatchThreads:MTLSizeMake(n, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
            [e endEncoding];
        }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fxp = fopen(argv[7], "wb");
        fwrite([bXp contents], 1, nb, fxp); fclose(fxp);
        FILE *fvn = fopen(argv[8], "wb");
        fwrite([bVn contents], 1, nb, fvn); fclose(fvn);
    }
    free(x_old); free(v_old);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M7.C — step_simple_bwd: backward of the simplified step.
// Inputs: grad_x_pred, grad_v_new (both [n,3] fp32 arrays)
// Outputs: grad_x_old, grad_v_old (both [n,3] fp32), grad_g_y_per (n fp32)
// Host can sum grad_g_y_per to get scalar grad_g_y.
// ──────────────────────────────────────────────────────────────────────
static int run_step_simple_bwd(int argc, char **argv) {
    // op + 7 args = 9 total (argv[0]=binary, argv[1]=op,
    //                        argv[2..8] = n, dt, 2 inputs, 3 outputs).
    if (argc != 9) {
        fprintf(stderr,
            "usage: sib_metal step_simple_bwd "
            "<n> <dt> <grad_x_pred.bin> <grad_v_new.bin> "
            "<grad_x_old_out.bin> <grad_v_old_out.bin> "
            "<grad_g_y_out.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    float dt = (float)atof(argv[3]);

    FILE *fxp = fopen(argv[4], "rb");
    FILE *fvn = fopen(argv[5], "rb");
    float *g_xp = (float *)malloc((size_t)n * 3 * sizeof(float));
    float *g_vn = (float *)malloc((size_t)n * 3 * sizeof(float));
    fread(g_xp, sizeof(float), n * 3, fxp);
    fread(g_vn, sizeof(float), n * 3, fvn);
    fclose(fxp); fclose(fvn);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pbw  = make_pso(ctx, "predict_positions_backward");
        id<MTLComputePipelineState> pso_uvbw = make_pso(ctx, "update_velocities_backward");

        size_t nb = (size_t)n * 3 * sizeof(float);
        size_t ns = (size_t)n * sizeof(float);
        // Output gradients start at zero (kernels accumulate via +=).
        id<MTLBuffer> bGxp = [ctx.device newBufferWithBytes:g_xp length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGvn = [ctx.device newBufferWithBytes:g_vn length:nb
            options:MTLResourceStorageModeShared];
        // Initialize grad outputs to zero (alloc with newBufferWithLength
        // is undefined per Metal docs — explicitly memset for safety).
        id<MTLBuffer> bGxo = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGvo = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGgy = [ctx.device newBufferWithLength:ns
            options:MTLResourceStorageModeShared];
        memset([bGxo contents], 0, nb);
        memset([bGvo contents], 0, nb);
        memset([bGgy contents], 0, ns);

        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN  = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // Reverse order: backward through update_velocities first, then predict.
        {
            id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
            [e setComputePipelineState:pso_uvbw];
            [e setBuffer:bGvn offset:0 atIndex:0];
            [e setBuffer:bGxp offset:0 atIndex:1];   // accumulate into existing grad_x_pred
            [e setBuffer:bGxo offset:0 atIndex:2];   // accumulate into grad_x_old
            [e setBuffer:bDt  offset:0 atIndex:3];
            [e setBuffer:bN   offset:0 atIndex:4];
            [e dispatchThreads:MTLSizeMake(n, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
            [e endEncoding];
        }
        {
            id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
            [e setComputePipelineState:pso_pbw];
            [e setBuffer:bGxp offset:0 atIndex:0];
            [e setBuffer:bGxo offset:0 atIndex:1];
            [e setBuffer:bGvo offset:0 atIndex:2];
            [e setBuffer:bGgy offset:0 atIndex:3];
            [e setBuffer:bDt  offset:0 atIndex:4];
            [e setBuffer:bN   offset:0 atIndex:5];
            [e dispatchThreads:MTLSizeMake(n, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
            [e endEncoding];
        }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *o1 = fopen(argv[6], "wb");
        fwrite([bGxo contents], 1, nb, o1); fclose(o1);
        FILE *o2 = fopen(argv[7], "wb");
        fwrite([bGvo contents], 1, nb, o2); fclose(o2);
        FILE *o3 = fopen(argv[8], "wb");
        fwrite([bGgy contents], 1, ns, o3); fclose(o3);
    }
    free(g_xp); free(g_vn);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M7.C — step_floor_fwd: predict + floor_with_mask + update_velocities
// over K steps. Saves per-step trajectory of positions and clamp masks
// so the backward can replay through them.
//
// Argv layout:
//   [2]=n  [3]=K  [4]=dt  [5]=gravity_y  [6]=floor_y
//   [7]=pos0.bin  [8]=vel0.bin
//   [9]=traj_pos.bin (out, [K+1, n, 3] fp32 — start-of-step positions)
//   [10]=traj_clamp.bin (out, [K, n] int32 — per-step masks)
//   [11]=vel_final.bin (out, [n, 3] fp32 — velocity at end of K steps)
// argc = 12 (op + 10 positional args).
// ──────────────────────────────────────────────────────────────────────
static int run_step_floor_fwd(int argc, char **argv) {
    if (argc != 12) {
        fprintf(stderr,
            "usage: sib_metal step_floor_fwd "
            "<n> <K> <dt> <gravity_y> <floor_y> "
            "<pos0.bin> <vel0.bin> "
            "<traj_pos.bin> <traj_clamp.bin> <vel_final.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    uint32_t K = (uint32_t)atoi(argv[3]);
    float dt   = (float)atof(argv[4]);
    float g_y  = (float)atof(argv[5]);
    float fl_y = (float)atof(argv[6]);
    const char *path_x0    = argv[7];
    const char *path_v0    = argv[8];
    const char *path_tpos  = argv[9];
    const char *path_tclmp = argv[10];
    const char *path_vfinal= argv[11];

    size_t nb = (size_t)n * 3 * sizeof(float);
    size_t nm = (size_t)n * sizeof(int32_t);
    FILE *fp = fopen(path_x0, "rb");
    FILE *fv = fopen(path_v0, "rb");
    if (!fp || !fv) { fprintf(stderr, "cannot open input\n"); return 1; }
    float *x0 = (float *)malloc(nb);
    float *v0 = (float *)malloc(nb);
    fread(x0, 1, nb, fp); fread(v0, 1, nb, fv);
    fclose(fp); fclose(fv);

    float *traj_pos = (float *)malloc((size_t)(K + 1) * nb);
    int32_t *traj_clamp = (int32_t *)malloc((size_t)K * nm);
    memcpy(traj_pos, x0, nb);  // frame 0 = initial positions

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pred  = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_floor = make_pso(ctx, "solve_floor_constraint_with_mask");
        id<MTLComputePipelineState> pso_uv    = make_pso(ctx, "update_velocities");

        id<MTLBuffer> bX  = [ctx.device newBufferWithBytes:x0 length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bV  = [ctx.device newBufferWithBytes:v0 length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bXp = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bClmp = [ctx.device newBufferWithLength:nm
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGy = [ctx.device newBufferWithBytes:&g_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bFy = [ctx.device newBufferWithBytes:&fl_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN  = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        for (uint32_t k = 0; k < K; k++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            // predict
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pred];
              [e setBuffer:bX offset:0 atIndex:0]; [e setBuffer:bV offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bGy offset:0 atIndex:4]; [e setBuffer:bN offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n,1,1)
                  threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // floor with mask
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_floor];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bClmp offset:0 atIndex:1];
              [e setBuffer:bFy offset:0 atIndex:2]; [e setBuffer:bN offset:0 atIndex:3];
              [e dispatchThreads:MTLSizeMake(n,1,1)
                  threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // update_vel
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uv];
              [e setBuffer:bV offset:0 atIndex:0]; [e setBuffer:bX offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bN offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n,1,1)
                  threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // Save this step's mask + post-floor positions, advance state.
            memcpy(&traj_clamp[k * n], [bClmp contents], nm);
            memcpy((char *)traj_pos + (size_t)(k + 1) * nb, [bXp contents], nb);
            memcpy([bX contents], [bXp contents], nb);
        }

        FILE *fvf = fopen(path_vfinal, "wb");
        if (!fvf) { fprintf(stderr, "cannot open vel_final out\n"); return 1; }
        fwrite([bV contents], 1, nb, fvf); fclose(fvf);
    }
    FILE *ftp = fopen(path_tpos, "wb");
    fwrite(traj_pos, 1, (size_t)(K + 1) * nb, ftp); fclose(ftp);
    FILE *ftc = fopen(path_tclmp, "wb");
    fwrite(traj_clamp, 1, (size_t)K * nm, ftc); fclose(ftc);

    free(x0); free(v0); free(traj_pos); free(traj_clamp);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M7.C — step_floor_bwd: backward through K steps of (predict + floor +
// update_velocities), given trajectory and masks from forward.
//
// Inputs:  ∂L/∂x_final (gradient on final positions, [n, 3] fp32)
// Outputs: ∂L/∂x_init, ∂L/∂v_init, ∂L/∂g_y (scalar), ∂L/∂floor_y (scalar)
//
// Argv layout:
//   [2]=n  [3]=K  [4]=dt
//   [5]=traj_pos.bin  [6]=traj_clamp.bin
//   [7]=grad_x_final.bin
//   [8]=grad_x_init.bin (out)  [9]=grad_v_init.bin (out)
//   [10]=grad_g_y.bin (out, single fp32)  [11]=grad_floor_y.bin (out, single fp32)
// argc = 12.
// ──────────────────────────────────────────────────────────────────────
static int run_step_floor_bwd(int argc, char **argv) {
    if (argc != 12) {
        fprintf(stderr,
            "usage: sib_metal step_floor_bwd "
            "<n> <K> <dt> <traj_pos.bin> <traj_clamp.bin> "
            "<grad_x_final.bin> "
            "<grad_x_init.bin> <grad_v_init.bin> "
            "<grad_g_y.bin> <grad_floor_y.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    uint32_t K = (uint32_t)atoi(argv[3]);
    float dt = (float)atof(argv[4]);

    size_t nb = (size_t)n * 3 * sizeof(float);
    size_t nm = (size_t)n * sizeof(int32_t);
    size_t ns = (size_t)n * sizeof(float);

    FILE *ftp = fopen(argv[5], "rb");
    FILE *ftc = fopen(argv[6], "rb");
    FILE *fgx = fopen(argv[7], "rb");
    if (!ftp || !ftc || !fgx) { fprintf(stderr, "cannot open input\n"); return 1; }
    float   *traj_pos   = (float *)malloc((size_t)(K + 1) * nb);
    int32_t *traj_clamp = (int32_t *)malloc((size_t)K * nm);
    float   *grad_x_fin = (float *)malloc(nb);
    fread(traj_pos,   1, (size_t)(K + 1) * nb, ftp);
    fread(traj_clamp, 1, (size_t)K * nm, ftc);
    fread(grad_x_fin, 1, nb, fgx);
    fclose(ftp); fclose(ftc); fclose(fgx);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pbw   = make_pso(ctx, "predict_positions_backward");
        id<MTLComputePipelineState> pso_uvbw  = make_pso(ctx, "update_velocities_backward");
        id<MTLComputePipelineState> pso_flbw  = make_pso(ctx, "solve_floor_constraint_backward");

        // Persistent gradient buffers (accumulate across timesteps).
        id<MTLBuffer> bGx_old = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];  // ∂L/∂x_old (running)
        id<MTLBuffer> bGv_old = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];  // ∂L/∂v_old (running)
        // Per-step working buffers.
        id<MTLBuffer> bGx_pred = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGv_new  = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGx_pre_floor = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        // Per-particle parameter gradients (summed host-side).
        id<MTLBuffer> bGgy_per = [ctx.device newBufferWithLength:ns
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGfy_per = [ctx.device newBufferWithLength:ns
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bClmp = [ctx.device newBufferWithLength:nm
            options:MTLResourceStorageModeShared];

        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN  = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        // Initialize: at the END of the simulation, gradient on the final
        // post-floor position is grad_x_fin. v_new gradient = 0 (loss only
        // depends on final positions in this demo). x_old grad = 0 too.
        memcpy([bGx_old contents], grad_x_fin, nb);  // initial x_old grad = grad_x_fin
        memset([bGv_old contents], 0, nb);
        // Accumulators for parameter gradients.
        float total_g_y = 0.0f;
        float total_fy = 0.0f;

        // Walk backwards through K steps.
        for (int32_t k = (int32_t)K - 1; k >= 0; k--) {
            // At this step's BACKWARD entry: bGx_old holds ∂L/∂x_post_floor
            // for step k (i.e., ∂L/∂position at end of step k = start of step k+1).
            // We need to chain through: update_vel ← floor ← predict.

            // Reset per-step working grads.
            memset([bGx_pred contents], 0, nb);
            memset([bGv_new contents], 0, nb);
            memcpy([bClmp contents], &traj_clamp[k * n], nm);

            // Note: in the forward, after this step's update_vel:
            //   v_new[k] is the "vel" at end of step (input to next step's predict)
            //   But our loss only depends on FINAL positions (step K). So
            //   ∂L/∂v_new[K-1] = 0. For interior steps, the gradient flows
            //   from next step's ∂L/∂v_old (which we'd accumulate into here).
            //   In our simplified driver: we maintain bGv_old as the running
            //   "what was ∂L/∂v at next-step-start" i.e. ∂L/∂v_new[k].
            // So: ∂L/∂v_new[k] = bGv_old (which is the running v gradient from
            //                              prior backward steps, or 0 at start).

            // bGx_old at this point = ∂L/∂x_post_floor[k] from prior steps.
            // We need to also incorporate any ∂L/∂x_old that came from this
            // step's update_velocities forward (which used x_old[k] as input).
            // → handled by update_velocities_backward.

            // ── Backward through update_velocities ──
            // Forward: v_new = (x_post_floor - x_old) / dt
            // Inputs to grad: ∂L/∂v_new (in bGv_old), need to flow into
            //   ∂L/∂x_post_floor (add to bGx_old) and ∂L/∂x_old (output).
            // Note: bGx_old accumulates the "x_post_floor" gradient from
            // both the next step's x_old grad AND this step's update_vel
            // contribution. That's exactly what we want.
            // After update_vel backward, we will use bGx_old as
            // ∂L/∂x_post_floor and accumulate ∂L/∂x_old separately into
            // a fresh buffer (then promote it to bGx_old at end of step).
            id<MTLBuffer> bGx_old_new = [ctx.device newBufferWithLength:nb
                options:MTLResourceStorageModeShared];
            memset([bGx_old_new contents], 0, nb);

            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uvbw];
              [e setBuffer:bGv_old   offset:0 atIndex:0];   // ∂L/∂v_new
              [e setBuffer:bGx_old   offset:0 atIndex:1];   // accumulate into ∂L/∂x_post_floor
              [e setBuffer:bGx_old_new offset:0 atIndex:2]; // accumulate into ∂L/∂x_old (this step)
              [e setBuffer:bDt offset:0 atIndex:3]; [e setBuffer:bN offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }

            // ── Backward through floor ──
            // Inputs: ∂L/∂x_post_floor (= bGx_old). Outputs: ∂L/∂x_pre_floor
            // (which we'll call bGx_pred since that's pos before floor).
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_flbw];
              [e setBuffer:bGx_old offset:0 atIndex:0];   // ∂L/∂x_post_floor
              [e setBuffer:bGx_pred offset:0 atIndex:1];  // accumulate into ∂L/∂x_pre_floor (= ∂L/∂x_pred)
              [e setBuffer:bGfy_per offset:0 atIndex:2];
              [e setBuffer:bClmp offset:0 atIndex:3]; [e setBuffer:bN offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }

            // ── Backward through predict ──
            // Inputs: ∂L/∂x_pred (bGx_pred). Outputs: ∂L/∂x_old (accumulate
            // into bGx_old_new), ∂L/∂v (accumulate into a fresh bGv_old_new),
            // ∂L/∂g_y (per-particle into bGgy_per).
            id<MTLBuffer> bGv_old_new = [ctx.device newBufferWithLength:nb
                options:MTLResourceStorageModeShared];
            memset([bGv_old_new contents], 0, nb);
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pbw];
              [e setBuffer:bGx_pred offset:0 atIndex:0];
              [e setBuffer:bGx_old_new offset:0 atIndex:1];   // accumulate
              [e setBuffer:bGv_old_new offset:0 atIndex:2];   // accumulate
              [e setBuffer:bGgy_per offset:0 atIndex:3];
              [e setBuffer:bDt offset:0 atIndex:4]; [e setBuffer:bN offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }

            [cmd commit]; [cmd waitUntilCompleted];

            // Sum per-particle parameter grads into running totals.
            float *gy = (float *)[bGgy_per contents];
            float *fy = (float *)[bGfy_per contents];
            for (uint32_t i = 0; i < n; i++) { total_g_y += gy[i]; total_fy += fy[i]; }

            // Promote per-step grads to running grads for the next (earlier) step.
            memcpy([bGx_old contents], [bGx_old_new contents], nb);
            memcpy([bGv_old contents], [bGv_old_new contents], nb);
        }

        // Write outputs.
        FILE *o1 = fopen(argv[8], "wb");
        fwrite([bGx_old contents], 1, nb, o1); fclose(o1);
        FILE *o2 = fopen(argv[9], "wb");
        fwrite([bGv_old contents], 1, nb, o2); fclose(o2);
        FILE *o3 = fopen(argv[10], "wb");
        fwrite(&total_g_y, 1, sizeof(float), o3); fclose(o3);
        FILE *o4 = fopen(argv[11], "wb");
        fwrite(&total_fy, 1, sizeof(float), o4); fclose(o4);
    }
    free(traj_pos); free(traj_clamp); free(grad_x_fin);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M8 — step_bond_fwd: K-step forward of (predict + bonds_with_save +
// update_velocities). NO floor for simplicity — focus on isolating the
// bond physics so backward correctness is unambiguous. n_iters fixed at 1.
//
// Saves:
//   traj_pos.bin    — [K+1, n, 3] fp32 positions at start of each step
//   traj_state.bin  — [K, n_bonds*7] fp32 per-bond pre-state per step
//   traj_lambda.bin — [K, n_bonds] fp32 lambda before each step's projection
// ──────────────────────────────────────────────────────────────────────
static int run_step_bond_fwd(int argc, char **argv) {
    if (argc != 14) {
        fprintf(stderr,
            "usage: sib_metal step_bond_fwd "
            "<n> <K> <n_bonds> <dt> <gravity_y> <mass> <alpha_dist> "
            "<pos0.bin> <vel0.bin> <bonds.bin> "
            "<traj_pos.bin> <traj_state.bin> <vel_final.bin>\n");
        return 1;
    }
    uint32_t n        = (uint32_t)atoi(argv[2]);
    uint32_t K        = (uint32_t)atoi(argv[3]);
    uint32_t n_bonds  = (uint32_t)atoi(argv[4]);
    float dt          = (float)atof(argv[5]);
    float g_y         = (float)atof(argv[6]);
    float mass        = (float)atof(argv[7]);
    float alpha_dist  = (float)atof(argv[8]);
    const char *p_x0    = argv[9];
    const char *p_v0    = argv[10];
    const char *p_bonds = argv[11];
    const char *p_tpos  = argv[12];
    const char *p_tstate = argv[13];
    const char *p_vfin  = "/tmp/sibm_step_bond_vfinal.bin";  // fixed for now

    float alpha_inv_dt2 = alpha_dist / (dt * dt);
    float mass_inv = 1.0f / mass;

    size_t nb = (size_t)n * 3 * sizeof(float);
    size_t state_bytes_per_step = (size_t)n_bonds * 7 * sizeof(float);

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        if (fread(buf, sizeof(float), n_floats, f) != n_floats) {
            fprintf(stderr, "short read on %s\n", path); exit(1);
        }
        fclose(f); return buf;
    };
    float *x0 = read_floats(p_x0, (size_t)n * 3);
    float *v0 = read_floats(p_v0, (size_t)n * 3);

    // Read bonds and unpack to (i,j) and rest_len arrays.
    FILE *fb = fopen(p_bonds, "rb");
    if (!fb) { fprintf(stderr, "cannot open %s\n", p_bonds); return 1; }
    void *bonds_raw = malloc((size_t)n_bonds * 16);
    fread(bonds_raw, 16, n_bonds, fb); fclose(fb);
    int32_t *bond_ij = (int32_t *)malloc((size_t)n_bonds * 2 * sizeof(int32_t));
    float *bond_rest = (float *)malloc((size_t)n_bonds * sizeof(float));
    for (uint32_t b = 0; b < n_bonds; b++) {
        memcpy(&bond_ij[b * 2], (uint8_t *)bonds_raw + b * 16, 8);
        memcpy(&bond_rest[b],   (uint8_t *)bonds_raw + b * 16 + 8, 4);
    }
    free(bonds_raw);

    float *traj_pos = (float *)malloc((size_t)(K + 1) * nb);
    float *traj_state = (float *)malloc((size_t)K * state_bytes_per_step);
    memcpy(traj_pos, x0, nb);  // frame 0

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pred = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_bond = make_pso(ctx, "solve_distance_constraints_seq_with_save");
        id<MTLComputePipelineState> pso_uv   = make_pso(ctx, "update_velocities");

        id<MTLBuffer> bX = [ctx.device newBufferWithBytes:x0 length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bV = [ctx.device newBufferWithBytes:v0 length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bXp = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bBondIJ = [ctx.device newBufferWithBytes:bond_ij
            length:(size_t)n_bonds * 2 * sizeof(int32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bBondRest = [ctx.device newBufferWithBytes:bond_rest
            length:(size_t)n_bonds * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bLambda = [ctx.device newBufferWithLength:(size_t)n_bonds * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bState = [ctx.device newBufferWithLength:state_bytes_per_step
            options:MTLResourceStorageModeShared];

        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGy = [ctx.device newBufferWithBytes:&g_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bAdt2 = [ctx.device newBufferWithBytes:&alpha_inv_dt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bMinv = [ctx.device newBufferWithBytes:&mass_inv length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNb = [ctx.device newBufferWithBytes:&n_bonds length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        for (uint32_t k = 0; k < K; k++) {
            // Reset lambdas to zero at start of each step.
            memset([bLambda contents], 0, (size_t)n_bonds * sizeof(float));

            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            // predict
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pred];
              [e setBuffer:bX  offset:0 atIndex:0]; [e setBuffer:bV  offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bGy offset:0 atIndex:4]; [e setBuffer:bN  offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n,1,1)
                  threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // bonds (n_iters=1)
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_bond];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bLambda offset:0 atIndex:1];
              [e setBuffer:bBondIJ offset:0 atIndex:2]; [e setBuffer:bBondRest offset:0 atIndex:3];
              [e setBuffer:bState offset:0 atIndex:4]; [e setBuffer:bAdt2 offset:0 atIndex:5];
              [e setBuffer:bMinv offset:0 atIndex:6]; [e setBuffer:bNb offset:0 atIndex:7];
              [e dispatchThreads:MTLSizeMake(1,1,1)
                  threadsPerThreadgroup:MTLSizeMake(1,1,1)];
              [e endEncoding]; }
            // update_vel
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uv];
              [e setBuffer:bV offset:0 atIndex:0]; [e setBuffer:bX offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bN offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n,1,1)
                  threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // Save state for this step, advance positions.
            memcpy((char *)traj_state + (size_t)k * state_bytes_per_step,
                   [bState contents], state_bytes_per_step);
            memcpy((char *)traj_pos + (size_t)(k + 1) * nb,
                   [bXp contents], nb);
            memcpy([bX contents], [bXp contents], nb);
        }

        FILE *fvf = fopen(p_vfin, "wb");
        fwrite([bV contents], 1, nb, fvf); fclose(fvf);
    }
    FILE *ftp = fopen(p_tpos, "wb");
    fwrite(traj_pos, 1, (size_t)(K + 1) * nb, ftp); fclose(ftp);
    FILE *fts = fopen(p_tstate, "wb");
    fwrite(traj_state, 1, (size_t)K * state_bytes_per_step, fts); fclose(fts);

    free(x0); free(v0); free(bond_ij); free(bond_rest);
    free(traj_pos); free(traj_state);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M8 — step_bond_bwd: backward through K-step bond simulation. Walks
// trajectory in REVERSE, calling per-step backward kernels in order:
// update_vel_backward → bond_backward → predict_backward.
//
// Inputs:
//   grad_x_final.bin — [n, 3] fp32 ∂L/∂(positions at end)
//   plus traj_pos, traj_state, bonds from forward
//
// Outputs:
//   grad_x_init, grad_v_init, grad_alpha (single-element fp32)
// ──────────────────────────────────────────────────────────────────────
static int run_step_bond_bwd(int argc, char **argv) {
    if (argc != 16) {
        fprintf(stderr,
            "usage: sib_metal step_bond_bwd "
            "<n> <K> <n_bonds> <dt> <mass> <alpha_dist> "
            "<bonds.bin> <traj_pos.bin> <traj_state.bin> <grad_x_final.bin> "
            "<grad_x_init.bin> <grad_v_init.bin> <grad_alpha.bin>\n");
        return 1;
    }
    uint32_t n       = (uint32_t)atoi(argv[2]);
    uint32_t K       = (uint32_t)atoi(argv[3]);
    uint32_t n_bonds = (uint32_t)atoi(argv[4]);
    float dt         = (float)atof(argv[5]);
    float mass       = (float)atof(argv[6]);
    float alpha_dist = (float)atof(argv[7]);
    const char *p_bonds  = argv[8];
    const char *p_tpos   = argv[9];
    const char *p_tstate = argv[10];
    const char *p_gxf    = argv[11];
    const char *p_gxi    = argv[12];
    const char *p_gvi    = argv[13];
    const char *p_galpha = argv[14];
    // argv[15] currently unused — placeholder for grad_L_rest if added.

    float alpha_inv_dt2 = alpha_dist / (dt * dt);
    float dt2 = dt * dt;
    float mass_inv = 1.0f / mass;

    size_t nb = (size_t)n * 3 * sizeof(float);
    size_t state_bytes_per_step = (size_t)n_bonds * 7 * sizeof(float);

    // Read bonds → (ij, rest)
    FILE *fb = fopen(p_bonds, "rb");
    void *bonds_raw = malloc((size_t)n_bonds * 16);
    fread(bonds_raw, 16, n_bonds, fb); fclose(fb);
    int32_t *bond_ij = (int32_t *)malloc((size_t)n_bonds * 2 * sizeof(int32_t));
    float *bond_rest = (float *)malloc((size_t)n_bonds * sizeof(float));
    for (uint32_t b = 0; b < n_bonds; b++) {
        memcpy(&bond_ij[b * 2], (uint8_t *)bonds_raw + b * 16, 8);
        memcpy(&bond_rest[b],   (uint8_t *)bonds_raw + b * 16 + 8, 4);
    }
    free(bonds_raw);

    float *traj_state = (float *)malloc((size_t)K * state_bytes_per_step);
    float *grad_x_fin = (float *)malloc(nb);
    FILE *fts = fopen(p_tstate, "rb");
    FILE *fgx = fopen(p_gxf, "rb");
    fread(traj_state, 1, (size_t)K * state_bytes_per_step, fts);
    fread(grad_x_fin, 1, nb, fgx);
    fclose(fts); fclose(fgx);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pbw  = make_pso(ctx, "predict_positions_backward");
        id<MTLComputePipelineState> pso_uvbw = make_pso(ctx, "update_velocities_backward");
        id<MTLComputePipelineState> pso_bbw  = make_pso(ctx, "solve_distance_constraints_seq_backward");

        id<MTLBuffer> bGx_old = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGv_old = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGx_pred = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGlam = [ctx.device newBufferWithLength:(size_t)n_bonds * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGalpha = [ctx.device newBufferWithLength:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGgy_per = [ctx.device newBufferWithLength:(size_t)n * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bState = [ctx.device newBufferWithLength:state_bytes_per_step
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bBondIJ = [ctx.device newBufferWithBytes:bond_ij
            length:(size_t)n_bonds * 2 * sizeof(int32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bBondRest = [ctx.device newBufferWithBytes:bond_rest
            length:(size_t)n_bonds * sizeof(float)
            options:MTLResourceStorageModeShared];

        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDt2 = [ctx.device newBufferWithBytes:&dt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bAdt2 = [ctx.device newBufferWithBytes:&alpha_inv_dt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bAlpha = [ctx.device newBufferWithBytes:&alpha_dist length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bMinv = [ctx.device newBufferWithBytes:&mass_inv length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNb = [ctx.device newBufferWithBytes:&n_bonds length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        // Initial: ∂L/∂x_final → ∂L/∂x_old (running across reverse steps).
        memcpy([bGx_old contents], grad_x_fin, nb);
        memset([bGv_old contents], 0, nb);
        memset([bGalpha contents], 0, sizeof(float));

        // Walk backward through K steps.
        for (int32_t k = (int32_t)K - 1; k >= 0; k--) {
            // Per-step lambda gradient resets at start (forward resets λ to 0).
            memset([bGlam contents], 0, (size_t)n_bonds * sizeof(float));
            // Per-step working buffers.
            memset([bGx_pred contents], 0, nb);

            // Load this step's saved bond state.
            memcpy([bState contents],
                   (char *)traj_state + (size_t)k * state_bytes_per_step,
                   state_bytes_per_step);

            // bGx_old_new accumulates ∂L/∂x_old for THIS step (separate
            // buffer so we don't mix with the running bGx_old that holds
            // the post-step gradient).
            id<MTLBuffer> bGx_old_new = [ctx.device newBufferWithLength:nb
                options:MTLResourceStorageModeShared];
            id<MTLBuffer> bGv_old_new = [ctx.device newBufferWithLength:nb
                options:MTLResourceStorageModeShared];
            memset([bGx_old_new contents], 0, nb);
            memset([bGv_old_new contents], 0, nb);

            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            // (1) update_vel backward: ∂L/∂v_new (=bGv_old running) flows into
            //     ∂L/∂x_post_bonds (accumulate into bGx_old) and ∂L/∂x_old
            //     (accumulate into bGx_old_new).
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uvbw];
              [e setBuffer:bGv_old offset:0 atIndex:0];
              [e setBuffer:bGx_old offset:0 atIndex:1];
              [e setBuffer:bGx_old_new offset:0 atIndex:2];
              [e setBuffer:bDt offset:0 atIndex:3]; [e setBuffer:bN offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // bGx_old now = ∂L/∂x_post_bonds for this step. The bond_backward
            // overwrites it in place (read pos_grad[i,j], write pos_grad[i,j]).
            cmd = [ctx.queue commandBuffer];
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_bbw];
              [e setBuffer:bGx_old offset:0 atIndex:0];   // pos_grad in/out
              [e setBuffer:bGlam offset:0 atIndex:1];     // lambda_grad in/out
              [e setBuffer:bGalpha offset:0 atIndex:2];   // alpha grad accum
              [e setBuffer:bBondIJ offset:0 atIndex:3];
              [e setBuffer:bBondRest offset:0 atIndex:4];
              [e setBuffer:bState offset:0 atIndex:5];
              [e setBuffer:bAdt2 offset:0 atIndex:6];
              [e setBuffer:bAlpha offset:0 atIndex:7];
              [e setBuffer:bDt2 offset:0 atIndex:8];
              [e setBuffer:bMinv offset:0 atIndex:9];
              [e setBuffer:bNb offset:0 atIndex:10];
              [e dispatchThreads:MTLSizeMake(1,1,1) threadsPerThreadgroup:MTLSizeMake(1,1,1)];
              [e endEncoding]; }
            // bGx_old now = ∂L/∂x_pre_bonds = ∂L/∂x_pred (post-predict).
            // Move into bGx_pred for predict_backward.
            // (We could just pass bGx_old to predict_backward directly,
            // but keep separate for clarity.)
            { id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
              [blit copyFromBuffer:bGx_old sourceOffset:0
                          toBuffer:bGx_pred destinationOffset:0
                              size:nb];
              [blit endEncoding]; }
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pbw];
              [e setBuffer:bGx_pred offset:0 atIndex:0];
              [e setBuffer:bGx_old_new offset:0 atIndex:1];   // accumulate
              [e setBuffer:bGv_old_new offset:0 atIndex:2];   // accumulate
              [e setBuffer:bGgy_per offset:0 atIndex:3];
              [e setBuffer:bDt offset:0 atIndex:4]; [e setBuffer:bN offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // Promote per-step gradients to running gradients for previous step.
            memcpy([bGx_old contents], [bGx_old_new contents], nb);
            memcpy([bGv_old contents], [bGv_old_new contents], nb);
        }

        FILE *o1 = fopen(p_gxi, "wb");
        fwrite([bGx_old contents], 1, nb, o1); fclose(o1);
        FILE *o2 = fopen(p_gvi, "wb");
        fwrite([bGv_old contents], 1, nb, o2); fclose(o2);
        FILE *o3 = fopen(p_galpha, "wb");
        fwrite([bGalpha contents], 1, sizeof(float), o3); fclose(o3);
    }
    free(bond_ij); free(bond_rest); free(traj_state); free(grad_x_fin);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.A — density_as_fwd: compute density at each active particle from
// active×static neighbor SPH kernel evaluation. No active-active term;
// kept minimal to validate the density backward chain in isolation.
//
// Pipeline: dist_active_static → wpoly6_inplace → rowsum_density.
// Saves r² (post-distance, pre-Wpoly6) for backward.
// ──────────────────────────────────────────────────────────────────────
static int run_density_as_fwd(int argc, char **argv) {
    if (argc != 10) {
        fprintf(stderr,
            "usage: sib_metal density_as_fwd "
            "<n_active> <n_static> <h> <mass> "
            "<active.bin> <static.bin> <density_out.bin> <r2_saved_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h    = (float)atof(argv[4]);
    float mass = (float)atof(argv[5]);

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    uint32_t n_r2 = n_active * n_static;

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        if (fread(buf, sizeof(float), n_floats, f) != n_floats) {
            fprintf(stderr, "short read on %s\n", path); exit(1);
        }
        fclose(f); return buf;
    };
    float *active   = read_floats(argv[6], (size_t)n_active * 3);
    float *static_p = read_floats(argv[7], (size_t)n_static * 3);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_dist = make_pso(ctx, "dist_active_static");
        id<MTLComputePipelineState> pso_wp   = make_pso(ctx, "wpoly6_inplace");
        id<MTLComputePipelineState> pso_rs   = make_pso(ctx, "rowsum_density");

        size_t r2_bytes = (size_t)n_r2 * sizeof(float);
        size_t dens_bytes = (size_t)n_active * sizeof(float);

        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:active
            length:(size_t)n_active * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bS = [ctx.device newBufferWithBytes:static_p
            length:(size_t)n_static * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2 = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bW = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];

        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNs = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNr2 = [ctx.device newBufferWithBytes:&n_r2 length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // distance
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_dist];
          [e setBuffer:bA offset:0 atIndex:0]; [e setBuffer:bS offset:0 atIndex:1];
          [e setBuffer:bR2 offset:0 atIndex:2];
          [e setBuffer:bNa offset:0 atIndex:3]; [e setBuffer:bNs offset:0 atIndex:4];
          [e dispatchThreads:MTLSizeMake(n_active, n_static, 1)
              threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
          [e endEncoding]; }
        // GPU-side blit r2 → W (so we can later transform W in-place
        // while keeping r2 saved for backward).
        { id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
          [blit copyFromBuffer:bR2 sourceOffset:0
                      toBuffer:bW destinationOffset:0
                          size:r2_bytes];
          [blit endEncoding]; }
        // wpoly6 in-place on W
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_wp];
          [e setBuffer:bW offset:0 atIndex:0]; [e setBuffer:bH2 offset:0 atIndex:1];
          [e setBuffer:bP6 offset:0 atIndex:2]; [e setBuffer:bNr2 offset:0 atIndex:3];
          [e dispatchThreads:MTLSizeMake(n_r2, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        // rowsum density
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_rs];
          [e setBuffer:bW offset:0 atIndex:0]; [e setBuffer:bD offset:0 atIndex:1];
          [e setBuffer:bM offset:0 atIndex:2];
          [e setBuffer:bNs offset:0 atIndex:3]; [e setBuffer:bNa offset:0 atIndex:4];
          [e dispatchThreads:MTLSizeMake(256, n_active, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fd = fopen(argv[8], "wb");
        fwrite([bD contents], 1, dens_bytes, fd); fclose(fd);
        FILE *fr = fopen(argv[9], "wb");
        fwrite([bR2 contents], 1, r2_bytes, fr); fclose(fr);
    }
    free(active); free(static_p);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.A — density_as_bwd: backward through density_as_fwd.
//
// Inputs: ∂L/∂density [n_active], saved r², positions
// Outputs: ∂L/∂active [n_active × 3]  (static positions are frozen)
//
// Reverse pipeline:
//   ∂L/∂W   = rowsum_density_backward(∂L/∂density)         (broadcast)
//   ∂L/∂r²  = wpoly6_inplace_backward(∂L/∂W; saved r²)     (in-place)
//   ∂L/∂act = dist_active_static_backward(∂L/∂r²)
// ──────────────────────────────────────────────────────────────────────
static int run_density_as_bwd(int argc, char **argv) {
    if (argc != 11) {
        fprintf(stderr,
            "usage: sib_metal density_as_bwd "
            "<n_active> <n_static> <h> <mass> "
            "<active.bin> <static.bin> <r2_saved.bin> "
            "<grad_density.bin> <grad_active_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h    = (float)atof(argv[4]);
    float mass = (float)atof(argv[5]);

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    uint32_t n_r2 = n_active * n_static;

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        fread(buf, sizeof(float), n_floats, f); fclose(f); return buf;
    };
    float *active   = read_floats(argv[6], (size_t)n_active * 3);
    float *static_p = read_floats(argv[7], (size_t)n_static * 3);
    float *r2_saved = read_floats(argv[8], (size_t)n_r2);
    float *grad_d   = read_floats(argv[9], (size_t)n_active);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_rs_bw = make_pso(ctx, "rowsum_density_backward");
        id<MTLComputePipelineState> pso_wp_bw = make_pso(ctx, "wpoly6_inplace_backward");
        id<MTLComputePipelineState> pso_d_bw  = make_pso(ctx, "dist_active_static_backward");

        size_t r2_bytes = (size_t)n_r2 * sizeof(float);
        size_t a_bytes  = (size_t)n_active * 3 * sizeof(float);

        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:active length:a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bS = [ctx.device newBufferWithBytes:static_p
            length:(size_t)n_static * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2 = [ctx.device newBufferWithBytes:r2_saved length:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGd = [ctx.device newBufferWithBytes:grad_d
            length:(size_t)n_active * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGW_or_Gr2 = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGa = [ctx.device newBufferWithLength:a_bytes
            options:MTLResourceStorageModeShared];
        memset([bGa contents], 0, a_bytes);

        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNs = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNr2 = [ctx.device newBufferWithBytes:&n_r2 length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // (1) rowsum_bwd: broadcast grad_density → grad_W (write)
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_rs_bw];
          [e setBuffer:bGd offset:0 atIndex:0];
          [e setBuffer:bGW_or_Gr2 offset:0 atIndex:1];
          [e setBuffer:bM offset:0 atIndex:2];
          [e setBuffer:bNs offset:0 atIndex:3];
          [e dispatchThreads:MTLSizeMake(n_static, n_active, 1)
              threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
          [e endEncoding]; }
        // (2) wpoly6_bwd in-place: grad_W → grad_r²
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_wp_bw];
          [e setBuffer:bR2 offset:0 atIndex:0];
          [e setBuffer:bGW_or_Gr2 offset:0 atIndex:1];
          [e setBuffer:bH2 offset:0 atIndex:2];
          [e setBuffer:bP6 offset:0 atIndex:3];
          [e setBuffer:bNr2 offset:0 atIndex:4];
          [e dispatchThreads:MTLSizeMake(n_r2, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        // (3) dist_bwd: grad_r² → grad_active
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_d_bw];
          [e setBuffer:bA offset:0 atIndex:0];
          [e setBuffer:bS offset:0 atIndex:1];
          [e setBuffer:bGW_or_Gr2 offset:0 atIndex:2];
          [e setBuffer:bGa offset:0 atIndex:3];
          [e setBuffer:bNa offset:0 atIndex:4];
          [e setBuffer:bNs offset:0 atIndex:5];
          [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
          [e endEncoding]; }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fo = fopen("/tmp/density_as_grad_active.bin", "wb");
        fwrite([bGa contents], 1, a_bytes, fo); fclose(fo);
        // Allow caller to override output path via argv[9] if it's
        // different — for simplicity fixed path.
    }
    free(active); free(static_p); free(r2_saved); free(grad_d);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.B — density_aa_fwd: density from active×active SPH neighbors only.
// (No static/boundary contribution; for testing dist_active_active_bwd
// in isolation. Self-contribution at r=0 is included.)
// ──────────────────────────────────────────────────────────────────────
static int run_density_aa_fwd(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr,
            "usage: sib_metal density_aa_fwd "
            "<n_active> <h> <mass> "
            "<active.bin> <density_out.bin> <r2_saved_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    float h    = (float)atof(argv[3]);
    float mass = (float)atof(argv[4]);
    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    uint32_t n_r2 = n_active * n_active;

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        fread(buf, sizeof(float), n_floats, f); fclose(f); return buf;
    };
    float *active = read_floats(argv[5], (size_t)n_active * 3);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_dist = make_pso(ctx, "dist_active_active");
        id<MTLComputePipelineState> pso_wp   = make_pso(ctx, "wpoly6_inplace");
        id<MTLComputePipelineState> pso_rs   = make_pso(ctx, "rowsum_density");

        size_t r2_bytes = (size_t)n_r2 * sizeof(float);
        size_t dens_bytes = (size_t)n_active * sizeof(float);
        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:active
            length:(size_t)n_active * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2 = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bW = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNr2 = [ctx.device newBufferWithBytes:&n_r2 length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // dist_active_active
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_dist];
          [e setBuffer:bA offset:0 atIndex:0]; [e setBuffer:bR2 offset:0 atIndex:1];
          [e setBuffer:bNa offset:0 atIndex:2];
          [e dispatchThreads:MTLSizeMake(n_active, n_active, 1)
              threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
          [e endEncoding]; }
        // blit r2 → W (preserve r2 for backward)
        { id<MTLBlitCommandEncoder> b = [cmd blitCommandEncoder];
          [b copyFromBuffer:bR2 sourceOffset:0 toBuffer:bW destinationOffset:0
                       size:r2_bytes];
          [b endEncoding]; }
        // wpoly6 in-place on W
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_wp];
          [e setBuffer:bW offset:0 atIndex:0]; [e setBuffer:bH2 offset:0 atIndex:1];
          [e setBuffer:bP6 offset:0 atIndex:2]; [e setBuffer:bNr2 offset:0 atIndex:3];
          [e dispatchThreads:MTLSizeMake(n_r2, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        // rowsum density (n_cols = n_rows = n_active)
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_rs];
          [e setBuffer:bW offset:0 atIndex:0]; [e setBuffer:bD offset:0 atIndex:1];
          [e setBuffer:bM offset:0 atIndex:2];
          [e setBuffer:bNa offset:0 atIndex:3]; [e setBuffer:bNa offset:0 atIndex:4];
          [e dispatchThreads:MTLSizeMake(256, n_active, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fd = fopen(argv[6], "wb"); fwrite([bD contents], 1, dens_bytes, fd); fclose(fd);
        FILE *fr = fopen(argv[7], "wb"); fwrite([bR2 contents], 1, r2_bytes, fr); fclose(fr);
    }
    free(active);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.B — density_aa_bwd: backward through density_aa_fwd.
// Reverse: rowsum_bwd → wpoly6_bwd → dist_active_active_bwd.
// ──────────────────────────────────────────────────────────────────────
static int run_density_aa_bwd(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr,
            "usage: sib_metal density_aa_bwd "
            "<n_active> <h> <mass> "
            "<active.bin> <r2_saved.bin> <grad_density.bin> <grad_active_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    float h    = (float)atof(argv[3]);
    float mass = (float)atof(argv[4]);
    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    uint32_t n_r2 = n_active * n_active;

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        fread(buf, sizeof(float), n_floats, f); fclose(f); return buf;
    };
    float *active   = read_floats(argv[5], (size_t)n_active * 3);
    float *r2_saved = read_floats(argv[6], (size_t)n_r2);
    float *grad_d   = read_floats(argv[7], (size_t)n_active);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_rs_bw = make_pso(ctx, "rowsum_density_backward");
        id<MTLComputePipelineState> pso_wp_bw = make_pso(ctx, "wpoly6_inplace_backward");
        id<MTLComputePipelineState> pso_d_bw  = make_pso(ctx, "dist_active_active_backward");

        size_t r2_bytes = (size_t)n_r2 * sizeof(float);
        size_t a_bytes  = (size_t)n_active * 3 * sizeof(float);
        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:active length:a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2 = [ctx.device newBufferWithBytes:r2_saved length:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGd = [ctx.device newBufferWithBytes:grad_d
            length:(size_t)n_active * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGW_or_Gr2 = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGa = [ctx.device newBufferWithLength:a_bytes
            options:MTLResourceStorageModeShared];
        memset([bGa contents], 0, a_bytes);
        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNr2 = [ctx.device newBufferWithBytes:&n_r2 length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // (1) rowsum_bwd: broadcast grad_density → grad_W
        // n_cols=n_active here.
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_rs_bw];
          [e setBuffer:bGd offset:0 atIndex:0];
          [e setBuffer:bGW_or_Gr2 offset:0 atIndex:1];
          [e setBuffer:bM offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
          [e dispatchThreads:MTLSizeMake(n_active, n_active, 1)
              threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
          [e endEncoding]; }
        // (2) wpoly6_bwd in-place: grad_W → grad_r²
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_wp_bw];
          [e setBuffer:bR2 offset:0 atIndex:0]; [e setBuffer:bGW_or_Gr2 offset:0 atIndex:1];
          [e setBuffer:bH2 offset:0 atIndex:2]; [e setBuffer:bP6 offset:0 atIndex:3];
          [e setBuffer:bNr2 offset:0 atIndex:4];
          [e dispatchThreads:MTLSizeMake(n_r2, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        // (3) dist_active_active_bwd: grad_r² → grad_active
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_d_bw];
          [e setBuffer:bA offset:0 atIndex:0]; [e setBuffer:bGW_or_Gr2 offset:0 atIndex:1];
          [e setBuffer:bGa offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
          [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
          [e endEncoding]; }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fo = fopen(argv[8], "wb");
        fwrite([bGa contents], 1, a_bytes, fo); fclose(fo);
    }
    free(active); free(r2_saved); free(grad_d);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.C — density_constraint_grad_backward host op.
// CLI:
//   density_constraint_grad_bwd <n_active> <n_static> <h> <mass> <rho_rest>
//     <active.bin> <static.bin> <r2_aa.bin> <r2_as.bin>
//     <grad_grad_C.bin> <grad_denom_helper.bin> <grad_active_out.bin>
// ──────────────────────────────────────────────────────────────────────
static int run_density_constraint_grad_bwd(int argc, char **argv) {
    if (argc != 14) {
        fprintf(stderr,
            "usage: sib_metal density_constraint_grad_bwd "
            "<n_active> <n_static> <h> <mass> <rho_rest> "
            "<active.bin> <static.bin> <r2_aa.bin> <r2_as.bin> "
            "<grad_grad_C.bin> <grad_denom_helper.bin> <grad_active_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h        = (float)atof(argv[4]);
    float mass     = (float)atof(argv[5]);
    float rho_rest = (float)atof(argv[6]);
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        fread(buf, sizeof(float), n_floats, f); fclose(f); return buf;
    };
    float *active   = read_floats(argv[7], (size_t)n_active * 3);
    float *static_p = read_floats(argv[8], (size_t)n_static * 3);
    float *r2_aa    = read_floats(argv[9], (size_t)n_active * n_active);
    float *r2_as    = read_floats(argv[10], (size_t)n_active * n_static);
    float *g_gC     = read_floats(argv[11], (size_t)n_active * 3);
    float *g_dh     = read_floats(argv[12], (size_t)n_active);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso = make_pso(ctx, "density_constraint_grad_backward");

        size_t a_bytes = (size_t)n_active * 3 * sizeof(float);
        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:active length:a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bS = [ctx.device newBufferWithBytes:static_p
            length:(size_t)n_static * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2aa = [ctx.device newBufferWithBytes:r2_aa
            length:(size_t)n_active * n_active * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2as = [ctx.device newBufferWithBytes:r2_as
            length:(size_t)n_active * n_static * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGgC = [ctx.device newBufferWithBytes:g_gC length:a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGdh = [ctx.device newBufferWithBytes:g_dh
            length:(size_t)n_active * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGa = [ctx.device newBufferWithLength:a_bytes
            options:MTLResourceStorageModeShared];
        memset([bGa contents], 0, a_bytes);

        id<MTLBuffer> bH = [ctx.device newBufferWithBytes:&h length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSp = [ctx.device newBufferWithBytes:&spiky_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNs = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
        [e setComputePipelineState:pso];
        [e setBuffer:bA offset:0 atIndex:0];     [e setBuffer:bS offset:0 atIndex:1];
        [e setBuffer:bR2aa offset:0 atIndex:2];  [e setBuffer:bR2as offset:0 atIndex:3];
        [e setBuffer:bGgC offset:0 atIndex:4];   [e setBuffer:bGdh offset:0 atIndex:5];
        [e setBuffer:bGa offset:0 atIndex:6];
        [e setBuffer:bH offset:0 atIndex:7];     [e setBuffer:bSp offset:0 atIndex:8];
        [e setBuffer:bM offset:0 atIndex:9];     [e setBuffer:bR offset:0 atIndex:10];
        [e setBuffer:bNa offset:0 atIndex:11];   [e setBuffer:bNs offset:0 atIndex:12];
        [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
        [e endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fo = fopen(argv[13], "wb");
        fwrite([bGa contents], 1, a_bytes, fo); fclose(fo);
    }
    free(active); free(static_p); free(r2_aa); free(r2_as);
    free(g_gC); free(g_dh);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.D — solve_density_constraint isolated forward + backward ops.
//
// solve_density_constraint_fwd:
//   inputs: pos_pre, lambda_pre, density, grad_C, denom_helper,
//           constants (rho_rest, mass, alpha_inv_dt2)
//   outputs: pos_post, lambda_post (both modified in place)
//   This is the same kernel xpbd_step uses internally, exposed for
//   isolated testing.
//
// solve_density_constraint_bwd:
//   inputs: density, grad_C, denom_helper, lambda_pre (FROM FORWARD),
//           grad_pos_post, grad_lambda_post (gradient seeds)
//   outputs: grad_pos_pre, grad_lambda_pre, grad_density, grad_grad_C,
//            grad_denom_h, grad_rho_rest (per-particle, host sums)
// ──────────────────────────────────────────────────────────────────────
static int run_solve_density_constraint_fwd(int argc, char **argv) {
    if (argc != 13) {
        fprintf(stderr,
            "usage: sib_metal solve_density_constraint_fwd "
            "<n_active> <rho_rest> <mass> <alpha_inv_dt2> "
            "<pos_pre.bin> <lambda_pre.bin> <density.bin> "
            "<grad_C.bin> <denom_helper.bin> "
            "<pos_post_out.bin> <lambda_post_out.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    float rho_rest = (float)atof(argv[3]);
    float mass     = (float)atof(argv[4]);
    float adt2     = (float)atof(argv[5]);

    auto rd = ^(const char *p, size_t n_floats) {
        FILE *f = fopen(p, "rb");
        if (!f) { fprintf(stderr, "open %s\n", p); exit(1); }
        float *b = (float *)malloc(n_floats * sizeof(float));
        fread(b, sizeof(float), n_floats, f); fclose(f); return b;
    };
    float *pos = rd(argv[6], (size_t)n * 3);
    float *lam = rd(argv[7], (size_t)n);
    float *dens = rd(argv[8], (size_t)n);
    float *gC = rd(argv[9], (size_t)n * 3);
    float *dh = rd(argv[10], (size_t)n);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso = make_pso(ctx, "solve_density_constraint");

        size_t pos_b = (size_t)n * 3 * sizeof(float);
        size_t s_b = (size_t)n * sizeof(float);
        id<MTLBuffer> bP = [ctx.device newBufferWithBytes:pos length:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bL = [ctx.device newBufferWithBytes:lam length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD = [ctx.device newBufferWithBytes:dens length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bG = [ctx.device newBufferWithBytes:gC length:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH = [ctx.device newBufferWithBytes:dh length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:&adt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
        [e setComputePipelineState:pso];
        [e setBuffer:bP offset:0 atIndex:0]; [e setBuffer:bL offset:0 atIndex:1];
        [e setBuffer:bD offset:0 atIndex:2]; [e setBuffer:bG offset:0 atIndex:3];
        [e setBuffer:bH offset:0 atIndex:4]; [e setBuffer:bR offset:0 atIndex:5];
        [e setBuffer:bM offset:0 atIndex:6]; [e setBuffer:bA offset:0 atIndex:7];
        [e setBuffer:bN offset:0 atIndex:8];
        [e dispatchThreads:MTLSizeMake(n,1,1)
            threadsPerThreadgroup:MTLSizeMake(64,1,1)];
        [e endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *o1 = fopen(argv[11], "wb"); fwrite([bP contents], 1, pos_b, o1); fclose(o1);
        FILE *o2 = fopen(argv[12], "wb"); fwrite([bL contents], 1, s_b, o2); fclose(o2);
    }
    free(pos); free(lam); free(dens); free(gC); free(dh);
    return 0;
}

static int run_solve_density_constraint_bwd(int argc, char **argv) {
    if (argc != 18) {
        fprintf(stderr,
            "usage: sib_metal solve_density_constraint_bwd "
            "<n_active> <rho_rest> <mass> <alpha_inv_dt2> "
            "<density.bin> <grad_C.bin> <denom_helper.bin> <lambda_pre.bin> "
            "<grad_pos_post.bin> <grad_lambda_post.bin> "
            "<grad_pos_pre_out.bin> <grad_lambda_pre_out.bin> "
            "<grad_density_out.bin> <grad_grad_C_out.bin> "
            "<grad_denom_helper_out.bin> <grad_rho_rest_per_particle.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    float rho_rest = (float)atof(argv[3]);
    float mass     = (float)atof(argv[4]);
    float adt2     = (float)atof(argv[5]);

    auto rd = ^(const char *p, size_t n_floats) {
        FILE *f = fopen(p, "rb");
        if (!f) { fprintf(stderr, "open %s\n", p); exit(1); }
        float *b = (float *)malloc(n_floats * sizeof(float));
        fread(b, sizeof(float), n_floats, f); fclose(f); return b;
    };
    float *dens = rd(argv[6], (size_t)n);
    float *gC   = rd(argv[7], (size_t)n * 3);
    float *dh   = rd(argv[8], (size_t)n);
    float *lam  = rd(argv[9], (size_t)n);
    float *gP   = rd(argv[10], (size_t)n * 3);
    float *gL   = rd(argv[11], (size_t)n);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso = make_pso(ctx, "solve_density_constraint_backward");

        size_t pos_b = (size_t)n * 3 * sizeof(float);
        size_t s_b = (size_t)n * sizeof(float);
        id<MTLBuffer> bD = [ctx.device newBufferWithBytes:dens length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGc = [ctx.device newBufferWithBytes:gC length:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH = [ctx.device newBufferWithBytes:dh length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bL = [ctx.device newBufferWithBytes:lam length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGp = [ctx.device newBufferWithBytes:gP length:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGl = [ctx.device newBufferWithBytes:gL length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGpout = [ctx.device newBufferWithLength:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGlout = [ctx.device newBufferWithLength:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGdout = [ctx.device newBufferWithLength:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGgcout = [ctx.device newBufferWithLength:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGhout = [ctx.device newBufferWithLength:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGrout = [ctx.device newBufferWithLength:s_b
            options:MTLResourceStorageModeShared];
        memset([bGpout contents], 0, pos_b);
        // Other outputs are written (not accumulated) so no need to zero.

        id<MTLBuffer> bR = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:&adt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
        [e setComputePipelineState:pso];
        [e setBuffer:bD offset:0 atIndex:0];     [e setBuffer:bGc offset:0 atIndex:1];
        [e setBuffer:bH offset:0 atIndex:2];     [e setBuffer:bL offset:0 atIndex:3];
        [e setBuffer:bGp offset:0 atIndex:4];    [e setBuffer:bGl offset:0 atIndex:5];
        [e setBuffer:bGpout offset:0 atIndex:6]; [e setBuffer:bGlout offset:0 atIndex:7];
        [e setBuffer:bGdout offset:0 atIndex:8]; [e setBuffer:bGgcout offset:0 atIndex:9];
        [e setBuffer:bGhout offset:0 atIndex:10]; [e setBuffer:bGrout offset:0 atIndex:11];
        [e setBuffer:bR offset:0 atIndex:12];    [e setBuffer:bM offset:0 atIndex:13];
        [e setBuffer:bA offset:0 atIndex:14];    [e setBuffer:bN offset:0 atIndex:15];
        [e dispatchThreads:MTLSizeMake(n,1,1)
            threadsPerThreadgroup:MTLSizeMake(64,1,1)];
        [e endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *o1 = fopen(argv[12], "wb"); fwrite([bGpout contents], 1, pos_b, o1); fclose(o1);
        FILE *o2 = fopen(argv[13], "wb"); fwrite([bGlout contents], 1, s_b, o2); fclose(o2);
        FILE *o3 = fopen(argv[14], "wb"); fwrite([bGdout contents], 1, s_b, o3); fclose(o3);
        FILE *o4 = fopen(argv[15], "wb"); fwrite([bGgcout contents], 1, pos_b, o4); fclose(o4);
        FILE *o5 = fopen(argv[16], "wb"); fwrite([bGhout contents], 1, s_b, o5); fclose(o5);
        FILE *o6 = fopen(argv[17], "wb"); fwrite([bGrout contents], 1, s_b, o6); fclose(o6);
    }
    free(dens); free(gC); free(dh); free(lam); free(gP); free(gL);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// Multi-step XPBD forward+backward (xpbd_full_fwd / _bwd).
//
// Pipeline per step k (n_iters=1 for now to keep state simple):
//   predict_positions
//   density chain: dist_aa + dist_as + wpoly6×2 + rowsum×2 + add
//   density_constraint_grad → grad_C, denom_helper
//   solve_density_constraint
//   update_velocities
//
// Saved state per step (all on host, per-step indexing):
//   x_old[k], v_old[k]          ← 6 floats per particle per step
//   density[k]                   ← 1 float per particle per step
//   grad_C[k]                    ← 3 floats per particle per step
//   denom_helper[k]              ← 1 float per particle per step
//
// r²_aa, r²_as are RECOMPUTED in backward (cheap distance ops) to
// avoid memory blowup at demo1 scale (r²_as = 24 MB per step).
//
// No bonds, no floor in this version — to be added as separate
// extensions once base multi-step backward validates.
// ──────────────────────────────────────────────────────────────────────
static int run_xpbd_full_fwd(int argc, char **argv) {
    if (argc != 15) {
        fprintf(stderr,
            "usage: sib_metal xpbd_full_fwd "
            "<n_active> <n_static> <K> <h> <mass> <rho_rest> <dt> "
            "<gravity_y> <alpha_density> "
            "<pos0.bin> <vel0.bin> <pos_static.bin> <state_out.bin>\n"
            "  state_out.bin contains the per-step trajectory: \n"
            "    K × [x_old(n*3) + v_old(n*3) + density(n) + grad_C(n*3) + denom_h(n)]\n"
            "    Plus K+1 frames of x_post for the trajectory.\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    uint32_t K        = (uint32_t)atoi(argv[4]);
    float h           = (float)atof(argv[5]);
    float mass        = (float)atof(argv[6]);
    float rho_rest    = (float)atof(argv[7]);
    float dt          = (float)atof(argv[8]);
    float g_y         = (float)atof(argv[9]);
    float alpha_dens  = (float)atof(argv[10]);

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));
    float alpha_inv_dt2 = alpha_dens / (dt * dt);

    auto rd = ^(const char *p, size_t n_floats) {
        FILE *f = fopen(p, "rb");
        if (!f) { fprintf(stderr, "open %s\n", p); exit(1); }
        float *b = (float *)malloc(n_floats * sizeof(float));
        fread(b, sizeof(float), n_floats, f); fclose(f); return b;
    };
    float *pos0 = rd(argv[11], (size_t)n_active * 3);
    float *vel0 = rd(argv[12], (size_t)n_active * 3);
    float *pos_static = rd(argv[13], (size_t)n_static * 3);

    size_t pos_b = (size_t)n_active * 3 * sizeof(float);
    size_t s_b = (size_t)n_active * sizeof(float);

    // Per-step state arrays (host).
    // Layout per step: [x_old(3n) + v_old(3n) + density(n) + grad_C(3n) + denom_h(n)]
    size_t per_step_floats = (size_t)n_active * (3 + 3 + 1 + 3 + 1);
    float *state = (float *)calloc((size_t)K * per_step_floats, sizeof(float));
    // Plus K+1 frames of x_post for the trajectory.
    float *traj = (float *)calloc((size_t)(K + 1) * (size_t)n_active * 3, sizeof(float));
    memcpy(traj, pos0, pos_b);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pred  = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_d_aa  = make_pso(ctx, "dist_active_active");
        id<MTLComputePipelineState> pso_d_as  = make_pso(ctx, "dist_active_static");
        id<MTLComputePipelineState> pso_wp    = make_pso(ctx, "wpoly6_inplace");
        id<MTLComputePipelineState> pso_rs    = make_pso(ctx, "rowsum_density");
        id<MTLComputePipelineState> pso_addin = make_pso(ctx, "add_inplace");
        id<MTLComputePipelineState> pso_dgrad = make_pso(ctx, "density_constraint_grad");
        id<MTLComputePipelineState> pso_solv  = make_pso(ctx, "solve_density_constraint");
        id<MTLComputePipelineState> pso_uv    = make_pso(ctx, "update_velocities");

        size_t r2aa_b = (size_t)n_active * n_active * sizeof(float);
        size_t r2as_b = (size_t)n_active * n_static * sizeof(float);

        id<MTLBuffer> bX  = [ctx.device newBufferWithBytes:pos0 length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bV  = [ctx.device newBufferWithBytes:vel0 length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bXp = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSp = [ctx.device newBufferWithBytes:pos_static
            length:(size_t)n_static * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2aa = [ctx.device newBufferWithLength:r2aa_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2as = [ctx.device newBufferWithLength:r2as_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bWaa  = [ctx.device newBufferWithLength:r2aa_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bWas  = [ctx.device newBufferWithLength:r2as_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD_aa = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD    = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGc   = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDh   = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bLam  = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];

        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGy = [ctx.device newBufferWithBytes:&g_y length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH  = [ctx.device newBufferWithBytes:&h length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSk = [ctx.device newBufferWithBytes:&spiky_const length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM  = [ctx.device newBufferWithBytes:&mass length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR  = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bA  = [ctx.device newBufferWithBytes:&alpha_inv_dt2 length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNs = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        uint32_t n_aa_total = n_active * n_active;
        uint32_t n_as_total = n_active * n_static;
        id<MTLBuffer> bNaaTot = [ctx.device newBufferWithBytes:&n_aa_total length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNasTot = [ctx.device newBufferWithBytes:&n_as_total length:sizeof(uint32_t) options:MTLResourceStorageModeShared];

        for (uint32_t k = 0; k < K; k++) {
            // SAVE x_old, v_old (current state at start of step)
            float *step_state = state + (size_t)k * per_step_floats;
            memcpy(step_state, [bX contents], pos_b);
            memcpy(step_state + 3 * n_active, [bV contents], pos_b);

            memset([bLam contents], 0, s_b);

            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            // (1) predict
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pred];
              [e setBuffer:bX offset:0 atIndex:0]; [e setBuffer:bV offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bGy offset:0 atIndex:4]; [e setBuffer:bNa offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (2) dist_aa
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_d_aa];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bR2aa offset:0 atIndex:1];
              [e setBuffer:bNa offset:0 atIndex:2];
              [e dispatchThreads:MTLSizeMake(n_active, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(16,16,1)];
              [e endEncoding]; }
            // (3) dist_as
            if (n_static > 0) {
              id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_d_as];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bSp offset:0 atIndex:1];
              [e setBuffer:bR2as offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
              [e setBuffer:bNs offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n_active, n_static, 1)
                  threadsPerThreadgroup:MTLSizeMake(16,16,1)];
              [e endEncoding];
            }
            // (4) Wpoly6 on aa
            { id<MTLBlitCommandEncoder> b = [cmd blitCommandEncoder];
              [b copyFromBuffer:bR2aa sourceOffset:0 toBuffer:bWaa destinationOffset:0 size:r2aa_b];
              [b endEncoding]; }
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_wp];
              [e setBuffer:bWaa offset:0 atIndex:0]; [e setBuffer:bH2 offset:0 atIndex:1];
              [e setBuffer:bP6 offset:0 atIndex:2]; [e setBuffer:bNaaTot offset:0 atIndex:3];
              [e dispatchThreads:MTLSizeMake(n_aa_total,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
              [e endEncoding]; }
            // (5) Wpoly6 on as
            if (n_static > 0) {
              { id<MTLBlitCommandEncoder> b = [cmd blitCommandEncoder];
                [b copyFromBuffer:bR2as sourceOffset:0 toBuffer:bWas destinationOffset:0 size:r2as_b];
                [b endEncoding]; }
              { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_wp];
                [e setBuffer:bWas offset:0 atIndex:0]; [e setBuffer:bH2 offset:0 atIndex:1];
                [e setBuffer:bP6 offset:0 atIndex:2]; [e setBuffer:bNasTot offset:0 atIndex:3];
                [e dispatchThreads:MTLSizeMake(n_as_total,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
                [e endEncoding]; }
            }
            // (6) rowsum aa → density_aa
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_rs];
              [e setBuffer:bWaa offset:0 atIndex:0]; [e setBuffer:bD_aa offset:0 atIndex:1];
              [e setBuffer:bM offset:0 atIndex:2];
              [e setBuffer:bNa offset:0 atIndex:3]; [e setBuffer:bNa offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(256, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(256,1,1)];
              [e endEncoding]; }
            // (7) rowsum as → density (will accumulate via add)
            if (n_static > 0) {
              id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_rs];
              [e setBuffer:bWas offset:0 atIndex:0]; [e setBuffer:bD offset:0 atIndex:1];
              [e setBuffer:bM offset:0 atIndex:2];
              [e setBuffer:bNs offset:0 atIndex:3]; [e setBuffer:bNa offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(256, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(256,1,1)];
              [e endEncoding];
              // density += density_aa
              { id<MTLComputeCommandEncoder> e2 = [cmd computeCommandEncoder];
                [e2 setComputePipelineState:pso_addin];
                [e2 setBuffer:bD offset:0 atIndex:0]; [e2 setBuffer:bD_aa offset:0 atIndex:1];
                [e2 setBuffer:bNa offset:0 atIndex:2];
                [e2 dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
                [e2 endEncoding]; }
            } else {
              // density = density_aa
              id<MTLBlitCommandEncoder> b = [cmd blitCommandEncoder];
              [b copyFromBuffer:bD_aa sourceOffset:0 toBuffer:bD destinationOffset:0 size:s_b];
              [b endEncoding];
            }
            // (8) density_constraint_grad
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_dgrad];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bSp offset:0 atIndex:1];
              [e setBuffer:bR2aa offset:0 atIndex:2]; [e setBuffer:bR2as offset:0 atIndex:3];
              [e setBuffer:bGc offset:0 atIndex:4]; [e setBuffer:bDh offset:0 atIndex:5];
              [e setBuffer:bH offset:0 atIndex:6]; [e setBuffer:bSk offset:0 atIndex:7];
              [e setBuffer:bM offset:0 atIndex:8]; [e setBuffer:bR offset:0 atIndex:9];
              [e setBuffer:bNa offset:0 atIndex:10]; [e setBuffer:bNs offset:0 atIndex:11];
              [e dispatchThreads:MTLSizeMake(256, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(256,1,1)];
              [e endEncoding]; }
            // (9) solve_density_constraint
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_solv];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bLam offset:0 atIndex:1];
              [e setBuffer:bD offset:0 atIndex:2]; [e setBuffer:bGc offset:0 atIndex:3];
              [e setBuffer:bDh offset:0 atIndex:4]; [e setBuffer:bR offset:0 atIndex:5];
              [e setBuffer:bM offset:0 atIndex:6]; [e setBuffer:bA offset:0 atIndex:7];
              [e setBuffer:bNa offset:0 atIndex:8];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (10) update_vel: v_new = (xp - x_old)/dt
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uv];
              [e setBuffer:bV offset:0 atIndex:0]; [e setBuffer:bX offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bNa offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // SAVE density, grad_C, denom_helper for backward
            float *p = step_state + 6 * n_active;
            memcpy(p, [bD contents], s_b);                            // density
            memcpy(p + n_active, [bGc contents], pos_b);              // grad_C
            memcpy(p + n_active + 3 * n_active, [bDh contents], s_b); // denom_h
            // Save x_post in trajectory
            memcpy(traj + (size_t)(k + 1) * n_active * 3,
                   [bXp contents], pos_b);
            // Advance: x_old ← x_post
            memcpy([bX contents], [bXp contents], pos_b);
        }

        // Write outputs: state buffer, trajectory, final velocity
        FILE *fs = fopen(argv[14], "wb");
        fwrite(state, sizeof(float), (size_t)K * per_step_floats, fs);
        fwrite(traj, sizeof(float), (size_t)(K + 1) * n_active * 3, fs);
        fwrite([bV contents], 1, pos_b, fs);
        fclose(fs);
    }
    free(pos0); free(vel0); free(pos_static); free(state); free(traj);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// Multi-step backward driver. Walks K steps in reverse, calling all
// the M9 backward kernels per step, accumulating ∂L/∂(rho_rest, x_init,
// v_init).
//
// Inputs:
//   state file (from xpbd_full_fwd)
//   ∂L/∂x_final (input gradient seed)
// Outputs:
//   ∂L/∂x_init, ∂L/∂v_init, ∂L/∂rho_rest (scalar)
// ──────────────────────────────────────────────────────────────────────
static int run_xpbd_full_bwd(int argc, char **argv) {
    if (argc != 17) {
        fprintf(stderr,
            "usage: sib_metal xpbd_full_bwd "
            "<n_active> <n_static> <K> <h> <mass> <rho_rest> <dt> "
            "<gravity_y> <alpha_density> "
            "<state_in.bin> <pos_static.bin> <grad_x_final.bin> "
            "<grad_x_init_out.bin> <grad_v_init_out.bin> <grad_rho_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    uint32_t K        = (uint32_t)atoi(argv[4]);
    float h           = (float)atof(argv[5]);
    float mass        = (float)atof(argv[6]);
    float rho_rest    = (float)atof(argv[7]);
    float dt          = (float)atof(argv[8]);
    float g_y         = (float)atof(argv[9]);
    float alpha_dens  = (float)atof(argv[10]);

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));
    float alpha_inv_dt2 = alpha_dens / (dt * dt);

    auto rd = ^(const char *p, size_t n_floats) {
        FILE *f = fopen(p, "rb");
        if (!f) { fprintf(stderr, "open %s\n", p); exit(1); }
        float *b = (float *)malloc(n_floats * sizeof(float));
        fread(b, sizeof(float), n_floats, f); fclose(f); return b;
    };
    // State file layout (from xpbd_full_fwd):
    //   [K × per_step_floats] state
    //   [(K+1) × n*3] traj
    //   [n*3] vel_final
    size_t per_step_floats = (size_t)n_active * (3 + 3 + 1 + 3 + 1);
    size_t state_size = (size_t)K * per_step_floats
                      + (size_t)(K + 1) * n_active * 3
                      + (size_t)n_active * 3;
    float *all = rd(argv[11], state_size);
    float *state = all;
    float *traj = all + (size_t)K * per_step_floats;
    // (vel_final at all + ... + (K+1)*n*3 — not used here)

    float *pos_static = rd(argv[12], (size_t)n_static * 3);
    float *grad_x_fin = rd(argv[13], (size_t)n_active * 3);

    size_t pos_b = (size_t)n_active * 3 * sizeof(float);
    size_t s_b = (size_t)n_active * sizeof(float);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        // Forward kernels (recompute distances and predict x_pre):
        id<MTLComputePipelineState> pso_pred  = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_d_aa  = make_pso(ctx, "dist_active_active");
        id<MTLComputePipelineState> pso_d_as  = make_pso(ctx, "dist_active_static");
        // Backward kernels:
        id<MTLComputePipelineState> pso_uvbw  = make_pso(ctx, "update_velocities_backward");
        id<MTLComputePipelineState> pso_solv_bw = make_pso(ctx, "solve_density_constraint_backward");
        id<MTLComputePipelineState> pso_dgrad_bw = make_pso(ctx, "density_constraint_grad_backward");
        id<MTLComputePipelineState> pso_rs_bw = make_pso(ctx, "rowsum_density_backward");
        id<MTLComputePipelineState> pso_wp_bw = make_pso(ctx, "wpoly6_inplace_backward");
        id<MTLComputePipelineState> pso_d_aa_bw = make_pso(ctx, "dist_active_active_backward");
        id<MTLComputePipelineState> pso_d_as_bw = make_pso(ctx, "dist_active_static_backward");
        id<MTLComputePipelineState> pso_pred_bw = make_pso(ctx, "predict_positions_backward");

        // Persistent buffers.
        id<MTLBuffer> bX  = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bV  = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bXp = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSp = [ctx.device newBufferWithBytes:pos_static
            length:(size_t)n_static * 3 * sizeof(float) options:MTLResourceStorageModeShared];
        size_t r2aa_b = (size_t)n_active * n_active * sizeof(float);
        size_t r2as_b = (size_t)n_active * n_static * sizeof(float);
        id<MTLBuffer> bR2aa = [ctx.device newBufferWithLength:r2aa_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2as = [ctx.device newBufferWithLength:r2as_b options:MTLResourceStorageModeShared];

        // Per-step inputs from state:
        id<MTLBuffer> bD  = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGc = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDh = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bLam = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        memset([bLam contents], 0, s_b);  // λ_pre always 0

        // Running gradients:
        id<MTLBuffer> bGx_running = [ctx.device newBufferWithBytes:grad_x_fin length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGv_running = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        memset([bGv_running contents], 0, pos_b);

        // Per-step working gradients:
        id<MTLBuffer> bGv_in = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGx_pred = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGx_old_new = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGv_old_new = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        // Outputs of solve_density_constraint_backward:
        id<MTLBuffer> bGlam_pre = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGdens = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGgC = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGdh = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGrho_per = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        // Density chain backward intermediates:
        id<MTLBuffer> bGW_or_Gr2_aa = [ctx.device newBufferWithLength:r2aa_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGW_or_Gr2_as = [ctx.device newBufferWithLength:r2as_b options:MTLResourceStorageModeShared];
        // For predict_backward:
        id<MTLBuffer> bGgy_per = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];

        // Constants.
        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGy = [ctx.device newBufferWithBytes:&g_y length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH  = [ctx.device newBufferWithBytes:&h length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSk = [ctx.device newBufferWithBytes:&spiky_const length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM  = [ctx.device newBufferWithBytes:&mass length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR  = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bA  = [ctx.device newBufferWithBytes:&alpha_inv_dt2 length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNs = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        uint32_t n_aa_total = n_active * n_active;
        uint32_t n_as_total = n_active * n_static;
        id<MTLBuffer> bNaaTot = [ctx.device newBufferWithBytes:&n_aa_total length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNasTot = [ctx.device newBufferWithBytes:&n_as_total length:sizeof(uint32_t) options:MTLResourceStorageModeShared];

        // Total ρ gradient accumulator (host-side scalar).
        float total_grad_rho = 0.0f;

        // Walk K steps backward.
        for (int32_t k = (int32_t)K - 1; k >= 0; k--) {
            float *step_state = state + (size_t)k * per_step_floats;
            // Load saved state.
            memcpy([bX contents], step_state, pos_b);                  // x_old
            memcpy([bV contents], step_state + 3 * n_active, pos_b);   // v_old
            memcpy([bD contents], step_state + 6 * n_active, s_b);     // density
            memcpy([bGc contents], step_state + 7 * n_active, pos_b);  // grad_C
            memcpy([bDh contents], step_state + 10 * n_active, s_b);   // denom_h

            // Recompute x_pre = x_old + dt²·g + dt·v_old (same as predict).
            id<MTLCommandBuffer> cmd1 = [ctx.queue commandBuffer];
            { id<MTLComputeCommandEncoder> e = [cmd1 computeCommandEncoder];
              [e setComputePipelineState:pso_pred];
              [e setBuffer:bX offset:0 atIndex:0]; [e setBuffer:bV offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bGy offset:0 atIndex:4]; [e setBuffer:bNa offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // Recompute r2_aa, r2_as for the density chain backward.
            { id<MTLComputeCommandEncoder> e = [cmd1 computeCommandEncoder];
              [e setComputePipelineState:pso_d_aa];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bR2aa offset:0 atIndex:1];
              [e setBuffer:bNa offset:0 atIndex:2];
              [e dispatchThreads:MTLSizeMake(n_active, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(16,16,1)];
              [e endEncoding]; }
            if (n_static > 0) {
              id<MTLComputeCommandEncoder> e = [cmd1 computeCommandEncoder];
              [e setComputePipelineState:pso_d_as];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bSp offset:0 atIndex:1];
              [e setBuffer:bR2as offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
              [e setBuffer:bNs offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n_active, n_static, 1)
                  threadsPerThreadgroup:MTLSizeMake(16,16,1)];
              [e endEncoding];
            }
            [cmd1 commit]; [cmd1 waitUntilCompleted];

            memset([bGx_old_new contents], 0, pos_b);
            memset([bGv_old_new contents], 0, pos_b);
            memset([bGx_pred contents], 0, pos_b);
            // Move running v gradient into bGv_in for update_vel_bwd.
            memcpy([bGv_in contents], [bGv_running contents], pos_b);

            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            // (a) update_vel backward: ∂L/∂v_new (=bGv_in) → +bGx_running (∂L/∂x_post), +bGx_old_new (∂L/∂x_old)
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uvbw];
              [e setBuffer:bGv_in offset:0 atIndex:0];
              [e setBuffer:bGx_running offset:0 atIndex:1];
              [e setBuffer:bGx_old_new offset:0 atIndex:2];
              [e setBuffer:bDt offset:0 atIndex:3]; [e setBuffer:bNa offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (b) solve_density_constraint backward
            // ∂L/∂λ_post = 0 (single-iter, no propagation)
            // We use bGdens as a SCRATCH for the zeros.
            // Actually solve bw expects buf at index 5 to be ∂L/∂λ_post — pass bDh re-purposed? No, need a zero buffer.
            // Re-use bGdh for this purpose (we'll overwrite it with output ∂L/∂dh anyway).
            // But we need a SEPARATE zero buffer. Use a small temp.
            // Quick fix: just zero bGdh first and pass it. But bGdh gets overwritten by the kernel.
            // Cleaner: allocate a separate zero buffer for ∂L/∂λ_post.
            // Hack: bLam itself is zeros (we set above) — reuse it.
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_solv_bw];
              [e setBuffer:bD offset:0 atIndex:0];     [e setBuffer:bGc offset:0 atIndex:1];
              [e setBuffer:bDh offset:0 atIndex:2];    [e setBuffer:bLam offset:0 atIndex:3];
              [e setBuffer:bGx_running offset:0 atIndex:4];  // ∂L/∂x_post (input)
              [e setBuffer:bLam offset:0 atIndex:5];   // ∂L/∂λ_post = 0 (reuse bLam as zeros)
              [e setBuffer:bGx_pred offset:0 atIndex:6];  // ∂L/∂x_pre (accum)
              [e setBuffer:bGlam_pre offset:0 atIndex:7];
              [e setBuffer:bGdens offset:0 atIndex:8]; [e setBuffer:bGgC offset:0 atIndex:9];
              [e setBuffer:bGdh offset:0 atIndex:10];  [e setBuffer:bGrho_per offset:0 atIndex:11];
              [e setBuffer:bR offset:0 atIndex:12];    [e setBuffer:bM offset:0 atIndex:13];
              [e setBuffer:bA offset:0 atIndex:14];    [e setBuffer:bNa offset:0 atIndex:15];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (c) density_constraint_grad backward: contributes ∂L/∂x via grad_C and denom_helper
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_dgrad_bw];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bSp offset:0 atIndex:1];
              [e setBuffer:bR2aa offset:0 atIndex:2]; [e setBuffer:bR2as offset:0 atIndex:3];
              [e setBuffer:bGgC offset:0 atIndex:4]; [e setBuffer:bGdh offset:0 atIndex:5];
              [e setBuffer:bGx_pred offset:0 atIndex:6];   // accumulate
              [e setBuffer:bH offset:0 atIndex:7]; [e setBuffer:bSk offset:0 atIndex:8];
              [e setBuffer:bM offset:0 atIndex:9]; [e setBuffer:bR offset:0 atIndex:10];
              [e setBuffer:bNa offset:0 atIndex:11]; [e setBuffer:bNs offset:0 atIndex:12];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (d) density chain backward: ∂L/∂density → ∂L/∂x_pre (additional contribution)
            //     Pipeline: rowsum_bwd → wpoly6_bwd → dist_aa_bwd
            //               rowsum_bwd → wpoly6_bwd → dist_as_bwd
            //     Both contribute to bGx_pred (same buffer, multiple paths add up via dist_*_backward kernels).
            // density_aa branch:
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_rs_bw];
              [e setBuffer:bGdens offset:0 atIndex:0];
              [e setBuffer:bGW_or_Gr2_aa offset:0 atIndex:1];
              [e setBuffer:bM offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
              [e dispatchThreads:MTLSizeMake(n_active, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(16,16,1)];
              [e endEncoding]; }
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_wp_bw];
              [e setBuffer:bR2aa offset:0 atIndex:0]; [e setBuffer:bGW_or_Gr2_aa offset:0 atIndex:1];
              [e setBuffer:bH2 offset:0 atIndex:2]; [e setBuffer:bP6 offset:0 atIndex:3];
              [e setBuffer:bNaaTot offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n_aa_total,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
              [e endEncoding]; }
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_d_aa_bw];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bGW_or_Gr2_aa offset:0 atIndex:1];
              [e setBuffer:bGx_pred offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // density_as branch:
            if (n_static > 0) {
              { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_rs_bw];
                [e setBuffer:bGdens offset:0 atIndex:0];
                [e setBuffer:bGW_or_Gr2_as offset:0 atIndex:1];
                [e setBuffer:bM offset:0 atIndex:2]; [e setBuffer:bNs offset:0 atIndex:3];
                [e dispatchThreads:MTLSizeMake(n_static, n_active, 1)
                    threadsPerThreadgroup:MTLSizeMake(16,16,1)];
                [e endEncoding]; }
              { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_wp_bw];
                [e setBuffer:bR2as offset:0 atIndex:0]; [e setBuffer:bGW_or_Gr2_as offset:0 atIndex:1];
                [e setBuffer:bH2 offset:0 atIndex:2]; [e setBuffer:bP6 offset:0 atIndex:3];
                [e setBuffer:bNasTot offset:0 atIndex:4];
                [e dispatchThreads:MTLSizeMake(n_as_total,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
                [e endEncoding]; }
              { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_d_as_bw];
                [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bSp offset:0 atIndex:1];
                [e setBuffer:bGW_or_Gr2_as offset:0 atIndex:2];
                [e setBuffer:bGx_pred offset:0 atIndex:3];
                [e setBuffer:bNa offset:0 atIndex:4]; [e setBuffer:bNs offset:0 atIndex:5];
                [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
                [e endEncoding]; }
            }
            // (e) predict_positions backward: ∂L/∂x_pred → +bGx_old_new, +bGv_old_new, ∂L/∂g_y per particle
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pred_bw];
              [e setBuffer:bGx_pred offset:0 atIndex:0];
              [e setBuffer:bGx_old_new offset:0 atIndex:1];
              [e setBuffer:bGv_old_new offset:0 atIndex:2];
              [e setBuffer:bGgy_per offset:0 atIndex:3];
              [e setBuffer:bDt offset:0 atIndex:4]; [e setBuffer:bNa offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // Sum ρ gradient (kernel partial + implicit via grad_C).
            //   implicit = -(grad_grad_C · grad_C) / ρ_rest
            //   Note: bGgC currently holds ∂L/∂grad_C from solve_dens_bwd; bGc holds grad_C (forward saved).
            float kernel_rho = 0.0f;
            float *gr = (float *)[bGrho_per contents];
            for (uint32_t i = 0; i < n_active; i++) kernel_rho += gr[i];
            float implicit_rho = 0.0f;
            float *ggC = (float *)[bGgC contents];
            float *gC_fwd = (float *)[bGc contents];
            for (uint32_t i = 0; i < n_active; i++) {
                float dot = 0.0f;
                for (int ax = 0; ax < 3; ax++)
                    dot += ggC[i * 3 + ax] * gC_fwd[i * 3 + ax];
                implicit_rho -= dot / rho_rest;
            }
            total_grad_rho += kernel_rho + implicit_rho;

            // Promote per-step grads to running for previous step.
            memcpy([bGx_running contents], [bGx_old_new contents], pos_b);
            memcpy([bGv_running contents], [bGv_old_new contents], pos_b);
        }

        FILE *o1 = fopen(argv[14], "wb"); fwrite([bGx_running contents], 1, pos_b, o1); fclose(o1);
        FILE *o2 = fopen(argv[15], "wb"); fwrite([bGv_running contents], 1, pos_b, o2); fclose(o2);
        FILE *o3 = fopen(argv[16], "wb"); fwrite(&total_grad_rho, 1, sizeof(float), o3); fclose(o3);
    }
    free(all); free(pos_static); free(grad_x_fin);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: sib_metal <op> [args...]\n"
                        "ops: dist_active_static, dist_active_active, "
                        "wpoly6_inplace, rowsum_density, "
                        "density_constraint_grad, density_constraint_grad_bwd, "
                        "solve_density_constraint_fwd, solve_density_constraint_bwd, "
                        "xpbd_step, xpbd_full_fwd, xpbd_full_bwd, "
                        "step_simple_fwd, step_simple_bwd, "
                        "step_floor_fwd, step_floor_bwd, "
                        "step_bond_fwd, step_bond_bwd, "
                        "density_as_fwd, density_as_bwd, "
                        "density_aa_fwd, density_aa_bwd\n");
        return 1;
    }
    if (strcmp(argv[1], "dist_active_static") == 0)
        return run_dist_active_static(argc, argv);
    if (strcmp(argv[1], "dist_active_active") == 0)
        return run_dist_active_active(argc, argv);
    if (strcmp(argv[1], "wpoly6_inplace") == 0)
        return run_wpoly6_inplace(argc, argv);
    if (strcmp(argv[1], "rowsum_density") == 0)
        return run_rowsum_density(argc, argv);
    if (strcmp(argv[1], "density_constraint_grad") == 0)
        return run_density_constraint_grad(argc, argv);
    if (strcmp(argv[1], "xpbd_step") == 0)
        return run_xpbd_step(argc, argv);
    if (strcmp(argv[1], "step_simple_fwd") == 0)
        return run_step_simple_fwd(argc, argv);
    if (strcmp(argv[1], "step_simple_bwd") == 0)
        return run_step_simple_bwd(argc, argv);
    if (strcmp(argv[1], "step_floor_fwd") == 0)
        return run_step_floor_fwd(argc, argv);
    if (strcmp(argv[1], "step_floor_bwd") == 0)
        return run_step_floor_bwd(argc, argv);
    if (strcmp(argv[1], "step_bond_fwd") == 0)
        return run_step_bond_fwd(argc, argv);
    if (strcmp(argv[1], "step_bond_bwd") == 0)
        return run_step_bond_bwd(argc, argv);
    if (strcmp(argv[1], "density_as_fwd") == 0)
        return run_density_as_fwd(argc, argv);
    if (strcmp(argv[1], "density_as_bwd") == 0)
        return run_density_as_bwd(argc, argv);
    if (strcmp(argv[1], "density_aa_fwd") == 0)
        return run_density_aa_fwd(argc, argv);
    if (strcmp(argv[1], "density_aa_bwd") == 0)
        return run_density_aa_bwd(argc, argv);
    if (strcmp(argv[1], "density_constraint_grad_bwd") == 0)
        return run_density_constraint_grad_bwd(argc, argv);
    if (strcmp(argv[1], "solve_density_constraint_fwd") == 0)
        return run_solve_density_constraint_fwd(argc, argv);
    if (strcmp(argv[1], "solve_density_constraint_bwd") == 0)
        return run_solve_density_constraint_bwd(argc, argv);
    if (strcmp(argv[1], "xpbd_full_fwd") == 0)
        return run_xpbd_full_fwd(argc, argv);
    if (strcmp(argv[1], "xpbd_full_bwd") == 0)
        return run_xpbd_full_bwd(argc, argv);
    fprintf(stderr, "unknown op: %s\n", argv[1]);
    return 1;
}
