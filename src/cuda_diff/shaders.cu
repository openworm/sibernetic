// AUTO-GENERATED from sib_cuda.cu via _split_sib_cuda.py
// All __global__ kernel definitions. Drivers in ops_*.cu launch these
// via the prototypes in shaders.h.

#include "cuda_common.h"
#include "shaders.h"


// ──────────────────────────────────────────────────────────────────────
// M6.1 — wpoly6_inplace
//
// Apply Müller 2003 Wpoly6 SPH smoothing kernel elementwise on a squared
// distance buffer, in-place: input is r², output is W(r²) overwriting
// the same buffer.
//
//     W_poly6(r²) = poly6_const · (h² - r²)³   if r < h else 0
//     poly6_const = 315 / (64π h⁹)             (computed host-side)
//
// One thread per matrix element. Trivially parallel, no atomics. Direct
// translation of src/metal_diff/shaders.metal::wpoly6_inplace.
// ──────────────────────────────────────────────────────────────────────
__global__ void wpoly6_inplace(float *r2_or_W,
                               float h2,
                               float poly6_const,
                               unsigned int n_total)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n_total) return;
    float r2 = r2_or_W[gid];
    if (r2 >= h2) {
        r2_or_W[gid] = 0.0f;
        return;
    }
    float diff = h2 - r2;
    r2_or_W[gid] = poly6_const * diff * diff * diff;
}

// ──────────────────────────────────────────────────────────────────────
// M6.0 — dist_active_static
//   r²[i, j] = ||active[i] - static_p[j]||²
// One thread per (active, static) pair. 2-D grid.
// Direct port of src/metal_diff/shaders.metal::dist_active_static.
// ──────────────────────────────────────────────────────────────────────
__global__ void dist_active_static(const float3 *active,
                                   const float3 *static_p,
                                   float *dist,
                                   unsigned int n_active,
                                   unsigned int n_static)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= n_active || j >= n_static) return;
    float3 a = active[i];
    float3 s = static_p[j];
    float dx = a.x - s.x, dy = a.y - s.y, dz = a.z - s.z;
    dist[i * n_static + j] = dx*dx + dy*dy + dz*dz;
}

// ──────────────────────────────────────────────────────────────────────
// M6.3 — dist_active_active
//   r²[i, j] = ||active[i] - active[j]||²    (diagonal == 0)
// Same as M6.0 but single position buffer used twice. Downstream consumers
// (density_constraint_grad) skip i==j explicitly.
// ──────────────────────────────────────────────────────────────────────
__global__ void dist_active_active(const float3 *active,
                                   float *dist,
                                   unsigned int n_active)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= n_active || j >= n_active) return;
    float3 a = active[i];
    float3 b = active[j];
    float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    dist[i * n_active + j] = dx*dx + dy*dy + dz*dz;
}

// ──────────────────────────────────────────────────────────────────────
// M6.2 — rowsum_density
//   density[i] = mass · Σ_j W[i, j]
// One block per row; 256 threads cooperate via shared-mem tree reduction.
// ──────────────────────────────────────────────────────────────────────
__global__ void rowsum_density(const float *W,
                               float *density,
                               float mass,
                               unsigned int n_cols,
                               unsigned int n_rows)
{
    unsigned int i = blockIdx.x;
    if (i >= n_rows) return;
    __shared__ float partials[256];
    unsigned int t = threadIdx.x;
    unsigned int T = blockDim.x;

    float sum = 0.0f;
    for (unsigned int j = t; j < n_cols; j += T) {
        sum += W[i * n_cols + j];
    }
    partials[t] = sum;
    __syncthreads();

    for (unsigned int stride = T / 2; stride > 0; stride >>= 1) {
        if (t < stride) partials[t] += partials[t + stride];
        __syncthreads();
    }
    if (t == 0) density[i] = mass * partials[0];
}

// ──────────────────────────────────────────────────────────────────────
// M6.4 — density_constraint_grad
//   grad_C[i]      = (mass / ρ_rest) · Σ_j ∇W_spiky(p_i - p_j)
//   denom_helper[i] = Σ_j |∇W_spiky(p_i - p_j)|²
// where j ranges over (active∪static)∖{i}.
// ∇W_spiky vector at a pair: spiky_const · (h - r)² · (dir / r)
// One block per active particle; 256 threads cooperate.
// ──────────────────────────────────────────────────────────────────────
__global__ void density_constraint_grad(const float3 *active,
                                        const float3 *static_p,
                                        const float *r2_aa,
                                        const float *r2_as,
                                        float3 *grad_C,
                                        float *denom_helper,
                                        float h,
                                        float spiky_const,
                                        float mass,
                                        float rho_rest,
                                        unsigned int n_active,
                                        unsigned int n_static)
{
    unsigned int i = blockIdx.x;
    if (i >= n_active) return;

    __shared__ float3 partials_grad[256];
    __shared__ float  partials_denom[256];

    unsigned int t = threadIdx.x;
    unsigned int T = blockDim.x;

    float h2 = h * h;
    float3 p_i = active[i];
    float3 partial_grad = make_float3(0.0f, 0.0f, 0.0f);
    float partial_denom = 0.0f;
    unsigned int n_total = n_active + n_static;

    for (unsigned int k = t; k < n_total; k += T) {
        float r2;
        float3 dir;
        if (k < n_active) {
            unsigned int j = k;
            if (j == i) continue;
            r2 = r2_aa[i * n_active + j];
            if (r2 >= h2) continue;
            float3 pj = active[j];
            dir = make_float3(p_i.x - pj.x, p_i.y - pj.y, p_i.z - pj.z);
        } else {
            unsigned int j = k - n_active;
            r2 = r2_as[i * n_static + j];
            if (r2 >= h2) continue;
            float3 pj = static_p[j];
            dir = make_float3(p_i.x - pj.x, p_i.y - pj.y, p_i.z - pj.z);
        }
        float r = sqrtf(r2);
        float h_minus_r = h - r;
        float coef = spiky_const * h_minus_r * h_minus_r / (r + 1e-7f);
        float3 grad_W = make_float3(coef * dir.x, coef * dir.y, coef * dir.z);
        partial_grad.x += mass * grad_W.x;
        partial_grad.y += mass * grad_W.y;
        partial_grad.z += mass * grad_W.z;
        partial_denom += grad_W.x*grad_W.x + grad_W.y*grad_W.y
                       + grad_W.z*grad_W.z;
    }

    partials_grad[t]  = partial_grad;
    partials_denom[t] = partial_denom;
    __syncthreads();

    for (unsigned int stride = T / 2; stride > 0; stride >>= 1) {
        if (t < stride) {
            partials_grad[t].x  += partials_grad[t + stride].x;
            partials_grad[t].y  += partials_grad[t + stride].y;
            partials_grad[t].z  += partials_grad[t + stride].z;
            partials_denom[t]   += partials_denom[t + stride];
        }
        __syncthreads();
    }
    if (t == 0) {
        float inv_rho = 1.0f / rho_rest;
        grad_C[i] = make_float3(partials_grad[0].x * inv_rho,
                                partials_grad[0].y * inv_rho,
                                partials_grad[0].z * inv_rho);
        denom_helper[i] = partials_denom[0];
    }
}
// M7.0 — predict_positions: gravity + ballistic integration to x_pred.
__global__ void predict_positions(const float3 *pos_old,
                                  const float3 *vel,
                                  float3 *pos_pred,
                                  float dt,
                                  float gravity_y,
                                  unsigned int n,
                                  float sim_scale_inv)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n) return;
    float3 x = pos_old[gid];
    float3 v = vel[gid];
    float vy_pred = v.y + gravity_y * dt;
    pos_pred[gid] = make_float3(
        x.x + v.x * dt * sim_scale_inv,
        x.y + vy_pred * dt * sim_scale_inv,
        x.z + v.z * dt * sim_scale_inv);
}

// M7.1 — solve_density_constraint: project one-sided overcompression.
__global__ void solve_density_constraint(float3 *pos_pred,
                                         float *lambda,
                                         const float *density,
                                         const float3 *grad_C,
                                         const float *denom_helper,
                                         float rho_rest,
                                         float mass,
                                         float alpha_inv_dt2,
                                         unsigned int n_active)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n_active) return;
    float C = density[gid] / rho_rest - 1.0f;
    if (C <= 0.0f) return;
    float3 g = grad_C[gid];
    float g2 = g.x*g.x + g.y*g.y + g.z*g.z;
    float helper = denom_helper[gid];
    float denom = g2 / mass + (mass / (rho_rest * rho_rest)) * helper
                 + alpha_inv_dt2;
    if (denom < 1e-9f) return;
    float dlambda = -(C + alpha_inv_dt2 * lambda[gid]) / denom;
    float3 p = pos_pred[gid];
    pos_pred[gid] = make_float3(p.x + g.x * dlambda / mass,
                                p.y + g.y * dlambda / mass,
                                p.z + g.z * dlambda / mass);
    lambda[gid] += dlambda;
}

// M7.2 — solve_distance_constraints_seq: Gauss-Seidel bond projection.
// Single-threaded inside the GPU (one block, one thread) — bonds are
// processed sequentially so each sees the latest particle positions.
// For ~150 bonds × 25 ns/bond = 4 µs/iter — negligible.
__global__ void solve_distance_constraints_seq(float3 *pos_pred,
                                               float *lambdas,
                                               const int2 *bond_ij,
                                               const float *rest_len,
                                               float alpha_inv_dt2,
                                               float mass_inv,
                                               unsigned int n_bonds)
{
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    for (unsigned int b = 0; b < n_bonds; b++) {
        int i = bond_ij[b].x;
        int j = bond_ij[b].y;
        float3 pi = pos_pred[i];
        float3 pj = pos_pred[j];
        float dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
        float d = sqrtf(dx*dx + dy*dy + dz*dz);
        if (d < 1e-7f) continue;
        float C = d - rest_len[b];
        float inv_d = 1.0f / d;
        float gx = dx * inv_d, gy = dy * inv_d, gz = dz * inv_d;
        float denom = 2.0f * mass_inv + alpha_inv_dt2;
        float dlambda = -(C + alpha_inv_dt2 * lambdas[b]) / denom;
        float corr = dlambda * mass_inv;
        pos_pred[i] = make_float3(pi.x + gx * corr,
                                  pi.y + gy * corr,
                                  pi.z + gz * corr);
        pos_pred[j] = make_float3(pj.x - gx * corr,
                                  pj.y - gy * corr,
                                  pj.z - gz * corr);
        lambdas[b] += dlambda;
    }
}

// M7.3 — solve_floor_constraint: reflect with optional restitution.
__global__ void solve_floor_constraint(float3 *pos_pred,
                                       float floor_y,
                                       unsigned int n,
                                       float restitution)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n) return;
    float3 x = pos_pred[gid];
    if (x.y < floor_y) {
        float delta = floor_y - x.y;
        x.y = floor_y + restitution * delta;
        pos_pred[gid] = x;
    }
}

// M7.4 — update_velocities: v = (x_pred - x_old) · sim_scale / dt.
__global__ void update_velocities(float3 *vel,
                                  const float3 *pos_old,
                                  const float3 *pos_pred,
                                  float dt,
                                  unsigned int n,
                                  float sim_scale)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n) return;
    float3 a = pos_old[gid];
    float3 b = pos_pred[gid];
    float inv_dt = 1.0f / dt;
    vel[gid] = make_float3((b.x - a.x) * sim_scale * inv_dt,
                           (b.y - a.y) * sim_scale * inv_dt,
                           (b.z - a.z) * sim_scale * inv_dt);
}

// Utility — element-wise add: dst += src.
__global__ void add_inplace(float *dst, const float *src, unsigned int n)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n) return;
    dst[gid] += src[gid];
}

// Safety belt — clamp |v| ≤ v_max (CFL stability against teleport).
__global__ void clamp_velocity(float3 *vel, float v_max, unsigned int n)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n) return;
    float3 v = vel[gid];
    float speed = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
    if (speed > v_max && speed > 0.0f) {
        float s = v_max / speed;
        vel[gid] = make_float3(v.x * s, v.y * s, v.z * s);
    }
}

// Utility — sum a 1-D float array and atomicAdd the reduction to a scalar.
// Single-block launch <<<1, TPB>>>; TPB must equal the partials[] size below.
// Used to accumulate per-particle scalar grads (rho_rest, alpha_density)
// from the M7.1_bwd output into a running scalar accumulator across the
// multi-step reverse loop.
__global__ void sum_atomic_to_scalar(const float *src, float *dst,
                                     unsigned int n)
{
    __shared__ float partials[256];
    unsigned int t = threadIdx.x;
    unsigned int T = blockDim.x;
    float sum = 0.0f;
    for (unsigned int k = t; k < n; k += T) sum += src[k];
    partials[t] = sum;
    __syncthreads();
    for (unsigned int s = T / 2; s > 0; s >>= 1) {
        if (t < s) partials[t] += partials[t + s];
        __syncthreads();
    }
    if (t == 0) atomicAdd(dst, partials[0]);
}

// Adds the M6.4-side contribution to ∂L/∂ρ_rest.
// M6.4 forward computes grad_C[i] = (m/ρ_rest)·Σ_k ∇W, so grad_C scales as
// 1/ρ_rest. Given ω_i = ∂L/∂grad_C[i] (= grad_grad_C from M7.1_bwd) and
// the saved grad_C values, the correction is:
//   ∂L/∂ρ_rest_via_grad_C += -(1/ρ_rest) · Σ_i ⟨ω_i, grad_C[i]⟩
// This complements the direct ∂L/∂ρ_rest contribution from M7.1_bwd
// (which only captures C and D dependencies in solve_density_constraint).
// Single-block dot-product reduction → atomicAdd to scalar accumulator.
__global__ void rho_rest_grad_via_M6_4(const float3 *grad_grad_C,
                                       const float3 *grad_C_saved,
                                       float *grad_rho_dst,
                                       float rho_rest,
                                       unsigned int n)
{
    __shared__ float partials[256];
    unsigned int t = threadIdx.x;
    unsigned int T = blockDim.x;
    float sum = 0.0f;
    for (unsigned int i = t; i < n; i += T) {
        float3 a = grad_grad_C[i];
        float3 b = grad_C_saved[i];
        sum += a.x*b.x + a.y*b.y + a.z*b.z;
    }
    partials[t] = sum;
    __syncthreads();
    for (unsigned int s = T / 2; s > 0; s >>= 1) {
        if (t < s) partials[t] += partials[t + s];
        __syncthreads();
    }
    if (t == 0) atomicAdd(grad_rho_dst, -partials[0] / rho_rest);
}

// ──────────────────────────────────────────────────────────────────────
// M7 backward kernels — paired adjoint of predict/floor/update.
//
// Phase 4 demonstrates the differentiable architecture on the simplest
// non-trivial chain: predict_positions -> solve_floor_constraint ->
// update_velocities. Density + bonds backward kernels follow once their
// forward analog is exercised under SGD.
// ──────────────────────────────────────────────────────────────────────

// predict_positions_backward:
//   x_pred.x = x_old.x + v.x * dt * sim_scale_inv
//   x_pred.y = x_old.y + (v.y + g*dt) * dt * sim_scale_inv
//   x_pred.z = x_old.z + v.z * dt * sim_scale_inv
//
// Given grad_pos_pred:
//   grad_pos_old += grad_pos_pred                            (identity)
//   grad_vel.{x,y,z} = grad_pos_pred.{x,y,z} * dt * si
//   grad_gravity_y  += grad_pos_pred.y * dt^2 * si           (reduced over i)
__global__ void predict_positions_backward(const float3 *grad_pos_pred,
                                           float3 *grad_pos_old,
                                           float3 *grad_vel,
                                           float *grad_gravity_y,
                                           float dt,
                                           float sim_scale_inv,
                                           unsigned int n)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n) return;
    float3 g_xp = grad_pos_pred[gid];
    // grad_pos_old += grad_pos_pred (identity).
    float3 prev = grad_pos_old[gid];
    grad_pos_old[gid] = make_float3(prev.x + g_xp.x,
                                    prev.y + g_xp.y,
                                    prev.z + g_xp.z);
    // grad_vel = grad_pos_pred * dt * sim_scale_inv.
    float k = dt * sim_scale_inv;
    grad_vel[gid] = make_float3(g_xp.x * k, g_xp.y * k, g_xp.z * k);
    // grad_gravity_y += grad_pos_pred.y * dt * dt * sim_scale_inv (atomic).
    float contrib = g_xp.y * dt * dt * sim_scale_inv;
    atomicAdd(grad_gravity_y, contrib);
}

// solve_floor_constraint_backward:
//   m = (x_pred_saved.y < floor_y)
//   x_post.y = m ? (1+r)*floor_y - r*x_pred.y : x_pred.y
//   x_post.x = x_pred.x
//   x_post.z = x_pred.z
//
// Given grad_pos_post:
//   grad_pos_pred.y = m ? grad_pos_post.y * (-r) : grad_pos_post.y
//   grad_pos_pred.{x,z} = grad_pos_post.{x,z}
//   grad_floor_y += m ? grad_pos_post.y * (1+r) : 0
//   grad_restitution += m ? grad_pos_post.y * (floor_y - x_pred.y) : 0
__global__ void solve_floor_constraint_backward(const float3 *grad_pos_post,
                                                const float3 *pos_pred_saved,
                                                float3 *grad_pos_pred,
                                                float *grad_floor_y,
                                                float *grad_restitution,
                                                float floor_y,
                                                float restitution,
                                                unsigned int n)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n) return;
    float3 gpp = grad_pos_post[gid];
    float3 xp = pos_pred_saved[gid];
    bool m = xp.y < floor_y;
    float gy = m ? (-restitution) * gpp.y : gpp.y;
    grad_pos_pred[gid] = make_float3(gpp.x, gy, gpp.z);
    if (m) {
        atomicAdd(grad_floor_y, gpp.y * (1.0f + restitution));
        atomicAdd(grad_restitution, gpp.y * (floor_y - xp.y));
    }
}

// solve_distance_constraints_seq_with_save: same as the forward bond
// projector but also persists per-bond pre-state (pi, pj, λ_pre) to a
// `state` buffer of 7·n_bonds floats. The backward kernel needs these
// to reconstruct Δλ and the projection geometry.
__global__ void solve_distance_constraints_seq_with_save(
    float3 *pos_pred,
    float *lambdas,
    const int2 *bond_ij,
    const float *rest_len,
    float *state,
    float alpha_inv_dt2,
    float mass_inv,
    unsigned int n_bonds)
{
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    for (unsigned int b = 0; b < n_bonds; b++) {
        int i = bond_ij[b].x;
        int j = bond_ij[b].y;
        float3 pi = pos_pred[i];
        float3 pj = pos_pred[j];
        float lambda_pre = lambdas[b];

        unsigned int base = b * 7;
        state[base + 0] = pi.x; state[base + 1] = pi.y; state[base + 2] = pi.z;
        state[base + 3] = pj.x; state[base + 4] = pj.y; state[base + 5] = pj.z;
        state[base + 6] = lambda_pre;

        float dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
        float d = sqrtf(dx*dx + dy*dy + dz*dz);
        if (d < 1e-7f) continue;
        float C = d - rest_len[b];
        float inv_d = 1.0f / d;
        float gx = dx * inv_d, gy = dy * inv_d, gz = dz * inv_d;
        float D = 2.0f * mass_inv + alpha_inv_dt2;
        float dlambda = -(C + alpha_inv_dt2 * lambda_pre) / D;
        float corr = dlambda * mass_inv;
        pos_pred[i] = make_float3(pi.x + gx * corr,
                                  pi.y + gy * corr,
                                  pi.z + gz * corr);
        pos_pred[j] = make_float3(pj.x - gx * corr,
                                  pj.y - gy * corr,
                                  pj.z - gz * corr);
        lambdas[b] = lambda_pre + dlambda;
    }
}

// solve_distance_constraints_seq_backward: hand-derived adjoint of the
// Gauss-Seidel bond projection. Walks bonds in REVERSE order; each bond's
// backward sees pos_grad as updated by later bonds. Single thread; no
// races.
//
// alpha_param is the compliance α (NOT α/dt²); dt2 = dt²; the kernel
// scales internally per the chain-rule derivation in shaders.metal.
__global__ void solve_distance_constraints_seq_backward(
    float3 *pos_grad,
    float *lambda_grad,
    float *alpha_grad,
    const int2 *bond_ij,
    const float *rest_len,
    const float *state,
    float alpha_inv_dt2,
    float alpha_param,
    float dt2,
    float mass_inv,
    unsigned int n_bonds)
{
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    float local_alpha_grad = 0.0f;
    for (int bi = (int)n_bonds - 1; bi >= 0; bi--) {
        unsigned int b = (unsigned int)bi;
        int i = bond_ij[b].x;
        int j = bond_ij[b].y;
        unsigned int base = b * 7;

        float3 pi_pre = make_float3(state[base+0], state[base+1], state[base+2]);
        float3 pj_pre = make_float3(state[base+3], state[base+4], state[base+5]);
        float lambda_pre = state[base+6];

        float vx = pi_pre.x - pj_pre.x;
        float vy = pi_pre.y - pj_pre.y;
        float vz = pi_pre.z - pj_pre.z;
        float d = sqrtf(vx*vx + vy*vy + vz*vz);
        if (d < 1e-7f) continue;
        float inv_d = 1.0f / d;
        float gx = vx * inv_d, gy = vy * inv_d, gz = vz * inv_d;
        float L = rest_len[b];
        float C = d - L;
        float D = 2.0f * mass_inv + alpha_inv_dt2;
        float dlambda = -(C + alpha_inv_dt2 * lambda_pre) / D;

        float3 omega = pos_grad[i];
        float3 phi   = pos_grad[j];
        float psi    = lambda_grad[b];

        float deltax = omega.x - phi.x;
        float deltay = omega.y - phi.y;
        float deltaz = omega.z - phi.z;
        float delta_g = deltax * gx + deltay * gy + deltaz * gz;

        float coef1 = dlambda * mass_inv * inv_d;
        float coef2 = -delta_g * mass_inv / D;
        // delta_J = coef1 · (delta - delta_g · g) + coef2 · g
        float dJx = coef1 * (deltax - delta_g * gx) + coef2 * gx;
        float dJy = coef1 * (deltay - delta_g * gy) + coef2 * gy;
        float dJz = coef1 * (deltaz - delta_g * gz) + coef2 * gz;

        float psi_over_D = psi / D;
        pos_grad[i] = make_float3(omega.x + dJx - psi_over_D * gx,
                                  omega.y + dJy - psi_over_D * gy,
                                  omega.z + dJz - psi_over_D * gz);
        pos_grad[j] = make_float3(phi.x - dJx + psi_over_D * gx,
                                  phi.y - dJy + psi_over_D * gy,
                                  phi.z - dJz + psi_over_D * gz);

        lambda_grad[b] = -delta_g * alpha_param * mass_inv / (dt2 * D)
                       + psi * 2.0f * mass_inv / D;

        float dlambda_dalpha = (C - 2.0f * lambda_pre * mass_inv) / (dt2 * D * D);
        float chain_param = delta_g * mass_inv + psi;
        local_alpha_grad += chain_param * dlambda_dalpha;
    }
    atomicAdd(alpha_grad, local_alpha_grad);
}

// update_velocities_backward:
//   v_new[i] = (x_post[i] - x_old[i]) * sim_scale / dt
//
// Given grad_vel_new:
//   grad_pos_post[i] += grad_vel_new[i] * sim_scale / dt
//   grad_pos_old[i]  += -grad_vel_new[i] * sim_scale / dt
__global__ void update_velocities_backward(const float3 *grad_vel_new,
                                           float3 *grad_pos_post,
                                           float3 *grad_pos_old,
                                           float dt,
                                           float sim_scale,
                                           unsigned int n)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n) return;
    float3 gvn = grad_vel_new[gid];
    float k = sim_scale / dt;
    float3 prev_post = grad_pos_post[gid];
    grad_pos_post[gid] = make_float3(prev_post.x + gvn.x * k,
                                     prev_post.y + gvn.y * k,
                                     prev_post.z + gvn.z * k);
    float3 prev_old = grad_pos_old[gid];
    grad_pos_old[gid] = make_float3(prev_old.x - gvn.x * k,
                                    prev_old.y - gvn.y * k,
                                    prev_old.z - gvn.z * k);
}

// ══════════════════════════════════════════════════════════════════════
// M9: density-chain backward kernels (gradient chain Option 3).
//   dist_active_static_bwd       (M6.0_bwd)
//   dist_active_active_bwd       (M6.3_bwd)
//   wpoly6_inplace_bwd           (M6.1_bwd, in-place ∂L/∂W → ∂L/∂r²)
//   rowsum_density_bwd           (M6.2_bwd, broadcast)
//   density_constraint_grad_bwd  (M6.4_bwd a.k.a. M9.C)
//   solve_density_constraint_bwd (M7.1_bwd a.k.a. M9.D)
// Direct ports of src/metal_diff/shaders.metal M9 section. Derivations
// match the Metal comments; see those for math details.
// ══════════════════════════════════════════════════════════════════════

// M6.0_bwd — backward of dist_active_static.
// Static positions are frozen (no gradient written for them).
// ∂L/∂active[i] += Σ_j ∂L/∂r²[i,j] · 2·(active[i] - static_p[j])
// One thread per active i; no race (each thread owns its row of grad_active).
__global__ void dist_active_static_bwd(const float3 *active,
                                       const float3 *static_p,
                                       const float *grad_r2,
                                       float3 *grad_active,
                                       unsigned int n_active,
                                       unsigned int n_static)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_active) return;
    float3 ai = active[i];
    float gx = 0.0f, gy = 0.0f, gz = 0.0f;
    for (unsigned int j = 0; j < n_static; j++) {
        float3 sj = static_p[j];
        float c = 2.0f * grad_r2[i * n_static + j];
        gx += c * (ai.x - sj.x);
        gy += c * (ai.y - sj.y);
        gz += c * (ai.z - sj.z);
    }
    float3 prev = grad_active[i];
    grad_active[i] = make_float3(prev.x + gx, prev.y + gy, prev.z + gz);
}
//   ∂L/∂active[i] += Σ_{j≠i} (grad_r2[i,j] + grad_r2[j,i]) · 2·(active[i] - active[j])
__global__ void dist_active_active_bwd(const float3 *active,
                                       const float *grad_r2,
                                       float3 *grad_active,
                                       unsigned int n_active)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_active) return;
    float3 ai = active[i];
    float gx = 0.0f, gy = 0.0f, gz = 0.0f;
    for (unsigned int j = 0; j < n_active; j++) {
        if (j == i) continue;
        float3 aj = active[j];
        float c = 2.0f * (grad_r2[i * n_active + j]
                       + grad_r2[j * n_active + i]);
        gx += c * (ai.x - aj.x);
        gy += c * (ai.y - aj.y);
        gz += c * (ai.z - aj.z);
    }
    float3 prev = grad_active[i];
    grad_active[i] = make_float3(prev.x + gx, prev.y + gy, prev.z + gz);
}
// Requires saved r² from the forward pass (separate buffer).
__global__ void wpoly6_inplace_bwd(const float *r2_saved,
                                   float *grad_W_or_r2,
                                   float h2,
                                   float poly6_const,
                                   unsigned int n_total)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n_total) return;
    float r2 = r2_saved[gid];
    if (r2 >= h2) { grad_W_or_r2[gid] = 0.0f; return; }
    float diff = h2 - r2;
    float dW_dr2 = -3.0f * poly6_const * diff * diff;
    grad_W_or_r2[gid] = grad_W_or_r2[gid] * dW_dr2;
}

// 2-D dispatch (n_rows × n_cols), each thread writes one cell. Overwrites.
__global__ void rowsum_density_bwd(const float *grad_density,
                                   float *grad_W,
                                   float mass,
                                   unsigned int n_rows,
                                   unsigned int n_cols)
{
    unsigned int j = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int i = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= n_rows || j >= n_cols) return;
    grad_W[i * n_cols + j] = mass * grad_density[i];
}
// One thread per active i; no race (each thread writes only grad_active[i]).
__global__ void density_constraint_grad_bwd(const float3 *active,
                                            const float3 *static_p,
                                            const float *r2_aa,
                                            const float *r2_as,
                                            const float3 *grad_grad_C,   // ω
                                            const float *grad_denom_h,    // ψ
                                            float3 *grad_active,
                                            float h,
                                            float spiky_const,
                                            float mass,
                                            float rho_rest,
                                            unsigned int n_active,
                                            unsigned int n_static)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_active) return;
    const float R_MIN = 1e-4f;
    const float R_MIN2 = R_MIN * R_MIN;
    float h2 = h * h;
    float scale = mass / rho_rest;
    float3 p_i = active[i];
    float3 omega_i = grad_grad_C[i];
    float  psi_i   = grad_denom_h[i];
    float gx = 0.0f, gy = 0.0f, gz = 0.0f;

    // Active neighbors: both "self" (row i) and "as-neighbor" (row j) contributions.
    for (unsigned int j = 0; j < n_active; j++) {
        if (j == i) continue;
        float r2 = r2_aa[i * n_active + j];
        if (r2 >= h2) continue;
        if (r2 < R_MIN2) continue;
        float r = sqrtf(r2);
        float r_safe = r > R_MIN ? r : R_MIN;
        float h_minus_r = h - r;
        float3 pj = active[j];
        float vx = p_i.x - pj.x, vy = p_i.y - pj.y, vz = p_i.z - pj.z;
        float inv_r_safe = 1.0f / r_safe;
        float ghx = vx * inv_r_safe, ghy = vy * inv_r_safe, ghz = vz * inv_r_safe;
        float coef = spiky_const * h_minus_r * h_minus_r / (r + 1e-7f);
        float gWx = coef * vx, gWy = coef * vy, gWz = coef * vz;

        // self (row i): u_self = scale·ω_i + 2·ψ_i·grad_W
        float ux = scale * omega_i.x + 2.0f * psi_i * gWx;
        float uy = scale * omega_i.y + 2.0f * psi_i * gWy;
        float uz = scale * omega_i.z + 2.0f * psi_i * gWz;
        float ug = ux * ghx + uy * ghy + uz * ghz;
        float upx = ux - ug * ghx, upy = uy - ug * ghy, upz = uz - ug * ghz;
        float a1 = spiky_const * h_minus_r * (h_minus_r * inv_r_safe);
        float a2 = -2.0f * spiky_const * h_minus_r * ug;
        gx += a1 * upx + a2 * ghx;
        gy += a1 * upy + a2 * ghy;
        gz += a1 * upz + a2 * ghz;

        // neighbor (row j, pair (j,i)): u_neigh = scale·ω_j - 2·ψ_j·grad_W
        float3 omega_j = grad_grad_C[j];
        float  psi_j   = grad_denom_h[j];
        float unx = scale * omega_j.x - 2.0f * psi_j * gWx;
        float uny = scale * omega_j.y - 2.0f * psi_j * gWy;
        float unz = scale * omega_j.z - 2.0f * psi_j * gWz;
        float ugn = unx * ghx + uny * ghy + unz * ghz;
        float upnx = unx - ugn * ghx, upny = uny - ugn * ghy, upnz = unz - ugn * ghz;
        float b1 = spiky_const * h_minus_r * (h_minus_r * inv_r_safe);
        float b2 = -2.0f * spiky_const * h_minus_r * ugn;
        // ∂grad_W_ji/∂p_i = -J ⟹ subtract the neighbor's J·u term
        gx -= b1 * upnx + b2 * ghx;
        gy -= b1 * upny + b2 * ghy;
        gz -= b1 * upnz + b2 * ghz;
    }

    // Static neighbors: only "self" contribution (no gradient flows to frozen static).
    for (unsigned int k = 0; k < n_static; k++) {
        float r2 = r2_as[i * n_static + k];
        if (r2 >= h2) continue;
        if (r2 < R_MIN2) continue;
        float r = sqrtf(r2);
        float r_safe = r > R_MIN ? r : R_MIN;
        float h_minus_r = h - r;
        float3 sk = static_p[k];
        float vx = p_i.x - sk.x, vy = p_i.y - sk.y, vz = p_i.z - sk.z;
        float inv_r_safe = 1.0f / r_safe;
        float ghx = vx * inv_r_safe, ghy = vy * inv_r_safe, ghz = vz * inv_r_safe;
        float coef = spiky_const * h_minus_r * h_minus_r / (r + 1e-7f);
        float gWx = coef * vx, gWy = coef * vy, gWz = coef * vz;

        float ux = scale * omega_i.x + 2.0f * psi_i * gWx;
        float uy = scale * omega_i.y + 2.0f * psi_i * gWy;
        float uz = scale * omega_i.z + 2.0f * psi_i * gWz;
        float ug = ux * ghx + uy * ghy + uz * ghz;
        float upx = ux - ug * ghx, upy = uy - ug * ghy, upz = uz - ug * ghz;
        float a1 = spiky_const * h_minus_r * (h_minus_r * inv_r_safe);
        float a2 = -2.0f * spiky_const * h_minus_r * ug;
        gx += a1 * upx + a2 * ghx;
        gy += a1 * upy + a2 * ghy;
        gz += a1 * upz + a2 * ghz;
    }

    float3 prev = grad_active[i];
    grad_active[i] = make_float3(prev.x + gx, prev.y + gy, prev.z + gz);
}
// through pos_pre with zeros for all other gradients.
__global__ void solve_density_constraint_bwd(const float *density,
                                             const float3 *grad_C,
                                             const float *denom_helper,
                                             const float *lambda_pre,
                                             const float3 *grad_pos_post,
                                             const float *grad_lambda_post,
                                             float3 *grad_pos_pre,
                                             float *grad_lambda_pre,
                                             float *grad_density,
                                             float3 *grad_grad_C,
                                             float *grad_denom_h,
                                             float *grad_rho_rest,
                                             float *grad_alpha,
                                             float rho_rest,
                                             float mass,
                                             float alpha_inv_dt2,
                                             unsigned int n_active)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_active) return;

    float C = density[i] / rho_rest - 1.0f;
    float3 omega = grad_pos_post[i];
    float psi = grad_lambda_post[i];

    if (C <= 0.0f) {
        float3 prev = grad_pos_pre[i];
        grad_pos_pre[i] = make_float3(prev.x + omega.x, prev.y + omega.y, prev.z + omega.z);
        grad_lambda_pre[i] = psi;
        grad_density[i] = 0.0f;
        grad_grad_C[i] = make_float3(0.0f, 0.0f, 0.0f);
        grad_denom_h[i] = 0.0f;
        grad_rho_rest[i] = 0.0f;
        grad_alpha[i] = 0.0f;
        return;
    }

    float3 g = grad_C[i];
    float g2 = g.x*g.x + g.y*g.y + g.z*g.z;
    float helper = denom_helper[i];
    float lam = lambda_pre[i];
    float D = g2 / mass + (mass / (rho_rest * rho_rest)) * helper + alpha_inv_dt2;
    float dlambda = -(C + alpha_inv_dt2 * lam) / D;
    float lambda_post = lam + dlambda;

    float omega_dot_g = omega.x*g.x + omega.y*g.y + omega.z*g.z;
    float chain = omega_dot_g / mass + psi;

    float3 prev = grad_pos_pre[i];
    grad_pos_pre[i] = make_float3(prev.x + omega.x, prev.y + omega.y, prev.z + omega.z);

    grad_lambda_pre[i] = -omega_dot_g * alpha_inv_dt2 / (mass * D)
                       + psi * (1.0f - alpha_inv_dt2 / D);

    grad_density[i] = -chain / (D * rho_rest);

    float a_g = dlambda / mass;
    float b_g = -2.0f * dlambda / (mass * D) * chain;
    grad_grad_C[i] = make_float3(a_g * omega.x + b_g * g.x,
                                 a_g * omega.y + b_g * g.y,
                                 a_g * omega.z + b_g * g.z);

    grad_denom_h[i] = -dlambda * mass / (rho_rest * rho_rest * D) * chain;

    float dlambda_dr = density[i] / (rho_rest * rho_rest * D)
                     + 2.0f * dlambda * mass * helper
                       / (rho_rest * rho_rest * rho_rest * D);
    grad_rho_rest[i] = chain * dlambda_dr;

    grad_alpha[i] = chain * (-lambda_post / D);
}

// ══════════════════════════════════════════════════════════════════════
// Worm physics: force-based springs, anchors, and external acceleration.
// Direct ports of src/metal_diff/shaders.metal lines 618-773 + 942-1072.
//
// Unlike the constraint-based XPBD bonds (solve_distance_constraints_seq),
// these springs apply Hooke forces directly to an `ext_accel` buffer that
// gets fed into apply_ext_accel between density-projection and predict.
// Each per-particle thread loops over n_bonds to find its incident bonds
// (O(n·n_bonds), fine for the bond counts we use).
// ══════════════════════════════════════════════════════════════════════

// spring_bonds_force — Hooke spring force on each elastic→elastic bond.
// For each bond (i, j) with rest length L, the force on i is
//   a_i += -K · (||p_i - p_j|| - L) · (p_i - p_j) / ||p_i - p_j||
// (and the opposite sign for j; we let thread j compute its own contribution).
__global__ void spring_bonds_force(const float3 *active_pos,
                                   const int2 *bond_ij,
                                   const float *bond_rest,
                                   float3 *ext_accel,
                                   float spring_K,
                                   unsigned int n_bonds,
                                   unsigned int n_active)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_active) return;
    float3 p_i = active_pos[i];
    float ax = 0.0f, ay = 0.0f, az = 0.0f;
    for (unsigned int b = 0; b < n_bonds; b++) {
        int i_b = bond_ij[b].x;
        int j_b = bond_ij[b].y;
        int other;
        if (i_b == (int)i)      other = j_b;
        else if (j_b == (int)i) other = i_b;
        else continue;
        float L = bond_rest[b];
        float3 p_o = active_pos[other];
        float dx = p_i.x - p_o.x, dy = p_i.y - p_o.y, dz = p_i.z - p_o.z;
        float r = sqrtf(dx*dx + dy*dy + dz*dz);
        if (r < 1e-7f) continue;
        float coef = -spring_K * (r - L) / r;
        ax += coef * dx; ay += coef * dy; az += coef * dz;
    }
    float3 prev = ext_accel[i];
    ext_accel[i] = make_float3(prev.x + ax, prev.y + ay, prev.z + az);
}

// apply_ext_accel — v += dt · ext_accel
__global__ void apply_ext_accel(float3 *vel, const float3 *ext_accel,
                                float dt, unsigned int n)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n) return;
    float3 v = vel[gid];
    float3 a = ext_accel[gid];
    vel[gid] = make_float3(v.x + dt * a.x,
                           v.y + dt * a.y,
                           v.z + dt * a.z);
}

// spring_bonds_force_backward — gradient of ext_accel w.r.t. positions.
// Per bond (i, j), with G_diff = ga_j - ga_i:
//   ∂L/∂p_i += K · [(1 - L/r) · G_diff + (L/r³) · dir · ⟨dir, G_diff⟩]
// Each thread for particle i loops over its incident bonds and accumulates
// into its own grad_pos[i]; thread j picks up the symmetric contribution.
__global__ void spring_bonds_force_backward(const float3 *active_pos,
                                            const int2 *bond_ij,
                                            const float *bond_rest,
                                            const float3 *grad_ext_accel,
                                            float3 *grad_pos,
                                            float spring_K,
                                            unsigned int n_bonds,
                                            unsigned int n_active)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_active) return;
    float3 p_i = active_pos[i];
    float3 ga_i = grad_ext_accel[i];
    float dpx = 0.0f, dpy = 0.0f, dpz = 0.0f;
    for (unsigned int b = 0; b < n_bonds; b++) {
        int i_b = bond_ij[b].x;
        int j_b = bond_ij[b].y;
        int other;
        if (i_b == (int)i)      other = j_b;
        else if (j_b == (int)i) other = i_b;
        else continue;
        float L = bond_rest[b];
        float3 p_o = active_pos[other];
        float dx = p_i.x - p_o.x, dy = p_i.y - p_o.y, dz = p_i.z - p_o.z;
        float r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-14f) continue;
        float r = sqrtf(r2);
        float r3 = r * r2;
        float3 ga_o = grad_ext_accel[other];
        float Gx = ga_o.x - ga_i.x, Gy = ga_o.y - ga_i.y, Gz = ga_o.z - ga_i.z;
        float dot_dir_G = dx*Gx + dy*Gy + dz*Gz;
        float coef_iso  = 1.0f - L / r;     // identity-projection coefficient
        float coef_proj = L / r3;           // outer-product coefficient
        dpx += spring_K * (coef_iso * Gx + coef_proj * dx * dot_dir_G);
        dpy += spring_K * (coef_iso * Gy + coef_proj * dy * dot_dir_G);
        dpz += spring_K * (coef_iso * Gz + coef_proj * dz * dot_dir_G);
    }
    float3 prev = grad_pos[i];
    grad_pos[i] = make_float3(prev.x + dpx, prev.y + dpy, prev.z + dpz);
}

// spring_K_partial — per-particle partial of ∂L/∂spring_K. The per-bond
// contribution to a_i is -K·(r-L)·dir/r; differentiating w.r.t. K gives
// daK_i = -(r-L)·dir/r. The per-particle partial = ⟨daK_i, ga_i⟩, summed
// over incident bonds. Host sums per_particle[] to a scalar via
// sum_atomic_to_scalar to get ∂L/∂K.
__global__ void spring_K_partial(const float3 *active_pos,
                                 const int2 *bond_ij,
                                 const float *bond_rest,
                                 const float3 *grad_ext_accel,
                                 float *per_particle,
                                 unsigned int n_bonds,
                                 unsigned int n_active)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_active) return;
    float3 p_i = active_pos[i];
    float3 ga_i = grad_ext_accel[i];
    float partial = 0.0f;
    for (unsigned int b = 0; b < n_bonds; b++) {
        int i_b = bond_ij[b].x;
        int j_b = bond_ij[b].y;
        int other;
        if (i_b == (int)i)      other = j_b;
        else if (j_b == (int)i) other = i_b;
        else continue;
        float L = bond_rest[b];
        float3 p_o = active_pos[other];
        float dx = p_i.x - p_o.x, dy = p_i.y - p_o.y, dz = p_i.z - p_o.z;
        float r = sqrtf(dx*dx + dy*dy + dz*dz);
        if (r < 1e-7f) continue;
        float coef = -(r - L) / r;
        partial += coef * (dx * ga_i.x + dy * ga_i.y + dz * ga_i.z);
    }
    per_particle[i] = partial;
}

// apply_ext_accel_backward — straightforward chain through v += dt·a:
//   grad_v_old += grad_v_new
//   grad_ext_accel += dt · grad_v_new
__global__ void apply_ext_accel_backward(const float3 *grad_v_new,
                                         float3 *grad_v_old,
                                         float3 *grad_ext_accel,
                                         float dt, unsigned int n)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n) return;
    float3 g_new = grad_v_new[gid];
    float3 prev_v = grad_v_old[gid];
    grad_v_old[gid] = make_float3(prev_v.x + g_new.x,
                                  prev_v.y + g_new.y,
                                  prev_v.z + g_new.z);
    float3 prev_a = grad_ext_accel[gid];
    grad_ext_accel[gid] = make_float3(prev_a.x + dt * g_new.x,
                                      prev_a.y + dt * g_new.y,
                                      prev_a.z + dt * g_new.z);
}
