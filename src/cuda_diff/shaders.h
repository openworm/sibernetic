// shaders.h — forward declarations of all __global__ kernels.
//
// Kernels are defined in shaders.cu and launched from driver functions
// in ops_*.cu translation units. Separate compilation requires
// `nvcc -rdc=true` so device-side relocations across TUs are resolved
// at device-link time.
#pragma once

#include <cuda_runtime.h>

// ── M6 forward ────────────────────────────────────────────────────────
__global__ void wpoly6_inplace(float *r2_or_W, float h2, float poly6_const,
                               unsigned int n_total);
__global__ void dist_active_static(const float3 *active, const float3 *static_p,
                                   float *dist, unsigned int n_active,
                                   unsigned int n_static);
__global__ void dist_active_active(const float3 *active, float *dist,
                                   unsigned int n_active);
__global__ void rowsum_density(const float *W, float *density, float mass,
                               unsigned int n_cols, unsigned int n_rows);
__global__ void density_constraint_grad(const float3 *active,
                                        const float3 *static_p,
                                        const float *r2_aa,
                                        const float *r2_as,
                                        float3 *grad_C,
                                        float *denom_helper,
                                        float h, float spiky_const,
                                        float mass, float rho_rest,
                                        unsigned int n_active,
                                        unsigned int n_static);

// ── M7 forward (XPBD step pieces) ─────────────────────────────────────
__global__ void predict_positions(const float3 *pos_old, const float3 *vel,
                                  float3 *pos_pred, float dt, float gravity_y,
                                  unsigned int n, float sim_scale_inv);
__global__ void solve_density_constraint(float3 *pos_pred, float *lambda,
                                         const float *density,
                                         const float3 *grad_C,
                                         const float *denom_helper,
                                         float rho_rest, float mass,
                                         float alpha_inv_dt2,
                                         unsigned int n_active);
__global__ void solve_distance_constraints_seq(float3 *pos_pred, float *lambdas,
                                               const int2 *bond_ij,
                                               const float *rest_len,
                                               float alpha_inv_dt2,
                                               float mass_inv,
                                               unsigned int n_bonds);
__global__ void solve_floor_constraint(float3 *pos_pred, float floor_y,
                                       unsigned int n, float restitution);
__global__ void update_velocities(float3 *vel, const float3 *pos_old,
                                  const float3 *pos_pred, float dt,
                                  unsigned int n, float sim_scale);

// ── Utility ───────────────────────────────────────────────────────────
__global__ void add_inplace(float *dst, const float *src, unsigned int n);
__global__ void clamp_velocity(float3 *vel, float v_max, unsigned int n);
__global__ void sum_atomic_to_scalar(const float *src, float *dst,
                                     unsigned int n);
__global__ void rho_rest_grad_via_M6_4(const float3 *grad_grad_C,
                                       const float3 *grad_C_saved,
                                       float *grad_rho_dst,
                                       float rho_rest, unsigned int n);

// ── M7 backward ───────────────────────────────────────────────────────
__global__ void predict_positions_backward(const float3 *grad_pos_pred,
                                           float3 *grad_pos_old,
                                           float3 *grad_vel,
                                           float *grad_gravity_y,
                                           float dt, float sim_scale_inv,
                                           unsigned int n);
__global__ void solve_floor_constraint_backward(const float3 *grad_pos_post,
                                                const float3 *pos_pred_saved,
                                                float3 *grad_pos_pred,
                                                float *grad_floor_y,
                                                float *grad_restitution,
                                                float floor_y, float restitution,
                                                unsigned int n);
__global__ void solve_distance_constraints_seq_with_save(
    float3 *pos_pred, float *lambdas, const int2 *bond_ij,
    const float *rest_len, float *state,
    float alpha_inv_dt2, float mass_inv, unsigned int n_bonds);
__global__ void solve_distance_constraints_seq_backward(
    float3 *pos_grad, float *lambda_grad, float *alpha_grad,
    const int2 *bond_ij, const float *rest_len, const float *state,
    float alpha_inv_dt2, float alpha_param, float dt2,
    float mass_inv, unsigned int n_bonds);
__global__ void update_velocities_backward(const float3 *grad_vel_new,
                                           float3 *grad_pos_post,
                                           float3 *grad_pos_old,
                                           float dt, float sim_scale,
                                           unsigned int n);

// ── M9 density-chain backward ─────────────────────────────────────────
__global__ void dist_active_static_bwd(const float3 *active,
                                       const float3 *static_p,
                                       const float *grad_r2,
                                       float3 *grad_active,
                                       unsigned int n_active,
                                       unsigned int n_static);
__global__ void dist_active_active_bwd(const float3 *active,
                                       const float *grad_r2,
                                       float3 *grad_active,
                                       unsigned int n_active);
__global__ void wpoly6_inplace_bwd(const float *r2_saved,
                                   float *grad_W_or_r2,
                                   float h2, float poly6_const,
                                   unsigned int n_total);
__global__ void rowsum_density_bwd(const float *grad_density, float *grad_W,
                                   float mass, unsigned int n_rows,
                                   unsigned int n_cols);
__global__ void density_constraint_grad_bwd(const float3 *active,
                                            const float3 *static_p,
                                            const float *r2_aa,
                                            const float *r2_as,
                                            const float3 *grad_grad_C,
                                            const float *grad_denom_h,
                                            float3 *grad_active,
                                            float h, float spiky_const,
                                            float mass, float rho_rest,
                                            unsigned int n_active,
                                            unsigned int n_static);
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
                                             float rho_rest, float mass,
                                             float alpha_inv_dt2,
                                             unsigned int n_active);

// ── Worm physics: force-based springs + anchors + ext accel ──────────
__global__ void spring_bonds_force(const float3 *active_pos,
                                   const int2 *bond_ij,
                                   const float *bond_rest,
                                   float3 *ext_accel,    // accumulate
                                   float spring_K,
                                   unsigned int n_bonds,
                                   unsigned int n_active);
__global__ void apply_ext_accel(float3 *vel, const float3 *ext_accel,
                                float dt, unsigned int n);
__global__ void spring_bonds_force_backward(const float3 *active_pos,
                                            const int2 *bond_ij,
                                            const float *bond_rest,
                                            const float3 *grad_ext_accel,
                                            float3 *grad_pos,     // accumulate
                                            float spring_K,
                                            unsigned int n_bonds,
                                            unsigned int n_active);
__global__ void spring_K_partial(const float3 *active_pos,
                                 const int2 *bond_ij,
                                 const float *bond_rest,
                                 const float3 *grad_ext_accel,
                                 float *per_particle,  // write
                                 unsigned int n_bonds,
                                 unsigned int n_active);
__global__ void apply_ext_accel_backward(const float3 *grad_v_new,
                                         float3 *grad_v_old,
                                         float3 *grad_ext_accel,
                                         float dt, unsigned int n);
