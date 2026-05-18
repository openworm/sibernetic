/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013, 2026 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/

// ops_xpbd_full.cu — differentiable forward of the full XPBD pipeline.
//
// Per src/cuda/README.md Phase 4: this file exposes `xpbd_full_fwd` whose
// CLI signature matches metal_diff/sib_metal so the shared sgd_true.py
// runs against either substrate by only changing the BIN path. Force-based
// spring bonds via spring_K (NOT XPBD distance constraints), pair forces
// via visc_pair_coef, and the optional floor + restitution flags follow
// the Metal driver in ops_xpbd_full.mm.
//
// Tape layout per step k (matches Metal exactly):
//   x_old(3n) + v_old(3n) + density(n) + grad_C(3n) + denom_h(n)
//   [+ pair_density(n) if visc_pair_coef > 0]
//   [+ floor_clamped_mask(n) as int32 bit-cast into float if floor_y set]
// After K step blocks: (K+1)*n*3 trajectory frames, then final v(n*3).
//
// xpbd_full_bwd is sub-task F (not implemented here).
//
// Membrane args (n_membranes, n_elastic, r0, membranes.bin, pmem_index.bin)
// are PARSED but NO-OP'd — that's out of M10 scope per the spec scope rule.
// A stderr warning is emitted when n_membranes > 0 to alert the caller.

#include "cuda_common.h"
#include "shaders.h"
#include "ops.h"

// ──────────────────────────────────────────────────────────────────────
int run_xpbd_full_fwd(int argc, char **argv) {
    if (argc < 15 || argc > 26) {
        fprintf(stderr,
            "usage: sib_cuda xpbd_full_fwd "
            "<n_active> <n_static> <K> <h> <mass> <rho_rest> <dt> "
            "<gravity_y> <alpha_density> "
            "<pos0.bin> <vel0.bin> <pos_static.bin> <state_out.bin> "
            "[sim_scale] [visc_pair_coef] [spring_K] [bonds.bin] [floor_y] "
            "[restitution] [n_membranes] [n_elastic] [r0] "
            "[membranes.bin] [pmem_index.bin]\n"
            "  state_out.bin contains the per-step trajectory: \n"
            "    K x [x_old(n*3) + v_old(n*3) + density(n) + grad_C(n*3) + denom_h(n)"
            "    [+ pair_density(n) if visc_pair_coef>0]"
            "    [+ floor_clamped_mask(n) as int32 if floor_y is set]]\n"
            "    Plus K+1 frames of x_post for the trajectory.\n"
            "    Plus final v(n*3).\n"
            "  When spring_K > 0, bonds.bin must be provided (16 bytes per\n"
            "  bond: int32 i, int32 j, float32 rest_len, float32 _pad).\n"
            "  When floor_y is set, elastic floor at y=floor_y is applied\n"
            "  after density solve. restitution in [0, 1] (default 0).\n"
            "  Membrane args are accepted but no-op (out of scope; M10).\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)std::atoi(argv[2]);
    uint32_t n_static = (uint32_t)std::atoi(argv[3]);
    uint32_t K        = (uint32_t)std::atoi(argv[4]);
    float h           = (float)std::atof(argv[5]);
    float mass        = (float)std::atof(argv[6]);
    float rho_rest    = (float)std::atof(argv[7]);
    float dt          = (float)std::atof(argv[8]);
    float g_y         = (float)std::atof(argv[9]);
    float alpha_dens  = (float)std::atof(argv[10]);
    float sim_scale   = (argc >= 16) ? (float)std::atof(argv[15]) : 1.0f;
    float sim_scale_inv = 1.0f / sim_scale;
    float visc_pair_coef = (argc >= 17) ? (float)std::atof(argv[16]) : 0.0f;
    bool use_pair = visc_pair_coef != 0.0f;
    float spring_K = (argc >= 18) ? (float)std::atof(argv[17]) : 0.0f;
    const char *bonds_path = (argc >= 19) ? argv[18] : NULL;
    bool use_springs = spring_K != 0.0f;
    bool need_ext_accel = use_pair || use_springs;
    if (use_springs && !bonds_path) {
        fprintf(stderr, "spring_K > 0 requires bonds.bin path as 19th arg\n");
        return 1;
    }
    bool use_floor = (argc >= 20);
    float floor_y = use_floor ? (float)std::atof(argv[19]) : 0.0f;
    float restitution = (argc >= 21) ? (float)std::atof(argv[20]) : 0.0f;
    // Membrane args — parsed but no-op. Match Metal argc layout so the
    // shared SGD harness can pass either substrate the same argv.
    uint32_t n_membranes = (argc >= 22) ? (uint32_t)std::atoi(argv[21]) : 0u;
    if (n_membranes > 0) {
        fprintf(stderr,
                "warning: membranes ignored (out of scope)\n");
    }

    // Load bonds if springs enabled.
    uint32_t n_bonds = 0;
    int32_t *bond_ij_data = NULL;
    float *bond_rest_data = NULL;
    if (use_springs) {
        FILE *fb = std::fopen(bonds_path, "rb");
        if (!fb) { fprintf(stderr, "open %s\n", bonds_path); return 1; }
        std::fseek(fb, 0, SEEK_END);
        long sz = std::ftell(fb);
        std::fseek(fb, 0, SEEK_SET);
        n_bonds = (uint32_t)(sz / 16);
        uint8_t *raw = (uint8_t *)std::malloc((size_t)n_bonds * 16);
        std::fread(raw, 1, (size_t)n_bonds * 16, fb);
        std::fclose(fb);
        bond_ij_data = (int32_t *)std::malloc((size_t)n_bonds * 2 * sizeof(int32_t));
        bond_rest_data = (float *)std::malloc((size_t)n_bonds * sizeof(float));
        for (uint32_t b = 0; b < n_bonds; b++) {
            std::memcpy(&bond_ij_data[b * 2], raw + b * 16, 8);
            std::memcpy(&bond_rest_data[b], raw + b * 16 + 8, 4);
        }
        std::free(raw);
    }

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));
    float alpha_inv_dt2 = alpha_dens / (dt * dt);

    // Pair-force amps in fp64 to dodge fp32 underflow at small sim_scale.
    double h_scaled  = (double)h * (double)sim_scale;
    double h_s6      = std::pow(h_scaled, 6.0);
    double h_s9      = std::pow(h_scaled, 9.0);
    double divgradWvisco = 45.0 / (M_PI * h_s6);
    float visc_amp = (float)(1.5 * (double)mass * divgradWvisco
                              * std::pow((double)sim_scale, 3.0));
    double wpoly6_si = 315.0 / (64.0 * M_PI * h_s9);
    float surf_amp = (float)(-1.7e-9 * (double)mass * wpoly6_si
                              * (double)sim_scale / (double)mass);

    float *pos0 = read_floats_or_die(argv[11], (size_t)n_active * 3);
    float *vel0 = read_floats_or_die(argv[12], (size_t)n_active * 3);
    float *pos_static = (n_static > 0)
        ? read_floats_or_die(argv[13], (size_t)n_static * 3)
        : nullptr;

    size_t pos_b = (size_t)n_active * 3 * sizeof(float);
    size_t s_b   = (size_t)n_active * sizeof(float);

    // Per-step tape layout (matches Metal lines 174-178). Membranes are
    // intentionally absent from extra_per_step — they're no-op'd here.
    int extra_per_step = (use_pair ? 1 : 0) + (use_floor ? 1 : 0);
    size_t per_step_floats = (size_t)n_active * (3 + 3 + 1 + 3 + 1 + extra_per_step);
    float *state = (float *)std::calloc((size_t)K * per_step_floats, sizeof(float));

    // Trajectory: K+1 frames of x_post; frame 0 = pos0.
    float *traj = (float *)std::calloc((size_t)(K + 1) * (size_t)n_active * 3,
                                       sizeof(float));
    std::memcpy(traj, pos0, pos_b);

    // Build the static spatial grid one-time if pair forces are enabled.
    // Matches xpbd_step's grid setup and Metal forward.
    float *sorted_static_buf = nullptr;
    int *cell_start_buf = nullptr;
    GridDim3 grid_dim_struct = {0, 0, 0};
    GridOrigin3 grid_origin_struct = {0, 0, 0};
    int n_cells = 0;
    if (use_pair && n_static > 0) {
        build_static_grid(pos_static, n_static, h,
                          &sorted_static_buf, &cell_start_buf,
                          &grid_dim_struct, &grid_origin_struct, &n_cells);
    }

    // Device buffers.
    float3 *d_X = nullptr, *d_V = nullptr, *d_Xp = nullptr, *d_Sp = nullptr;
    float *d_R2aa = nullptr, *d_R2as = nullptr;
    float *d_Waa = nullptr, *d_Was = nullptr;
    float *d_D_aa = nullptr, *d_D = nullptr;
    float3 *d_Gc = nullptr;
    float *d_Dh = nullptr, *d_Lam = nullptr;
    float3 *d_ExtA = nullptr;
    int2 *d_BondIJ = nullptr;
    float *d_BondRest = nullptr;
    float3 *d_SortedS = nullptr;
    int *d_CellStart = nullptr;

    CUDA_CHECK(cudaMalloc(&d_X,  pos_b));
    CUDA_CHECK(cudaMalloc(&d_V,  pos_b));
    CUDA_CHECK(cudaMalloc(&d_Xp, pos_b));
    CUDA_CHECK(cudaMemcpy(d_X, pos0, pos_b, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_V, vel0, pos_b, cudaMemcpyHostToDevice));
    if (n_static > 0) {
        CUDA_CHECK(cudaMalloc(&d_Sp, (size_t)n_static * 3 * sizeof(float)));
        CUDA_CHECK(cudaMemcpy(d_Sp, pos_static,
                              (size_t)n_static * 3 * sizeof(float),
                              cudaMemcpyHostToDevice));
    }
    size_t r2aa_b = (size_t)n_active * n_active * sizeof(float);
    size_t r2as_b = (size_t)n_active * n_static * sizeof(float);
    CUDA_CHECK(cudaMalloc(&d_R2aa, r2aa_b));
    CUDA_CHECK(cudaMalloc(&d_Waa,  r2aa_b));
    CUDA_CHECK(cudaMalloc(&d_D_aa, s_b));
    CUDA_CHECK(cudaMalloc(&d_D,    s_b));
    CUDA_CHECK(cudaMalloc(&d_Gc,   pos_b));
    CUDA_CHECK(cudaMalloc(&d_Dh,   s_b));
    CUDA_CHECK(cudaMalloc(&d_Lam,  s_b));
    if (n_static > 0) {
        CUDA_CHECK(cudaMalloc(&d_R2as, r2as_b));
        CUDA_CHECK(cudaMalloc(&d_Was,  r2as_b));
    }
    if (need_ext_accel) {
        CUDA_CHECK(cudaMalloc(&d_ExtA, pos_b));
    }
    if (use_springs) {
        CUDA_CHECK(cudaMalloc(&d_BondIJ,   (size_t)n_bonds * sizeof(int2)));
        CUDA_CHECK(cudaMalloc(&d_BondRest, (size_t)n_bonds * sizeof(float)));
        // bond_ij_data is laid out as int32 pairs (i,j), matching int2's
        // layout exactly — same as ops_xpbd_step's bond load convention.
        CUDA_CHECK(cudaMemcpy(d_BondIJ, bond_ij_data,
                              (size_t)n_bonds * sizeof(int2),
                              cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_BondRest, bond_rest_data,
                              (size_t)n_bonds * sizeof(float),
                              cudaMemcpyHostToDevice));
    }
    if (use_pair && n_static > 0) {
        CUDA_CHECK(cudaMalloc(&d_SortedS,   (size_t)n_static * sizeof(float3)));
        CUDA_CHECK(cudaMalloc(&d_CellStart, (size_t)(n_cells + 1) * sizeof(int)));
        CUDA_CHECK(cudaMemcpy(d_SortedS, sorted_static_buf,
                              (size_t)n_static * 3 * sizeof(float),
                              cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_CellStart, cell_start_buf,
                              (size_t)(n_cells + 1) * sizeof(int),
                              cudaMemcpyHostToDevice));
    }

    // Seed bD with rho_rest so step 0's pair_forces has sensible 1/rho.
    // Matches Metal lines 341-344.
    if (use_pair) {
        std::vector<float> dens0((size_t)n_active, rho_rest);
        CUDA_CHECK(cudaMemcpy(d_D, dens0.data(), s_b,
                              cudaMemcpyHostToDevice));
    }

    const unsigned int TPB = 256;
    unsigned int gridA = (n_active + TPB - 1) / TPB;

    // Host scratch for tape readback. Avoid re-malloc per step.
    std::vector<float> h_pos_buf((size_t)n_active * 3);
    std::vector<float> h_vel_buf((size_t)n_active * 3);
    std::vector<float> h_dens_buf((size_t)n_active);
    std::vector<float> h_grad_buf((size_t)n_active * 3);
    std::vector<float> h_denomh_buf((size_t)n_active);
    std::vector<float> h_paird_buf((size_t)n_active);
    std::vector<float> h_xp_buf((size_t)n_active * 3);
    std::vector<int32_t> h_mask_buf((size_t)n_active);

    for (uint32_t k = 0; k < K; k++) {
        // SAVE x_old, v_old at start of step.
        float *step_state = state + (size_t)k * per_step_floats;
        CUDA_CHECK(cudaMemcpy(h_pos_buf.data(), d_X, pos_b,
                              cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_vel_buf.data(), d_V, pos_b,
                              cudaMemcpyDeviceToHost));
        std::memcpy(step_state,                  h_pos_buf.data(), pos_b);
        std::memcpy(step_state + 3 * n_active,   h_vel_buf.data(), pos_b);

        // SAVE pair_density (density seen by THIS step's pair_forces).
        if (use_pair) {
            CUDA_CHECK(cudaMemcpy(h_paird_buf.data(), d_D, s_b,
                                  cudaMemcpyDeviceToHost));
            std::memcpy(step_state + 11 * n_active, h_paird_buf.data(), s_b);
        }

        // Reset density Lagrange multipliers per step (XPBD convention).
        CUDA_CHECK(cudaMemset(d_Lam, 0, s_b));

        // (0) External forces accumulate into ext_accel, then advance v.
        if (need_ext_accel) {
            CUDA_CHECK(cudaMemset(d_ExtA, 0, pos_b));
        }
        if (use_pair) {
            // 32 threads per particle, one block per particle.
            // shaders.cu hardcodes a 5-round __shfl_down_sync (offsets
            // 16/8/4/2/1) that sums exactly one warp; anything other than
            // 32 silently miscomputes.
            constexpr unsigned int TPB_PAIR = 32;
            static_assert(TPB_PAIR == 32,
                          "pair_forces_grid_fwd: TPB must be 32 (single warp)");
            pair_forces_grid_fwd<<<n_active, TPB_PAIR>>>(
                d_X, d_V, d_SortedS, d_CellStart, d_D, d_ExtA,
                h, h2, sim_scale, visc_pair_coef, visc_amp, surf_amp,
                n_active,
                grid_dim_struct.x, grid_dim_struct.y, grid_dim_struct.z,
                grid_origin_struct.x, grid_origin_struct.y, grid_origin_struct.z);
        }
        if (use_springs) {
            const unsigned int TPB_S = 128;
            unsigned int grid_s = (n_active + TPB_S - 1) / TPB_S;
            spring_bonds_force<<<grid_s, TPB_S>>>(d_X, d_BondIJ, d_BondRest,
                                                  d_ExtA, spring_K, n_bonds,
                                                  n_active);
        }
        if (need_ext_accel) {
            apply_ext_accel<<<gridA, TPB>>>(d_V, d_ExtA, dt, n_active);
        }

        // (1) predict_positions: x_pred = x + dt*v + 0.5*dt^2*g (uses
        // gravity via vy_pred = v.y + g*dt, then ballistic step).
        predict_positions<<<gridA, TPB>>>(d_X, d_V, d_Xp,
                                          dt, g_y, n_active, sim_scale_inv);

        // (2) density chain. Metal's layout (lines 446-465):
        //   r2 -> Waa via blit, then wpoly6_inplace on Waa.
        //   R2aa keeps r^2 for density_constraint_grad below.
        {
            dim3 t2(16, 16);
            dim3 b2((n_active + 15) / 16, (n_active + 15) / 16);
            dist_active_active<<<b2, t2>>>(d_Xp, d_R2aa, n_active);
        }
        unsigned int n_aa_total = n_active * n_active;
        // Default-stream cudaMemcpyAsync would synchronize against the
        // next launch anyway; use plain cudaMemcpy so the error code
        // is checked synchronously.
        CUDA_CHECK(cudaMemcpy(d_Waa, d_R2aa, r2aa_b,
                              cudaMemcpyDeviceToDevice));
        wpoly6_inplace<<<(n_aa_total + TPB - 1) / TPB, TPB>>>(
            d_Waa, h2, poly6_const, n_aa_total);
        // rowsum_density / density_constraint_grad declare __shared__
        // float partials[256] and tree-reduce with stride = T/2. TPB > 256
        // writes OOB; TPB not a power of 2 <= 256 silently drops elements.
        constexpr unsigned int TPB_REDUCE = 256;
        static_assert(TPB_REDUCE == 256,
                      "rowsum_density/density_constraint_grad: TPB must be 256 "
                      "(shared partials[256] + power-of-2 tree reduce)");
        rowsum_density<<<n_active, TPB_REDUCE>>>(d_Waa, d_D_aa, mass,
                                                 n_active, n_active);
        if (n_static > 0) {
            {
                dim3 t2(16, 16);
                dim3 b2((n_active + 15) / 16, (n_static + 15) / 16);
                dist_active_static<<<b2, t2>>>(d_Xp, d_Sp, d_R2as,
                                               n_active, n_static);
            }
            unsigned int n_as_total = n_active * n_static;
            CUDA_CHECK(cudaMemcpy(d_Was, d_R2as, r2as_b,
                                  cudaMemcpyDeviceToDevice));
            wpoly6_inplace<<<(n_as_total + TPB - 1) / TPB, TPB>>>(
                d_Was, h2, poly6_const, n_as_total);
            rowsum_density<<<n_active, TPB_REDUCE>>>(d_Was, d_D, mass,
                                                     n_static, n_active);
            add_inplace<<<gridA, TPB>>>(d_D, d_D_aa, n_active);
        } else {
            // density = density_aa (memcpy).
            CUDA_CHECK(cudaMemcpy(d_D, d_D_aa, s_b,
                                  cudaMemcpyDeviceToDevice));
        }

        // (3) density_constraint_grad.
        density_constraint_grad<<<n_active, TPB_REDUCE>>>(
            d_Xp, d_Sp, d_R2aa, d_R2as,
            d_Gc, d_Dh, h, spiky_const, mass, rho_rest,
            n_active, n_static);

        // (4) solve_density_constraint.
        solve_density_constraint<<<gridA, TPB>>>(
            d_Xp, d_Lam, d_D, d_Gc, d_Dh,
            rho_rest, mass, alpha_inv_dt2, n_active);

        // (5) Floor mask + clamp (matches xpbd_step's solve_floor_constraint).
        // We compute the mask host-side BEFORE the clamp to avoid adding a
        // new mask-writing kernel. This trades one device→host copy of x_pred
        // for satisfying the "no new helper kernels" constraint.
        if (use_floor) {
            CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(h_xp_buf.data(), d_Xp, pos_b,
                                  cudaMemcpyDeviceToHost));
            for (uint32_t i = 0; i < n_active; i++) {
                h_mask_buf[i] = (h_xp_buf[3 * i + 1] < floor_y) ? 1 : 0;
            }
            solve_floor_constraint<<<gridA, TPB>>>(d_Xp, floor_y, n_active,
                                                   restitution);
        }

        // (6) update_velocities: v_new = (x_post - x_old) * sim_scale / dt.
        update_velocities<<<gridA, TPB>>>(d_V, d_X, d_Xp,
                                          dt, n_active, sim_scale);

        // Force sync before tape readback (mirrors Metal commit+wait).
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        // SAVE density, grad_C, denom_helper at offsets matching Metal.
        CUDA_CHECK(cudaMemcpy(h_dens_buf.data(),   d_D,  s_b,
                              cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_grad_buf.data(),   d_Gc, pos_b,
                              cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_denomh_buf.data(), d_Dh, s_b,
                              cudaMemcpyDeviceToHost));
        // Tape offset arithmetic (matches Metal lines 578-582):
        //   p = step_state + 6n
        //   p[0..n)      = density           (n floats)
        //   p[n..4n)     = grad_C            (3n floats)
        //   p[4n..5n)    = denom_h           (n floats)
        float *p = step_state + 6 * n_active;
        std::memcpy(p,              h_dens_buf.data(),   s_b);
        std::memcpy(p + n_active,   h_grad_buf.data(),   pos_b);
        std::memcpy(p + 4 * n_active, h_denomh_buf.data(), s_b);

        // Save floor mask (bit-cast int32 -> float, matches Metal line 586).
        if (use_floor) {
            size_t floor_off = (size_t)(11 + (use_pair ? 1 : 0)) * n_active;
            std::memcpy(step_state + floor_off, h_mask_buf.data(),
                        (size_t)n_active * sizeof(int32_t));
        }

        // Save x_post into trajectory frame k+1.
        CUDA_CHECK(cudaMemcpy(h_xp_buf.data(), d_Xp, pos_b,
                              cudaMemcpyDeviceToHost));
        std::memcpy(traj + (size_t)(k + 1) * n_active * 3,
                    h_xp_buf.data(), pos_b);

        // Advance: x_old <- x_post (device-side memcpy).
        CUDA_CHECK(cudaMemcpy(d_X, d_Xp, pos_b, cudaMemcpyDeviceToDevice));
    }

    // Read back final velocity.
    std::vector<float> final_v((size_t)n_active * 3);
    CUDA_CHECK(cudaMemcpy(final_v.data(), d_V, pos_b,
                          cudaMemcpyDeviceToHost));

    // Write outputs: state buffer, trajectory, final velocity (raw fwrite).
    FILE *fs = std::fopen(argv[14], "wb");
    if (!fs) { fprintf(stderr, "cannot open %s for write\n", argv[14]);
               return 1; }
    std::fwrite(state, sizeof(float), (size_t)K * per_step_floats, fs);
    std::fwrite(traj,  sizeof(float), (size_t)(K + 1) * n_active * 3, fs);
    std::fwrite(final_v.data(), sizeof(float), (size_t)n_active * 3, fs);
    std::fclose(fs);

    // Cleanup.
    cudaFree(d_X); cudaFree(d_V); cudaFree(d_Xp);
    if (d_Sp) cudaFree(d_Sp);
    cudaFree(d_R2aa); cudaFree(d_Waa);
    if (d_R2as) cudaFree(d_R2as);
    if (d_Was)  cudaFree(d_Was);
    cudaFree(d_D_aa); cudaFree(d_D);
    cudaFree(d_Gc);   cudaFree(d_Dh); cudaFree(d_Lam);
    if (d_ExtA) cudaFree(d_ExtA);
    if (d_BondIJ)   cudaFree(d_BondIJ);
    if (d_BondRest) cudaFree(d_BondRest);
    if (d_SortedS)   cudaFree(d_SortedS);
    if (d_CellStart) cudaFree(d_CellStart);

    std::free(pos0); std::free(vel0);
    if (pos_static) std::free(pos_static);
    std::free(state); std::free(traj);
    if (sorted_static_buf) std::free(sorted_static_buf);
    if (cell_start_buf)    std::free(cell_start_buf);
    if (bond_ij_data)   std::free(bond_ij_data);
    if (bond_rest_data) std::free(bond_rest_data);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// Differentiable XPBD reverse-mode pass. Walks K steps in reverse,
// running the matching backward kernel at each step, accumulating
// ∂L/∂{x_init, v_init, rho_rest, spring_K, visc_pair_coef, alpha_dens}.
//
// CLI signature matches Metal's xpbd_full_bwd exactly (see
// src/metal_diff/ops_xpbd_full.mm:631-645). Tape layout, optional arg
// positions, and reverse kernel sequence all mirror Metal so that the
// shared sgd_true.py operates on either binary by swapping BIN path.
//
// Membrane args are PARSED but NO-OP'd (out of scope per spec scope rule).
// A stderr warning is emitted when n_membranes > 0.
//
// BWD_CLIP_NORM env var: if set to positive float, per-step L2 clip of
// (||grad_x||² + ||grad_v||²)^0.5 after each k iteration. Matches
// --clip_norm in sgd_true.py.
// ──────────────────────────────────────────────────────────────────────
int run_xpbd_full_bwd(int argc, char **argv) {
    if (argc < 17 || argc > 31) {
        fprintf(stderr,
            "usage: sib_cuda xpbd_full_bwd "
            "<n_active> <n_static> <K> <h> <mass> <rho_rest> <dt> "
            "<gravity_y> <alpha_density> "
            "<state_in.bin> <pos_static.bin> <grad_x_final.bin> "
            "<grad_x_init_out.bin> <grad_v_init_out.bin> <grad_rho_out.bin> "
            "[sim_scale] [visc_pair_coef] [spring_K] [bonds.bin] "
            "[grad_spring_K_out.bin] [grad_visc_K_out.bin] [floor_y] "
            "[grad_alpha_dens_out.bin] [restitution] "
            "[n_membranes] [n_elastic] [r0] [membranes.bin] [pmem_index.bin]\n"
            "  (must match the xpbd_full_fwd args used to produce the "
            "state file)\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)std::atoi(argv[2]);
    uint32_t n_static = (uint32_t)std::atoi(argv[3]);
    uint32_t K        = (uint32_t)std::atoi(argv[4]);
    float h           = (float)std::atof(argv[5]);
    float mass        = (float)std::atof(argv[6]);
    float rho_rest    = (float)std::atof(argv[7]);
    float dt          = (float)std::atof(argv[8]);
    float g_y         = (float)std::atof(argv[9]);
    float alpha_dens  = (float)std::atof(argv[10]);
    // Optional positional args — match Metal layout exactly.
    float sim_scale   = (argc >= 18) ? (float)std::atof(argv[17]) : 1.0f;
    float sim_scale_inv = 1.0f / sim_scale;
    float visc_pair_coef = (argc >= 19) ? (float)std::atof(argv[18]) : 0.0f;
    bool use_pair = visc_pair_coef != 0.0f;
    float spring_K = (argc >= 20) ? (float)std::atof(argv[19]) : 0.0f;
    const char *bonds_path = (argc >= 21) ? argv[20] : NULL;
    bool use_springs = spring_K != 0.0f;
    bool need_ext_accel = use_pair || use_springs;
    if (use_springs && !bonds_path) {
        fprintf(stderr, "spring_K > 0 requires bonds.bin path as 21st arg\n");
        return 1;
    }
    const char *grad_spring_K_out_path = (argc >= 22) ? argv[21] : NULL;
    const char *grad_visc_K_out_path   = (argc >= 23) ? argv[22] : NULL;
    bool use_floor = (argc >= 24);
    float floor_y = use_floor ? (float)std::atof(argv[23]) : 0.0f;
    const char *grad_alpha_dens_out_path = (argc >= 25) ? argv[24] : NULL;
    float restitution = (argc >= 26) ? (float)std::atof(argv[25]) : 0.0f;
    uint32_t n_membranes = (argc >= 27) ? (uint32_t)std::atoi(argv[26]) : 0u;
    if (n_membranes > 0) {
        fprintf(stderr, "warning: membranes ignored (out of scope)\n");
    }

    // Load bonds (matches fwd loader).
    uint32_t n_bonds = 0;
    int32_t *bond_ij_data = NULL;
    float *bond_rest_data = NULL;
    if (use_springs) {
        FILE *fb = std::fopen(bonds_path, "rb");
        if (!fb) { fprintf(stderr, "open %s\n", bonds_path); return 1; }
        std::fseek(fb, 0, SEEK_END);
        long sz = std::ftell(fb);
        std::fseek(fb, 0, SEEK_SET);
        n_bonds = (uint32_t)(sz / 16);
        uint8_t *raw = (uint8_t *)std::malloc((size_t)n_bonds * 16);
        std::fread(raw, 1, (size_t)n_bonds * 16, fb);
        std::fclose(fb);
        bond_ij_data = (int32_t *)std::malloc((size_t)n_bonds * 2 * sizeof(int32_t));
        bond_rest_data = (float *)std::malloc((size_t)n_bonds * sizeof(float));
        for (uint32_t b = 0; b < n_bonds; b++) {
            std::memcpy(&bond_ij_data[b * 2], raw + b * 16, 8);
            std::memcpy(&bond_rest_data[b], raw + b * 16 + 8, 4);
        }
        std::free(raw);
    }

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));
    float alpha_inv_dt2 = alpha_dens / (dt * dt);
    // Pair-force amps (mirrors fwd).
    double h_scaled  = (double)h * (double)sim_scale;
    double h_s6      = std::pow(h_scaled, 6.0);
    double h_s9      = std::pow(h_scaled, 9.0);
    double divgradWvisco = 45.0 / (M_PI * h_s6);
    float visc_amp = (float)(1.5 * (double)mass * divgradWvisco
                              * std::pow((double)sim_scale, 3.0));
    double wpoly6_si = 315.0 / (64.0 * M_PI * h_s9);
    float surf_amp = (float)(-1.7e-9 * (double)mass * wpoly6_si
                              * (double)sim_scale / (double)mass);

    // Tape layout (mirrors fwd writes — see ops_xpbd_full.cu fwd section).
    // Membranes intentionally absent from extra_per_step.
    int extra_per_step = (use_pair ? 1 : 0) + (use_floor ? 1 : 0);
    size_t per_step_floats = (size_t)n_active * (3 + 3 + 1 + 3 + 1 + extra_per_step);
    size_t state_size = (size_t)K * per_step_floats
                      + (size_t)(K + 1) * n_active * 3
                      + (size_t)n_active * 3;
    float *all = read_floats_or_die(argv[11], state_size);
    float *state = all;
    // (traj at all + K*per_step_floats; vel_final after — not used here.)

    float *pos_static = (n_static > 0)
        ? read_floats_or_die(argv[12], (size_t)n_static * 3) : nullptr;
    float *grad_x_fin = read_floats_or_die(argv[13], (size_t)n_active * 3);

    size_t pos_b = (size_t)n_active * 3 * sizeof(float);
    size_t s_b   = (size_t)n_active * sizeof(float);

    // Build static spatial grid (matches fwd / pair_forces_grid_bwd convention).
    float *sorted_static_buf = nullptr;
    int *cell_start_buf = nullptr;
    GridDim3 grid_dim_struct = {0, 0, 0};
    GridOrigin3 grid_origin_struct = {0, 0, 0};
    int n_cells = 0;
    if (use_pair && n_static > 0) {
        build_static_grid(pos_static, n_static, h,
                          &sorted_static_buf, &cell_start_buf,
                          &grid_dim_struct, &grid_origin_struct, &n_cells);
    }

    // ── Device buffers ───────────────────────────────────────────────
    // Persistent across steps.
    float3 *d_X = nullptr, *d_V = nullptr, *d_Xp = nullptr;
    float3 *d_Sp = nullptr;
    float *d_R2aa = nullptr, *d_R2as = nullptr;
    float *d_GW_aa = nullptr, *d_GW_as = nullptr;
    float *d_D = nullptr;        // density / pair_density slot
    float3 *d_Gc = nullptr;
    float *d_Dh = nullptr, *d_Lam = nullptr;
    float3 *d_ExtA = nullptr;
    int2 *d_BondIJ = nullptr;
    float *d_BondRest = nullptr;
    float3 *d_SortedS = nullptr;
    int *d_CellStart = nullptr;
    // Gradient buffers (per-step working + running).
    float3 *d_Gv_in = nullptr;
    float3 *d_Gx_pred = nullptr, *d_Gx_old_new = nullptr, *d_Gv_old_new = nullptr;
    float3 *d_Gx_running = nullptr, *d_Gv_running = nullptr;
    float *d_Glam_pre = nullptr, *d_Gdens = nullptr;
    float3 *d_GgC = nullptr;
    float *d_Gdh = nullptr;
    float *d_Grho_per = nullptr, *d_Galpha_per = nullptr;
    float *d_Ggy_per = nullptr;          // grad_gravity_y atomicAdd target (unused output)
    float3 *d_Gext = nullptr;            // grad_ext_accel scratch
    float *d_SpringPart = nullptr, *d_ViscPart = nullptr;
    float3 *d_PosPredSaved = nullptr;    // synthetic positions for floor backward
    float *d_FloorYGrad = nullptr;       // atomicAdd target (unused output)
    float *d_RestitGrad = nullptr;       // atomicAdd target (unused output)

    CUDA_CHECK(cudaMalloc(&d_X,  pos_b));
    CUDA_CHECK(cudaMalloc(&d_V,  pos_b));
    CUDA_CHECK(cudaMalloc(&d_Xp, pos_b));
    if (n_static > 0) {
        CUDA_CHECK(cudaMalloc(&d_Sp, (size_t)n_static * 3 * sizeof(float)));
        CUDA_CHECK(cudaMemcpy(d_Sp, pos_static,
                              (size_t)n_static * 3 * sizeof(float),
                              cudaMemcpyHostToDevice));
    }
    size_t r2aa_b = (size_t)n_active * n_active * sizeof(float);
    size_t r2as_b = (size_t)n_active * n_static * sizeof(float);
    CUDA_CHECK(cudaMalloc(&d_R2aa,  r2aa_b));
    CUDA_CHECK(cudaMalloc(&d_GW_aa, r2aa_b));
    if (n_static > 0) {
        CUDA_CHECK(cudaMalloc(&d_R2as,  r2as_b));
        CUDA_CHECK(cudaMalloc(&d_GW_as, r2as_b));
    }
    CUDA_CHECK(cudaMalloc(&d_D,   s_b));
    CUDA_CHECK(cudaMalloc(&d_Gc,  pos_b));
    CUDA_CHECK(cudaMalloc(&d_Dh,  s_b));
    CUDA_CHECK(cudaMalloc(&d_Lam, s_b));
    CUDA_CHECK(cudaMemset(d_Lam, 0, s_b));  // λ_pre always 0 (single-iter XPBD)

    CUDA_CHECK(cudaMalloc(&d_Gv_in,      pos_b));
    CUDA_CHECK(cudaMalloc(&d_Gx_pred,    pos_b));
    CUDA_CHECK(cudaMalloc(&d_Gx_old_new, pos_b));
    CUDA_CHECK(cudaMalloc(&d_Gv_old_new, pos_b));
    CUDA_CHECK(cudaMalloc(&d_Gx_running, pos_b));
    CUDA_CHECK(cudaMalloc(&d_Gv_running, pos_b));
    // Seed gradient: grad_x_running = grad_x_final, grad_v_running = 0.
    CUDA_CHECK(cudaMemcpy(d_Gx_running, grad_x_fin, pos_b,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_Gv_running, 0, pos_b));

    CUDA_CHECK(cudaMalloc(&d_Glam_pre,   s_b));
    CUDA_CHECK(cudaMalloc(&d_Gdens,      s_b));
    CUDA_CHECK(cudaMalloc(&d_GgC,        pos_b));
    CUDA_CHECK(cudaMalloc(&d_Gdh,        s_b));
    CUDA_CHECK(cudaMalloc(&d_Grho_per,   s_b));
    CUDA_CHECK(cudaMalloc(&d_Galpha_per, s_b));
    CUDA_CHECK(cudaMalloc(&d_Ggy_per,    sizeof(float)));

    if (need_ext_accel) {
        CUDA_CHECK(cudaMalloc(&d_ExtA, pos_b));
        CUDA_CHECK(cudaMalloc(&d_Gext, pos_b));
    }
    if (use_springs) {
        CUDA_CHECK(cudaMalloc(&d_BondIJ,   (size_t)n_bonds * sizeof(int2)));
        CUDA_CHECK(cudaMalloc(&d_BondRest, (size_t)n_bonds * sizeof(float)));
        CUDA_CHECK(cudaMemcpy(d_BondIJ, bond_ij_data,
                              (size_t)n_bonds * sizeof(int2),
                              cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_BondRest, bond_rest_data,
                              (size_t)n_bonds * sizeof(float),
                              cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMalloc(&d_SpringPart, s_b));
    }
    if (use_pair) {
        CUDA_CHECK(cudaMalloc(&d_ViscPart, s_b));
    }
    if (use_pair && n_static > 0) {
        CUDA_CHECK(cudaMalloc(&d_SortedS,   (size_t)n_static * sizeof(float3)));
        CUDA_CHECK(cudaMalloc(&d_CellStart, (size_t)(n_cells + 1) * sizeof(int)));
        CUDA_CHECK(cudaMemcpy(d_SortedS, sorted_static_buf,
                              (size_t)n_static * 3 * sizeof(float),
                              cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_CellStart, cell_start_buf,
                              (size_t)(n_cells + 1) * sizeof(int),
                              cudaMemcpyHostToDevice));
    }
    if (use_floor) {
        CUDA_CHECK(cudaMalloc(&d_PosPredSaved, pos_b));
        CUDA_CHECK(cudaMalloc(&d_FloorYGrad, sizeof(float)));
        CUDA_CHECK(cudaMalloc(&d_RestitGrad, sizeof(float)));
    }

    // Host scratch buffers for per-step tape reads.
    std::vector<float> h_x_old((size_t)n_active * 3);
    std::vector<float> h_v_old((size_t)n_active * 3);
    std::vector<float> h_dens((size_t)n_active);
    std::vector<float> h_gradC((size_t)n_active * 3);
    std::vector<float> h_denomh((size_t)n_active);
    std::vector<float> h_paird((size_t)n_active);
    std::vector<int32_t> h_mask((size_t)n_active);
    std::vector<float> h_pos_synth((size_t)n_active * 3);

    // Host-side accumulators (matches Metal's host-side scalars).
    float total_grad_rho = 0.0f;
    float total_grad_spring_K = 0.0f;
    float total_grad_visc_K = 0.0f;
    float total_grad_alpha_inv_dt2 = 0.0f;

    // BWD_CLIP_NORM: per-step L2 clip on grad_x and grad_v INDEPENDENTLY.
    // Defaults to 1e3 (always on) — matches Metal at
    // src/metal_diff/ops_xpbd_full.mm:996-997. Required for chaotic SPH
    // dynamics: per-step amplification ~3x means K=50000 unrolled gradient
    // overflows fp32 within ~30 steps without clipping. NaN/Inf scrubbing
    // runs whenever the clip is on. Set BWD_CLIP_NORM=0 to disable.
    float clip_norm_x = 1e3f;
    float clip_norm_v = 1e3f;
    bool clip_on = true;
    const char *clip_env = std::getenv("BWD_CLIP_NORM");
    if (clip_env) {
        float v = (float)std::atof(clip_env);
        if (v <= 0.0f) {
            clip_on = false;
            fprintf(stderr, "[xpbd_full_bwd] per-step gradient clipping "
                            "DISABLED via BWD_CLIP_NORM=0\n");
        } else {
            clip_norm_x = v;
            clip_norm_v = v;
            fprintf(stderr, "[xpbd_full_bwd] per-step gradient clipping: "
                            "|grad_x|, |grad_v| L2 each capped at %.3e\n", v);
        }
    } else {
        fprintf(stderr, "[xpbd_full_bwd] per-step gradient clipping: "
                        "|grad_x|, |grad_v| L2 each capped at 1.000e+03 "
                        "(default)\n");
    }

    // BWD_TBPTT: truncated BPTT — chain back only the last `max_bw_steps`
    // steps. Default = K (no truncation). Mirror of Metal at
    // src/metal_diff/ops_xpbd_full.mm:974-986. Critical for K large
    // enough that the gradient chain saturates clipping bound; sgd_true.py
    // sets BWD_TBPTT=50 by default.
    int32_t max_bw_steps = (int32_t)K;
    const char *tbptt_env = std::getenv("BWD_TBPTT");
    if (tbptt_env) {
        max_bw_steps = std::atoi(tbptt_env);
        if (max_bw_steps <= 0 || max_bw_steps > (int32_t)K) {
            max_bw_steps = (int32_t)K;
        }
        fprintf(stderr, "[xpbd_full_bwd] truncated BPTT: %d / %u steps\n",
                max_bw_steps, K);
    }
    int32_t k_stop = (int32_t)K - max_bw_steps;

    const unsigned int TPB = 256;
    unsigned int gridA = (n_active + TPB - 1) / TPB;

    // Walk K steps in reverse.
    for (int32_t k = (int32_t)K - 1; k >= k_stop; k--) {
        float *step_state = state + (size_t)k * per_step_floats;
        // Load tape entries.
        std::memcpy(h_x_old.data(),  step_state,                pos_b);
        std::memcpy(h_v_old.data(),  step_state + 3 * n_active, pos_b);
        std::memcpy(h_dens.data(),   step_state + 6 * n_active, s_b);
        std::memcpy(h_gradC.data(),  step_state + 7 * n_active, pos_b);
        std::memcpy(h_denomh.data(), step_state + 10 * n_active, s_b);
        if (use_pair) {
            std::memcpy(h_paird.data(), step_state + 11 * n_active, s_b);
        }
        if (use_floor) {
            size_t floor_off = (size_t)(11 + (use_pair ? 1 : 0)) * n_active;
            std::memcpy(h_mask.data(), step_state + floor_off,
                        (size_t)n_active * sizeof(int32_t));
        }

        CUDA_CHECK(cudaMemcpy(d_X, h_x_old.data(), pos_b, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_V, h_v_old.data(), pos_b, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_D,  h_dens.data(),   s_b,   cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Gc, h_gradC.data(),  pos_b, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Dh, h_denomh.data(), s_b,   cudaMemcpyHostToDevice));

        // Recompute x_pre via the SAME fwd chain that produced it:
        // pair_forces → spring_bonds → apply_ext_accel → predict_positions.
        // Needed because density chain backwards want x_pre (pre-clamp).
        if (need_ext_accel) {
            CUDA_CHECK(cudaMemset(d_ExtA, 0, pos_b));
        }
        if (use_pair) {
            // Use this step's saved pair_density (= density BEFORE this
            // step's density solve, mirroring fwd ordering).
            CUDA_CHECK(cudaMemcpy(d_D, h_paird.data(), s_b,
                                  cudaMemcpyHostToDevice));
            // 32 threads per particle; single-warp __shfl_down_sync invariant.
            constexpr unsigned int TPB_PAIR = 32;
            static_assert(TPB_PAIR == 32,
                          "pair_forces_grid_fwd: TPB must be 32 (single warp)");
            pair_forces_grid_fwd<<<n_active, TPB_PAIR>>>(
                d_X, d_V, d_SortedS, d_CellStart, d_D, d_ExtA,
                h, h2, sim_scale, visc_pair_coef, visc_amp, surf_amp,
                n_active,
                grid_dim_struct.x, grid_dim_struct.y, grid_dim_struct.z,
                grid_origin_struct.x, grid_origin_struct.y, grid_origin_struct.z);
        }
        if (use_springs) {
            const unsigned int TPB_S = 128;
            unsigned int grid_s = (n_active + TPB_S - 1) / TPB_S;
            spring_bonds_force<<<grid_s, TPB_S>>>(d_X, d_BondIJ, d_BondRest,
                                                  d_ExtA, spring_K, n_bonds,
                                                  n_active);
        }
        if (need_ext_accel) {
            apply_ext_accel<<<gridA, TPB>>>(d_V, d_ExtA, dt, n_active);
        }
        predict_positions<<<gridA, TPB>>>(d_X, d_V, d_Xp,
                                          dt, g_y, n_active, sim_scale_inv);
        // Recompute r²_aa, r²_as for density chain backward.
        {
            dim3 t2(16, 16);
            dim3 b2((n_active + 15) / 16, (n_active + 15) / 16);
            dist_active_active<<<b2, t2>>>(d_Xp, d_R2aa, n_active);
        }
        if (n_static > 0) {
            dim3 t2(16, 16);
            dim3 b2((n_active + 15) / 16, (n_static + 15) / 16);
            dist_active_static<<<b2, t2>>>(d_Xp, d_Sp, d_R2as,
                                           n_active, n_static);
        }

        // Restore the SOLVED density (used by solve_density_constraint_bwd).
        // The pair_density we used above is the PRE-solve density for ext_accel.
        // solve_density_bwd needs the POST-solve density saved on tape.
        CUDA_CHECK(cudaMemcpy(d_D, h_dens.data(), s_b, cudaMemcpyHostToDevice));

        // ── Zero per-step working gradient buffers ──────────────────
        CUDA_CHECK(cudaMemset(d_Gx_old_new, 0, pos_b));
        CUDA_CHECK(cudaMemset(d_Gv_old_new, 0, pos_b));
        CUDA_CHECK(cudaMemset(d_Gx_pred,    0, pos_b));
        CUDA_CHECK(cudaMemset(d_Ggy_per,    0, sizeof(float)));
        if (need_ext_accel) {
            CUDA_CHECK(cudaMemset(d_Gext, 0, pos_b));
        }
        // Move running v gradient into bGv_in for update_velocities_backward.
        CUDA_CHECK(cudaMemcpy(d_Gv_in, d_Gv_running, pos_b,
                              cudaMemcpyDeviceToDevice));

        // ── (a) Floor backward (BEFORE update_velocities_backward) ───
        // The CUDA solve_floor_constraint_backward kernel takes positions
        // (not a precomputed mask). It computes the mask internally via
        // `xp.y < floor_y`. We synthesize positions from the saved mask:
        // clamped → y = floor_y - 1 (true), non-clamped → y = floor_y + 1.
        // The kernel WRITES (not accumulates) to its grad_pos_pred output,
        // and atomicAdd's to floor_y/restitution grad scalars (discarded).
        if (use_floor) {
            for (uint32_t i = 0; i < n_active; i++) {
                h_pos_synth[3 * i + 0] = 0.0f;
                h_pos_synth[3 * i + 1] = h_mask[i] ? (floor_y - 1.0f)
                                                   : (floor_y + 1.0f);
                h_pos_synth[3 * i + 2] = 0.0f;
            }
            CUDA_CHECK(cudaMemcpy(d_PosPredSaved, h_pos_synth.data(), pos_b,
                                  cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemset(d_FloorYGrad, 0, sizeof(float)));
            CUDA_CHECK(cudaMemset(d_RestitGrad, 0, sizeof(float)));
            // Write into Gx_pred as scratch — copy back to running afterwards.
            solve_floor_constraint_backward<<<gridA, TPB>>>(
                d_Gx_running,      // grad_pos_post (input)
                d_PosPredSaved,    // synthetic positions
                d_Gx_pred,         // grad_pos_pred (output, written)
                d_FloorYGrad, d_RestitGrad,
                floor_y, restitution, n_active);
            CUDA_CHECK(cudaGetLastError());
            // Replace Gx_running ← Gx_pred (floor-aware gradient), then reset Gx_pred.
            CUDA_CHECK(cudaMemcpy(d_Gx_running, d_Gx_pred, pos_b,
                                  cudaMemcpyDeviceToDevice));
            CUDA_CHECK(cudaMemset(d_Gx_pred, 0, pos_b));
        }

        // ── (b) update_velocities_backward ──────────────────────────
        // Accumulates: Gx_running += gvn * sim_scale/dt  (∂L/∂x_post)
        //              Gx_old_new -= gvn * sim_scale/dt  (∂L/∂x_old)
        //
        // KNOWN DEFECT (cross-backend parity workaround, tracked):
        //   The Metal backward at src/metal_diff/shaders.metal:2216
        //   computes g_v / dt instead of g_v * sim_scale / dt — the
        //   forward at :1942 multiplies by sim_scale but the backward
        //   drops it. Our CUDA kernel is mathematically correct
        //   (computes g_v * sim_scale / dt); we pass sim_scale = 1.0f
        //   here so that CUDA's output matches Metal's to within the
        //   README Phase-5 "1% parity" claim. Adam in log-space
        //   absorbs the uniform 1/sim_scale factor, so demo1 optima
        //   are unaffected.
        //
        //   To remove: patch Metal upstream (forward and backward both
        //   need sim_scale), revert this argument from 1.0f back to
        //   `sim_scale`, and tighten the Phase-5 parity test. The
        //   env-var escape hatch CUDA_BWD_TRUE_SIM_SCALE=1 below lets
        //   downstream consumers opt into the mathematically-correct
        //   gradient today, at the cost of failing Metal parity.
        float upd_sim_scale = 1.0f;
        const char *true_ss_env = std::getenv("CUDA_BWD_TRUE_SIM_SCALE");
        if (true_ss_env && std::atoi(true_ss_env) != 0) {
            upd_sim_scale = sim_scale;
        }
        update_velocities_backward<<<gridA, TPB>>>(
            d_Gv_in, d_Gx_running, d_Gx_old_new,
            dt, upd_sim_scale, n_active);

        // ── (c) solve_density_constraint_backward ───────────────────
        // ∂L/∂λ_post = 0 (single-iter). Reuse d_Lam as zeros.
        solve_density_constraint_bwd<<<gridA, TPB>>>(
            d_D, d_Gc, d_Dh, d_Lam,
            d_Gx_running,    // grad_pos_post (input)
            d_Lam,           // grad_lambda_post = 0 (reuse Lam as zeros)
            d_Gx_pred,       // grad_pos_pre (accum)
            d_Glam_pre, d_Gdens, d_GgC, d_Gdh, d_Grho_per, d_Galpha_per,
            rho_rest, mass, alpha_inv_dt2, n_active);

        // ── (d) density_constraint_grad_bwd ─────────────────────────
        density_constraint_grad_bwd<<<gridA, TPB>>>(
            d_Xp, d_Sp, d_R2aa, d_R2as,
            d_GgC, d_Gdh,
            d_Gx_pred,       // accumulates
            h, spiky_const, mass, rho_rest,
            n_active, n_static);

        // ── (e) density chain backward (rowsum + wpoly6 + dist) ─────
        // aa branch:
        {
            dim3 t2(16, 16);
            dim3 b2((n_active + 15) / 16, (n_active + 15) / 16);
            rowsum_density_bwd<<<b2, t2>>>(d_Gdens, d_GW_aa, mass,
                                           n_active, n_active);
        }
        unsigned int n_aa_total = n_active * n_active;
        wpoly6_inplace_bwd<<<(n_aa_total + TPB - 1) / TPB, TPB>>>(
            d_R2aa, d_GW_aa, h2, poly6_const, n_aa_total);
        dist_active_active_bwd<<<gridA, TPB>>>(
            d_Xp, d_GW_aa, d_Gx_pred, n_active);
        if (n_static > 0) {
            {
                dim3 t2(16, 16);
                dim3 b2((n_static + 15) / 16, (n_active + 15) / 16);
                rowsum_density_bwd<<<b2, t2>>>(d_Gdens, d_GW_as, mass,
                                               n_active, n_static);
            }
            unsigned int n_as_total = n_active * n_static;
            wpoly6_inplace_bwd<<<(n_as_total + TPB - 1) / TPB, TPB>>>(
                d_R2as, d_GW_as, h2, poly6_const, n_as_total);
            // dist_active_static_bwd: one block per active particle,
            // 256 threads cooperate over the n_static axis. TPB_REDUCE
            // must be 256 to match shared partials[256] + power-of-2
            // tree reduce.
            constexpr unsigned int TPB_REDUCE_DAS = 256;
            static_assert(TPB_REDUCE_DAS == 256,
                          "dist_active_static_bwd: TPB must be 256 "
                          "(shared partials[256] + power-of-2 tree reduce)");
            dist_active_static_bwd<<<n_active, TPB_REDUCE_DAS>>>(
                d_Xp, d_Sp, d_GW_as, d_Gx_pred, n_active, n_static);
        }

        // ── (f) predict_positions_backward ──────────────────────────
        //   Gx_old_new += Gx_pred  (accumulates)
        //   Gv_old_new = Gx_pred * dt * sim_scale_inv  (WRITES)
        // The kernel writes grad_vel — that's what we want for v_after_apply.
        predict_positions_backward<<<gridA, TPB>>>(
            d_Gx_pred, d_Gx_old_new, d_Gv_old_new, d_Ggy_per,
            dt, sim_scale_inv, n_active);

        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        // ── (g) ext_accel backward ──────────────────────────────────
        // Mirrors Metal lines 1246-1357. After predict_bw, Gv_old_new holds
        // ∂L/∂v_after_apply. We DERIVE Gext = dt · Gv_old_new (NOT an
        // apply_ext_accel_backward call — Metal does this on host as a
        // simple multiply). Then pair/spring backwards take Gext as input
        // and ACCUMULATE into Gx_old_new (position deps) and Gv_old_new
        // (velocity deps, pair only).
        if (need_ext_accel) {
            // Compute Gext = dt * Gv_old_new host-side (one memcpy round-trip).
            // Could be a kernel but Metal does it host-side too.
            std::vector<float> h_gv_new((size_t)n_active * 3);
            std::vector<float> h_gext((size_t)n_active * 3);
            CUDA_CHECK(cudaMemcpy(h_gv_new.data(), d_Gv_old_new, pos_b,
                                  cudaMemcpyDeviceToHost));
            for (uint32_t i = 0; i < 3u * n_active; i++) {
                h_gext[i] = dt * h_gv_new[i];
            }
            CUDA_CHECK(cudaMemcpy(d_Gext, h_gext.data(), pos_b,
                                  cudaMemcpyHostToDevice));

            if (use_pair) {
                // Restore this step's pair_density (the density seen by
                // pair_forces in fwd) — solve_density_bwd has overwritten
                // d_D usage but we re-load to be explicit.
                CUDA_CHECK(cudaMemcpy(d_D, h_paird.data(), s_b,
                                      cudaMemcpyHostToDevice));
                const unsigned int TPB_P = 64;
                unsigned int grid_p = (n_active + TPB_P - 1) / TPB_P;
                pair_forces_grid_bwd<<<grid_p, TPB_P>>>(
                    d_X, d_V, d_SortedS, d_CellStart, d_D,
                    d_Gext, d_Gx_old_new, d_Gv_old_new,
                    h, h2, sim_scale, visc_pair_coef, visc_amp, surf_amp,
                    n_active,
                    grid_dim_struct.x, grid_dim_struct.y, grid_dim_struct.z,
                    grid_origin_struct.x, grid_origin_struct.y, grid_origin_struct.z);
            }
            if (use_springs) {
                const unsigned int TPB_S = 128;
                unsigned int grid_s = (n_active + TPB_S - 1) / TPB_S;
                spring_bonds_force_backward<<<grid_s, TPB_S>>>(
                    d_X, d_BondIJ, d_BondRest,
                    d_Gext, d_Gx_old_new,
                    spring_K, n_bonds, n_active);
            }
            CUDA_CHECK(cudaGetLastError());
            CUDA_CHECK(cudaDeviceSynchronize());

            // Param-grad partials (∂L/∂spring_K, ∂L/∂visc_pair_coef).
            if (use_springs) {
                const unsigned int TPB_S = 128;
                unsigned int grid_s = (n_active + TPB_S - 1) / TPB_S;
                spring_K_partial<<<grid_s, TPB_S>>>(
                    d_X, d_BondIJ, d_BondRest, d_Gext, d_SpringPart,
                    n_bonds, n_active);
            }
            if (use_pair) {
                const unsigned int TPB_P = 64;
                unsigned int grid_p = (n_active + TPB_P - 1) / TPB_P;
                visc_K_partial<<<grid_p, TPB_P>>>(
                    d_X, d_V, d_SortedS, d_CellStart, d_D, d_Gext, d_ViscPart,
                    h, h2, sim_scale, visc_amp,
                    n_active,
                    grid_dim_struct.x, grid_dim_struct.y, grid_dim_struct.z,
                    grid_origin_struct.x, grid_origin_struct.y, grid_origin_struct.z);
            }
            CUDA_CHECK(cudaGetLastError());
            CUDA_CHECK(cudaDeviceSynchronize());

            // Host-side sum of per-particle param-grad partials.
            if (use_springs) {
                std::vector<float> h_sp((size_t)n_active);
                CUDA_CHECK(cudaMemcpy(h_sp.data(), d_SpringPart, s_b,
                                      cudaMemcpyDeviceToHost));
                float s = 0.0f;
                for (uint32_t i = 0; i < n_active; i++) s += h_sp[i];
                total_grad_spring_K += s;
            }
            if (use_pair) {
                std::vector<float> h_vp((size_t)n_active);
                CUDA_CHECK(cudaMemcpy(h_vp.data(), d_ViscPart, s_b,
                                      cudaMemcpyDeviceToHost));
                float s = 0.0f;
                for (uint32_t i = 0; i < n_active; i++) s += h_vp[i];
                total_grad_visc_K += s;
            }
        }

        // ── ρ gradient: kernel partial (per-particle) + implicit (via grad_C).
        //   implicit = -⟨grad_grad_C, grad_C⟩ / ρ_rest (per particle, summed).
        // Mirrors Metal lines 1363-1375.
        {
            std::vector<float> h_grho((size_t)n_active);
            std::vector<float> h_ggC ((size_t)n_active * 3);
            CUDA_CHECK(cudaMemcpy(h_grho.data(), d_Grho_per, s_b,
                                  cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(h_ggC.data(),  d_GgC,      pos_b,
                                  cudaMemcpyDeviceToHost));
            float kernel_rho = 0.0f, implicit_rho = 0.0f;
            for (uint32_t i = 0; i < n_active; i++) kernel_rho += h_grho[i];
            for (uint32_t i = 0; i < n_active; i++) {
                float dot = 0.0f;
                for (int ax = 0; ax < 3; ax++) {
                    dot += h_ggC[i * 3 + ax] * h_gradC[i * 3 + ax];
                }
                implicit_rho -= dot / rho_rest;
            }
            total_grad_rho += kernel_rho + implicit_rho;
        }

        // ── α gradient (per-particle ∂L/∂A from density solve backward).
        {
            std::vector<float> h_galp((size_t)n_active);
            CUDA_CHECK(cudaMemcpy(h_galp.data(), d_Galpha_per, s_b,
                                  cudaMemcpyDeviceToHost));
            for (uint32_t i = 0; i < n_active; i++) {
                total_grad_alpha_inv_dt2 += h_galp[i];
            }
        }

        // ── Per-step gradient clipping (BWD_CLIP_NORM) ──────────────
        // Compute the joint L2 norm sqrt(||grad_x||² + ||grad_v||²) and
        // scale both Gx_old_new and Gv_old_new by clip_norm / norm if it
        // exceeds clip_norm. Also zero NaN/Inf elements (degenerate-contact
        // particles can produce them despite kernel epsilons). Applied
        // BEFORE Gx_old_new/Gv_old_new get promoted to running gradients.
        if (clip_on) {
            std::vector<float> h_gx((size_t)n_active * 3);
            std::vector<float> h_gv((size_t)n_active * 3);
            CUDA_CHECK(cudaMemcpy(h_gx.data(), d_Gx_old_new, pos_b,
                                  cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(h_gv.data(), d_Gv_old_new, pos_b,
                                  cudaMemcpyDeviceToHost));
            size_t n_floats = 3u * (size_t)n_active;
            int n_bad_x = 0, n_bad_v = 0;
            for (size_t t = 0; t < n_floats; t++) {
                if (!std::isfinite(h_gx[t])) { h_gx[t] = 0.0f; n_bad_x++; }
                if (!std::isfinite(h_gv[t])) { h_gv[t] = 0.0f; n_bad_v++; }
            }
            // Independent grad_x and grad_v clipping (matches Metal at
            // src/metal_diff/ops_xpbd_full.mm:1397-1411).
            double sum_x = 0.0, sum_v = 0.0;
            for (size_t t = 0; t < n_floats; t++) {
                sum_x += (double)h_gx[t] * h_gx[t];
                sum_v += (double)h_gv[t] * h_gv[t];
            }
            float nx = (float)std::sqrt(sum_x);
            float nv = (float)std::sqrt(sum_v);
            if (nx > clip_norm_x && nx > 0.0f) {
                float s = clip_norm_x / nx;
                for (size_t t = 0; t < n_floats; t++) h_gx[t] *= s;
            }
            if (nv > clip_norm_v && nv > 0.0f) {
                float s = clip_norm_v / nv;
                for (size_t t = 0; t < n_floats; t++) h_gv[t] *= s;
            }
            CUDA_CHECK(cudaMemcpy(d_Gx_old_new, h_gx.data(), pos_b,
                                  cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_Gv_old_new, h_gv.data(), pos_b,
                                  cudaMemcpyHostToDevice));
            if ((n_bad_x || n_bad_v) && k % 50 == 0) {
                fprintf(stderr, "[xpbd_full_bwd k=%d] zeroed NaN/Inf: "
                                "%d in grad_x, %d in grad_v\n",
                        (int)k, n_bad_x, n_bad_v);
            }
        }

        // Promote per-step grads to running for previous step.
        CUDA_CHECK(cudaMemcpy(d_Gx_running, d_Gx_old_new, pos_b,
                              cudaMemcpyDeviceToDevice));
        CUDA_CHECK(cudaMemcpy(d_Gv_running, d_Gv_old_new, pos_b,
                              cudaMemcpyDeviceToDevice));
    }

    // ── Output files ─────────────────────────────────────────────────
    std::vector<float> h_gxi((size_t)n_active * 3);
    std::vector<float> h_gvi((size_t)n_active * 3);
    CUDA_CHECK(cudaMemcpy(h_gxi.data(), d_Gx_running, pos_b,
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_gvi.data(), d_Gv_running, pos_b,
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[14], h_gxi.data(), (size_t)n_active * 3);
    write_floats_or_die(argv[15], h_gvi.data(), (size_t)n_active * 3);
    write_floats_or_die(argv[16], &total_grad_rho, 1);
    // Optional param-grad outputs. Match Metal convention: write only if the
    // feature is enabled AND the caller provided the path.
    if (use_springs && grad_spring_K_out_path) {
        write_floats_or_die(grad_spring_K_out_path, &total_grad_spring_K, 1);
    }
    if (use_pair && grad_visc_K_out_path) {
        write_floats_or_die(grad_visc_K_out_path, &total_grad_visc_K, 1);
    }
    if (use_floor && grad_alpha_dens_out_path) {
        // ∂L/∂alpha_dens = ∂L/∂(alpha_inv_dt2) · ∂(alpha_inv_dt2)/∂alpha_dens
        //                = total_grad_alpha_inv_dt2 · (1/dt²)
        float total_grad_alpha_dens = total_grad_alpha_inv_dt2 / (dt * dt);
        write_floats_or_die(grad_alpha_dens_out_path, &total_grad_alpha_dens, 1);
    }

    // ── Cleanup ──────────────────────────────────────────────────────
    cudaFree(d_X); cudaFree(d_V); cudaFree(d_Xp);
    if (d_Sp) cudaFree(d_Sp);
    cudaFree(d_R2aa); cudaFree(d_GW_aa);
    if (d_R2as)  cudaFree(d_R2as);
    if (d_GW_as) cudaFree(d_GW_as);
    cudaFree(d_D); cudaFree(d_Gc); cudaFree(d_Dh); cudaFree(d_Lam);
    cudaFree(d_Gv_in); cudaFree(d_Gx_pred);
    cudaFree(d_Gx_old_new); cudaFree(d_Gv_old_new);
    cudaFree(d_Gx_running); cudaFree(d_Gv_running);
    cudaFree(d_Glam_pre); cudaFree(d_Gdens);
    cudaFree(d_GgC); cudaFree(d_Gdh);
    cudaFree(d_Grho_per); cudaFree(d_Galpha_per);
    cudaFree(d_Ggy_per);
    if (d_ExtA) cudaFree(d_ExtA);
    if (d_Gext) cudaFree(d_Gext);
    if (d_BondIJ)   cudaFree(d_BondIJ);
    if (d_BondRest) cudaFree(d_BondRest);
    if (d_SortedS)   cudaFree(d_SortedS);
    if (d_CellStart) cudaFree(d_CellStart);
    if (d_SpringPart) cudaFree(d_SpringPart);
    if (d_ViscPart)   cudaFree(d_ViscPart);
    if (d_PosPredSaved) cudaFree(d_PosPredSaved);
    if (d_FloorYGrad)   cudaFree(d_FloorYGrad);
    if (d_RestitGrad)   cudaFree(d_RestitGrad);

    std::free(all);
    if (pos_static) std::free(pos_static);
    std::free(grad_x_fin);
    if (sorted_static_buf) std::free(sorted_static_buf);
    if (cell_start_buf)    std::free(cell_start_buf);
    if (bond_ij_data)   std::free(bond_ij_data);
    if (bond_rest_data) std::free(bond_rest_data);
    return 0;
}
