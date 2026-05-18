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

// ops_xpbd_step.cu — M7 imperative XPBD orchestrator (xpbd_step op).
// Kernels in shaders.cu.

#include "cuda_common.h"
#include "shaders.h"
#include "ops.h"

// ──────────────────────────────────────────────────────────────────────
int run_xpbd_step(int argc, char **argv) {
    if (argc < 20) {
        fprintf(stderr,
            "usage: sib_cuda xpbd_step "
            "<n_active> <n_static> <h> <mass> <rho_rest> <dt> <gravity_y> "
            "<floor_y> <alpha_density> <n_iters> <n_steps> "
            "<pos_active.bin> <vel_active.bin> <pos_static.bin> "
            "<n_bonds> <bonds.bin> <alpha_dist> "
            "<pos_out.bin> <vel_out.bin> "
            "[traj_out.txt] [traj_every] [sim_scale] [vel_clamp] "
            "[restitution] [skip_density]\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_static = (unsigned int)std::atoi(argv[3]);
    float h           = (float)std::atof(argv[4]);
    float mass        = (float)std::atof(argv[5]);
    float rho_rest    = (float)std::atof(argv[6]);
    float dt          = (float)std::atof(argv[7]);
    float gravity_y   = (float)std::atof(argv[8]);
    float floor_y     = (float)std::atof(argv[9]);
    float alpha_dens  = (float)std::atof(argv[10]);
    unsigned int n_iters  = (unsigned int)std::atoi(argv[11]);
    unsigned int n_steps  = (unsigned int)std::atoi(argv[12]);
    const char *path_pos_a = argv[13];
    const char *path_vel_a = argv[14];
    const char *path_pos_s = argv[15];
    unsigned int n_bonds = (unsigned int)std::atoi(argv[16]);
    const char *path_bonds = argv[17];
    float alpha_dist = (float)std::atof(argv[18]);
    const char *path_pos_out = argv[19];
    const char *path_vel_out = argv[20];
    const char *path_traj   = (argc > 21) ? argv[21] : nullptr;
    int traj_every          = (argc > 22) ? std::atoi(argv[22]) : 0;
    float sim_scale         = (argc > 23) ? (float)std::atof(argv[23]) : 1.0f;
    float vel_clamp_v       = (argc > 24) ? (float)std::atof(argv[24]) : 0.0f;
    float restitution       = (argc > 25) ? (float)std::atof(argv[25]) : 0.0f;
    int skip_density        = (argc > 26) ? std::atoi(argv[26]) : 0;
    float sim_scale_inv = 1.0f / sim_scale;

    // Compliance constants for XPBD's α/dt² term.
    float inv_dt2 = 1.0f / (dt * dt);
    float alpha_dens_inv_dt2 = alpha_dens * inv_dt2;
    float alpha_dist_inv_dt2 = alpha_dist * inv_dt2;
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));
    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));

    // Host-side inputs.
    float *h_pos_a = read_floats_or_die(path_pos_a, (size_t)n_active * 3);
    float *h_vel_a = read_floats_or_die(path_vel_a, (size_t)n_active * 3);
    float *h_pos_s = (n_static > 0)
        ? read_floats_or_die(path_pos_s, (size_t)n_static * 3)
        : nullptr;
    // bonds.bin: 16 bytes/bond → unpack into int2 + float rest_len.
    int2 *h_bond_ij = nullptr;
    float *h_rest_len = nullptr;
    if (n_bonds > 0) {
        FILE *f = std::fopen(path_bonds, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path_bonds); return 1; }
        h_bond_ij = (int2 *)std::malloc((size_t)n_bonds * sizeof(int2));
        h_rest_len = (float *)std::malloc((size_t)n_bonds * sizeof(float));
        for (unsigned int b = 0; b < n_bonds; b++) {
            int32_t pair[2];
            float lenpad[2];
            if (std::fread(pair, sizeof(int32_t), 2, f) != 2 ||
                std::fread(lenpad, sizeof(float), 2, f) != 2) {
                fprintf(stderr, "short read on %s\n", path_bonds);
                return 1;
            }
            h_bond_ij[b].x = pair[0];
            h_bond_ij[b].y = pair[1];
            h_rest_len[b] = lenpad[0];
        }
        std::fclose(f);
    }

    // Device-side allocations.
    float3 *d_pos_old = nullptr, *d_pos_pred = nullptr;
    float3 *d_vel = nullptr, *d_pos_s = nullptr;
    float *d_r2_aa = nullptr, *d_r2_as = nullptr;
    float *d_W_aa = nullptr, *d_W_as = nullptr;
    float *d_density_aa = nullptr, *d_density_as = nullptr;
    float3 *d_grad_C = nullptr;
    float *d_denom_helper = nullptr;
    float *d_lambda_dens = nullptr;
    int2 *d_bond_ij = nullptr;
    float *d_rest_len = nullptr;
    float *d_lambda_dist = nullptr;

    size_t pos_bytes_a = (size_t)n_active * sizeof(float3);
    CUDA_CHECK(cudaMalloc(&d_pos_old,  pos_bytes_a));
    CUDA_CHECK(cudaMalloc(&d_pos_pred, pos_bytes_a));
    CUDA_CHECK(cudaMalloc(&d_vel,      pos_bytes_a));
    CUDA_CHECK(cudaMemcpy(d_pos_old, h_pos_a, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vel,     h_vel_a, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));

    if (n_static > 0) {
        CUDA_CHECK(cudaMalloc(&d_pos_s, (size_t)n_static * sizeof(float3)));
        CUDA_CHECK(cudaMemcpy(d_pos_s, h_pos_s,
                              (size_t)n_static * 3 * sizeof(float),
                              cudaMemcpyHostToDevice));
    }
    if (!skip_density) {
        CUDA_CHECK(cudaMalloc(&d_r2_aa, (size_t)n_active * n_active * sizeof(float)));
        CUDA_CHECK(cudaMalloc(&d_W_aa,  (size_t)n_active * n_active * sizeof(float)));
        CUDA_CHECK(cudaMalloc(&d_density_aa, (size_t)n_active * sizeof(float)));
        CUDA_CHECK(cudaMalloc(&d_grad_C, (size_t)n_active * sizeof(float3)));
        CUDA_CHECK(cudaMalloc(&d_denom_helper, (size_t)n_active * sizeof(float)));
        CUDA_CHECK(cudaMalloc(&d_lambda_dens, (size_t)n_active * sizeof(float)));
        CUDA_CHECK(cudaMemset(d_lambda_dens, 0, (size_t)n_active * sizeof(float)));
        if (n_static > 0) {
            CUDA_CHECK(cudaMalloc(&d_r2_as, (size_t)n_active * n_static * sizeof(float)));
            CUDA_CHECK(cudaMalloc(&d_W_as,  (size_t)n_active * n_static * sizeof(float)));
            CUDA_CHECK(cudaMalloc(&d_density_as, (size_t)n_active * sizeof(float)));
        }
    }
    if (n_bonds > 0) {
        CUDA_CHECK(cudaMalloc(&d_bond_ij, (size_t)n_bonds * sizeof(int2)));
        CUDA_CHECK(cudaMalloc(&d_rest_len, (size_t)n_bonds * sizeof(float)));
        CUDA_CHECK(cudaMalloc(&d_lambda_dist, (size_t)n_bonds * sizeof(float)));
        CUDA_CHECK(cudaMemcpy(d_bond_ij, h_bond_ij,
                              (size_t)n_bonds * sizeof(int2),
                              cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_rest_len, h_rest_len,
                              (size_t)n_bonds * sizeof(float),
                              cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemset(d_lambda_dist, 0, (size_t)n_bonds * sizeof(float)));
    }

    const unsigned int TPB = 256;
    unsigned int gridA = (n_active + TPB - 1) / TPB;

    // Optional trajectory writer (Sibernetic position_buffer.txt format).
    FILE *traj = nullptr;
    if (path_traj && traj_every > 0) {
        traj = std::fopen(path_traj, "w");
        if (!traj) { fprintf(stderr, "cannot open %s\n", path_traj); return 1; }
        std::fprintf(traj, "0\n0\n0\n0\n0\n0\n");  // box bounds placeholder
        std::fprintf(traj, "%u\n", n_active);       // numOfElasticP
        std::fprintf(traj, "0\n");                  // numOfLiquidP
        std::fprintf(traj, "%u\n", n_static);       // numOfBoundaryP
        std::fprintf(traj, "%g\n", dt);             // time_step
        std::fprintf(traj, "%d\n", traj_every);     // log_step
    }
    auto dump_frame = [&](float3 *d_pos) {
        if (!traj) return;
        std::vector<float> tmp(n_active * 3);
        CUDA_CHECK(cudaMemcpy(tmp.data(), d_pos,
                              (size_t)n_active * 3 * sizeof(float),
                              cudaMemcpyDeviceToHost));
        for (unsigned int k = 0; k < n_active; k++) {
            std::fprintf(traj, "%g %g %g 2\n",
                         tmp[3*k], tmp[3*k+1], tmp[3*k+2]);
        }
        if (n_static > 0) {
            for (unsigned int k = 0; k < n_static; k++) {
                std::fprintf(traj, "%g %g %g 3\n",
                             h_pos_s[3*k], h_pos_s[3*k+1], h_pos_s[3*k+2]);
            }
        }
    };

    dump_frame(d_pos_old);  // frame 0

    auto t0 = std::chrono::steady_clock::now();
    for (unsigned int step = 0; step < n_steps; step++) {
        // M7.0 predict
        predict_positions<<<gridA, TPB>>>(d_pos_old, d_vel, d_pos_pred,
                                          dt, gravity_y, n_active,
                                          sim_scale_inv);

        // Reset accumulated Lagrange multipliers each step (XPBD convention).
        if (!skip_density)
            CUDA_CHECK(cudaMemset(d_lambda_dens, 0,
                                  (size_t)n_active * sizeof(float)));
        if (n_bonds > 0)
            CUDA_CHECK(cudaMemset(d_lambda_dist, 0,
                                  (size_t)n_bonds * sizeof(float)));

        for (unsigned int it = 0; it < n_iters; it++) {
            if (!skip_density) {
                // rowsum_density / density_constraint_grad declare
                // __shared__ float partials[256] and tree-reduce with
                // stride = T/2. TPB > 256 writes OOB; TPB not a power of 2
                // <= 256 silently drops elements.
                constexpr unsigned int TPB_REDUCE = 256;
                static_assert(TPB_REDUCE == 256,
                              "rowsum_density/density_constraint_grad: "
                              "TPB must be 256 (shared partials[256] + "
                              "power-of-2 tree reduce)");
                // density chain
                {
                    dim3 t2(16, 16);
                    dim3 b2((n_active + 15) / 16, (n_active + 15) / 16);
                    dist_active_active<<<b2, t2>>>(d_pos_pred, d_r2_aa, n_active);
                }
                unsigned int n_aa = n_active * n_active;
                wpoly6_inplace<<<(n_aa + TPB - 1) / TPB, TPB>>>(
                    d_r2_aa, h2, poly6_const, n_aa);
                // Default stream: synchronous w.r.t. subsequent launches
                // on the same stream. Plain cudaMemcpy is the honest
                // spelling (cudaMemcpyAsync without a non-default stream
                // is misleading and easy to miss the error code on).
                CUDA_CHECK(cudaMemcpy(d_W_aa, d_r2_aa, n_aa * sizeof(float),
                                      cudaMemcpyDeviceToDevice));
                rowsum_density<<<n_active, TPB_REDUCE>>>(d_W_aa, d_density_aa,
                                                         mass,
                                                         n_active, n_active);
                if (n_static > 0) {
                    dim3 t2(16, 16);
                    dim3 b2((n_active + 15) / 16, (n_static + 15) / 16);
                    dist_active_static<<<b2, t2>>>(d_pos_pred, d_pos_s,
                                                   d_r2_as, n_active, n_static);
                    unsigned int n_as = n_active * n_static;
                    wpoly6_inplace<<<(n_as + TPB - 1) / TPB, TPB>>>(
                        d_r2_as, h2, poly6_const, n_as);
                    CUDA_CHECK(cudaMemcpy(d_W_as, d_r2_as,
                                          n_as * sizeof(float),
                                          cudaMemcpyDeviceToDevice));
                    rowsum_density<<<n_active, TPB_REDUCE>>>(d_W_as,
                                                             d_density_as,
                                                             mass,
                                                             n_static,
                                                             n_active);
                    add_inplace<<<gridA, TPB>>>(d_density_aa, d_density_as,
                                                n_active);
                }
                density_constraint_grad<<<n_active, TPB_REDUCE>>>(
                    d_pos_pred, d_pos_s, d_r2_aa, d_r2_as,
                    d_grad_C, d_denom_helper,
                    h, spiky_const, mass, rho_rest, n_active, n_static);
                solve_density_constraint<<<gridA, TPB>>>(
                    d_pos_pred, d_lambda_dens, d_density_aa,
                    d_grad_C, d_denom_helper,
                    rho_rest, mass, alpha_dens_inv_dt2, n_active);
            }
            if (n_bonds > 0) {
                solve_distance_constraints_seq<<<1, 1>>>(
                    d_pos_pred, d_lambda_dist, d_bond_ij, d_rest_len,
                    alpha_dist_inv_dt2, 1.0f / mass, n_bonds);
            }
            solve_floor_constraint<<<gridA, TPB>>>(d_pos_pred, floor_y,
                                                   n_active, restitution);
        }

        update_velocities<<<gridA, TPB>>>(d_vel, d_pos_old, d_pos_pred,
                                          dt, n_active, sim_scale);
        if (vel_clamp_v > 0.0f) {
            clamp_velocity<<<gridA, TPB>>>(d_vel, vel_clamp_v, n_active);
        }

        // swap: pos_old <- pos_pred (avoid extra copy by ping-pong'ing).
        std::swap(d_pos_old, d_pos_pred);

        if (traj_every > 0 && ((step + 1) % (unsigned)traj_every == 0)) {
            dump_frame(d_pos_old);
        }
    }
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    auto t1 = std::chrono::steady_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    fprintf(stderr, "xpbd_step: %u steps in %.1f ms (%.3f ms/step)\n",
            n_steps, total_ms, total_ms / std::max(1u, n_steps));

    // Write final pos/vel.
    float *out_pos = (float *)std::malloc((size_t)n_active * 3 * sizeof(float));
    float *out_vel = (float *)std::malloc((size_t)n_active * 3 * sizeof(float));
    CUDA_CHECK(cudaMemcpy(out_pos, d_pos_old,
                          (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(out_vel, d_vel,
                          (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(path_pos_out, out_pos, (size_t)n_active * 3);
    write_floats_or_die(path_vel_out, out_vel, (size_t)n_active * 3);

    if (traj) std::fclose(traj);

    cudaFree(d_pos_old); cudaFree(d_pos_pred); cudaFree(d_vel);
    if (d_pos_s) cudaFree(d_pos_s);
    if (d_r2_aa) cudaFree(d_r2_aa); if (d_r2_as) cudaFree(d_r2_as);
    if (d_W_aa) cudaFree(d_W_aa);   if (d_W_as) cudaFree(d_W_as);
    if (d_density_aa) cudaFree(d_density_aa);
    if (d_density_as) cudaFree(d_density_as);
    if (d_grad_C) cudaFree(d_grad_C);
    if (d_denom_helper) cudaFree(d_denom_helper);
    if (d_lambda_dens) cudaFree(d_lambda_dens);
    if (d_bond_ij) cudaFree(d_bond_ij);
    if (d_rest_len) cudaFree(d_rest_len);
    if (d_lambda_dist) cudaFree(d_lambda_dist);
    std::free(h_pos_a); std::free(h_vel_a);
    if (h_pos_s) std::free(h_pos_s);
    if (h_bond_ij) std::free(h_bond_ij);
    if (h_rest_len) std::free(h_rest_len);
    std::free(out_pos); std::free(out_vel);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// xpbd_step_diff_fwd / xpbd_step_diff_bwd
//
// Multi-step differentiable cube-drop pipeline:
//   predict -> (bonds + floor) × n_iters -> update_velocities
//
// Density / pair-forces / membranes are intentionally OUT of scope for
// this driver — they require their own backward kernels and substantial
// state. This is the differentiable cube-drop substrate that proves the
// multi-step gradient chain works end-to-end. Density follows once its
// backward chain lands.
//
// Tape layout per step:
//   offset                            content
//   ─────────────────────────────────────────────────────────────────
//   0                                 pos_old[3N]
//   3N                                pos_pred[3N]
//   for m in 0..n_iters-1:
//     6N + m*(7B + 3N)               bond_state[7B]
//     6N + m*(7B + 3N) + 7B          pos_pre_floor[3N]
//
// step_stride = 6N + n_iters * (7B + 3N) floats
// total tape  = n_steps * step_stride floats
// ──────────────────────────────────────────────────────────────────────
