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

// ops_pair_spring.cu — CLI drivers for the worm-physics kernel family.
//
// Kernels live in shaders.cu; this TU just handles I/O, device alloc,
// kernel launch, and writing results back. Bond format matches the rest
// of the codebase (16 bytes/bond: int i, int j, float rest_len, float pad).
// Anchor format is 32 bytes/anchor (8 floats):
//   slot 0: int particle index (bit-reinterpret of a float)
//   slot 1: pad
//   slot 2-4: world anchor position (x, y, z)
//   slot 5: rest length
//   slot 6-7: pad
//
// All "_inout" buffers expect the caller to pre-zero them (or pre-load
// with prior accumulated grads). All "_out" buffers are fresh outputs.

#include "cuda_common.h"
#include "shaders.h"
#include "ops.h"

// Helper: load the 16-byte/bond binary into int2 + rest_len arrays.
static void load_bonds(const char *path, unsigned int n_bonds,
                       int2 **out_ij, float **out_rest)
{
    FILE *f = std::fopen(path, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); std::exit(1); }
    *out_ij   = (int2 *)std::malloc((size_t)n_bonds * sizeof(int2));
    *out_rest = (float *)std::malloc((size_t)n_bonds * sizeof(float));
    for (unsigned int b = 0; b < n_bonds; b++) {
        int32_t pair[2]; float lenpad[2];
        if (std::fread(pair, sizeof(int32_t), 2, f) != 2 ||
            std::fread(lenpad, sizeof(float), 2, f) != 2) {
            fprintf(stderr, "short read on %s\n", path); std::exit(1);
        }
        (*out_ij)[b] = make_int2(pair[0], pair[1]);
        (*out_rest)[b] = lenpad[0];
    }
    std::fclose(f);
}

// Pair-force amplitude scalars (mirrors metal_diff/ops_pair_spring.mm:60-67).
// Computed host-side in double, downcast to float, then passed to the
// kernel. visc_amp drives the viscosity pair force; surf_amp drives the
// surface-tension cohesion force. Forward, backward, and visc_K_partial
// drivers (plus the integrated xpbd_full pipeline) all share this
// formula via the declaration in cuda_common.h.
void compute_pair_amps(float h, float sim_scale, float mass,
                       float *out_h2,
                       float *out_visc_amp, float *out_surf_amp)
{
    *out_h2 = h * h;
    double h_scaled = (double)h * (double)sim_scale;
    double h_s6 = std::pow(h_scaled, 6.0);
    double h_s9 = std::pow(h_scaled, 9.0);
    double divgradWvisco = 45.0 / (M_PI * h_s6);
    *out_visc_amp = (float)(1.5 * (double)mass * divgradWvisco *
                            std::pow((double)sim_scale, 3.0));
    double wpoly6_si = 315.0 / (64.0 * M_PI * h_s9);
    // mass cancels in surf_amp; kept as written to mirror the OpenCL /
    // Metal formula bit-for-bit (sphFluid.cl computes ·mass on the inner
    // sum and divides by ·mass when applying — preserving the same
    // grouping here keeps fp64-rounded parity with the other backends).
    *out_surf_amp = (float)(-1.7e-9 * (double)mass * wpoly6_si *
                            (double)sim_scale / (double)mass);
}

// ── spring_bonds_force ────────────────────────────────────────────────
int run_spring_bonds_force(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr,
            "usage: sib_cuda spring_bonds_force "
            "<n_active> <n_bonds> <spring_K> "
            "<active_pos.bin> <bonds.bin> <ext_accel_inout.bin>\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_bonds  = (unsigned int)std::atoi(argv[3]);
    float spring_K = (float)std::atof(argv[4]);

    float *h_pos = read_floats_or_die(argv[5], (size_t)n_active * 3);
    int2 *h_ij = nullptr; float *h_rest = nullptr;
    load_bonds(argv[6], n_bonds, &h_ij, &h_rest);
    float *h_ea = read_floats_or_die(argv[7], (size_t)n_active * 3);

    float3 *d_pos = nullptr, *d_ea = nullptr;
    int2 *d_ij = nullptr; float *d_rest = nullptr;
    CUDA_CHECK(cudaMalloc(&d_pos,  (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_ea,   (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_ij,   (size_t)n_bonds * sizeof(int2)));
    CUDA_CHECK(cudaMalloc(&d_rest, (size_t)n_bonds * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_pos, h_pos, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ea,  h_ea,  (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ij,  h_ij,  (size_t)n_bonds * sizeof(int2),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_rest, h_rest, (size_t)n_bonds * sizeof(float),
                          cudaMemcpyHostToDevice));

    const unsigned int TPB = 128;
    unsigned int grid = (n_active + TPB - 1) / TPB;
    spring_bonds_force<<<grid, TPB>>>(d_pos, d_ij, d_rest, d_ea,
                                      spring_K, n_bonds, n_active);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_ea, d_ea, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[7], h_ea, (size_t)n_active * 3);

    cudaFree(d_pos); cudaFree(d_ea); cudaFree(d_ij); cudaFree(d_rest);
    std::free(h_pos); std::free(h_ea); std::free(h_ij); std::free(h_rest);
    return 0;
}

// ── apply_ext_accel ───────────────────────────────────────────────────
int run_apply_ext_accel(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr,
            "usage: sib_cuda apply_ext_accel "
            "<n> <dt> <vel_inout.bin> <ext_accel.bin>\n");
        return 1;
    }
    unsigned int n = (unsigned int)std::atoi(argv[2]);
    float dt = (float)std::atof(argv[3]);

    float *h_vel = read_floats_or_die(argv[4], (size_t)n * 3);
    float *h_ea  = read_floats_or_die(argv[5], (size_t)n * 3);

    float3 *d_vel = nullptr, *d_ea = nullptr;
    CUDA_CHECK(cudaMalloc(&d_vel, (size_t)n * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_ea,  (size_t)n * sizeof(float3)));
    CUDA_CHECK(cudaMemcpy(d_vel, h_vel, (size_t)n * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ea,  h_ea,  (size_t)n * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));

    const unsigned int TPB = 256;
    unsigned int grid = (n + TPB - 1) / TPB;
    apply_ext_accel<<<grid, TPB>>>(d_vel, d_ea, dt, n);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_vel, d_vel, (size_t)n * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[4], h_vel, (size_t)n * 3);

    cudaFree(d_vel); cudaFree(d_ea);
    std::free(h_vel); std::free(h_ea);
    return 0;
}

// ── spring_bonds_force_backward ───────────────────────────────────────
int run_spring_bonds_force_bwd(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr,
            "usage: sib_cuda spring_bonds_force_bwd "
            "<n_active> <n_bonds> <spring_K> "
            "<active_pos.bin> <bonds.bin> <grad_ext_accel.bin> "
            "<grad_pos_inout.bin>\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_bonds  = (unsigned int)std::atoi(argv[3]);
    float spring_K = (float)std::atof(argv[4]);

    float *h_pos = read_floats_or_die(argv[5], (size_t)n_active * 3);
    int2 *h_ij = nullptr; float *h_rest = nullptr;
    load_bonds(argv[6], n_bonds, &h_ij, &h_rest);
    float *h_gea = read_floats_or_die(argv[7], (size_t)n_active * 3);
    float *h_gp  = read_floats_or_die(argv[8], (size_t)n_active * 3);

    float3 *d_pos = nullptr, *d_gea = nullptr, *d_gp = nullptr;
    int2 *d_ij = nullptr; float *d_rest = nullptr;
    CUDA_CHECK(cudaMalloc(&d_pos,  (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gea,  (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gp,   (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_ij,   (size_t)n_bonds * sizeof(int2)));
    CUDA_CHECK(cudaMalloc(&d_rest, (size_t)n_bonds * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_pos, h_pos, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gea, h_gea, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gp,  h_gp,  (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ij,  h_ij,  (size_t)n_bonds * sizeof(int2),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_rest, h_rest, (size_t)n_bonds * sizeof(float),
                          cudaMemcpyHostToDevice));

    const unsigned int TPB = 128;
    unsigned int grid = (n_active + TPB - 1) / TPB;
    spring_bonds_force_backward<<<grid, TPB>>>(d_pos, d_ij, d_rest, d_gea,
                                               d_gp, spring_K, n_bonds, n_active);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_gp, d_gp, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[8], h_gp, (size_t)n_active * 3);

    cudaFree(d_pos); cudaFree(d_gea); cudaFree(d_gp);
    cudaFree(d_ij); cudaFree(d_rest);
    std::free(h_pos); std::free(h_gea); std::free(h_gp);
    std::free(h_ij); std::free(h_rest);
    return 0;
}

// ── spring_K_partial ──────────────────────────────────────────────────
// Outputs per-particle scalar; host sums to get ∂L/∂spring_K.
int run_spring_K_partial(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr,
            "usage: sib_cuda spring_K_partial "
            "<n_active> <n_bonds> "
            "<active_pos.bin> <bonds.bin> <grad_ext_accel.bin> "
            "<per_particle_out.bin>\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_bonds  = (unsigned int)std::atoi(argv[3]);

    float *h_pos = read_floats_or_die(argv[4], (size_t)n_active * 3);
    int2 *h_ij = nullptr; float *h_rest = nullptr;
    load_bonds(argv[5], n_bonds, &h_ij, &h_rest);
    float *h_gea = read_floats_or_die(argv[6], (size_t)n_active * 3);

    float3 *d_pos = nullptr, *d_gea = nullptr;
    int2 *d_ij = nullptr; float *d_rest = nullptr, *d_pp = nullptr;
    CUDA_CHECK(cudaMalloc(&d_pos,  (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gea,  (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_ij,   (size_t)n_bonds * sizeof(int2)));
    CUDA_CHECK(cudaMalloc(&d_rest, (size_t)n_bonds * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_pp,   (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_pos, h_pos, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gea, h_gea, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ij,  h_ij,  (size_t)n_bonds * sizeof(int2),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_rest, h_rest, (size_t)n_bonds * sizeof(float),
                          cudaMemcpyHostToDevice));

    const unsigned int TPB = 128;
    unsigned int grid = (n_active + TPB - 1) / TPB;
    spring_K_partial<<<grid, TPB>>>(d_pos, d_ij, d_rest, d_gea, d_pp,
                                    n_bonds, n_active);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float *h_pp = (float *)std::malloc((size_t)n_active * sizeof(float));
    CUDA_CHECK(cudaMemcpy(h_pp, d_pp, (size_t)n_active * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[7], h_pp, (size_t)n_active);

    cudaFree(d_pos); cudaFree(d_gea); cudaFree(d_ij);
    cudaFree(d_rest); cudaFree(d_pp);
    std::free(h_pos); std::free(h_gea); std::free(h_ij);
    std::free(h_rest); std::free(h_pp);
    return 0;
}

// ── pair_forces_grid_fwd ──────────────────────────────────────────────
// Reads cell_start (n_cells+1 int32) and grid_origin (3 float32) from
// explicit file paths — Windows has no /tmp so the Metal hardcoded path
// trick doesn't carry over. Computes visc_amp and surf_amp host-side
// (in double, downcast to float) using the same formulas as the Metal
// driver, and passes those scalars through to the kernel.
static int *read_ints_or_die(const char *path, size_t n_ints) {
    FILE *f = std::fopen(path, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); std::exit(1); }
    int *buf = (int *)std::malloc(n_ints * sizeof(int));
    if (std::fread(buf, sizeof(int), n_ints, f) != n_ints) {
        fprintf(stderr, "short read on %s\n", path); std::exit(1);
    }
    std::fclose(f);
    return buf;
}

int run_pair_forces_grid_fwd(int argc, char **argv) {
    if (argc != 18) {
        fprintf(stderr,
            "usage: sib_cuda pair_forces_grid_fwd "
            "<n_active> <n_static> <h> <mass> <sim_scale> <visc_pair_coef> "
            "<pos_active.bin> <vel_active.bin> <sorted_static.bin> "
            "<cell_start.bin> <density.bin> "
            "<grid_dim_x> <grid_dim_y> <grid_dim_z> "
            "<grid_origin.bin> <ext_accel_out.bin>\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_static = (unsigned int)std::atoi(argv[3]);
    float h              = (float)std::atof(argv[4]);
    float mass           = (float)std::atof(argv[5]);
    float sim_scale      = (float)std::atof(argv[6]);
    float visc_pair_coef = (float)std::atof(argv[7]);
    const char *path_pa = argv[8];
    const char *path_va = argv[9];
    const char *path_ss = argv[10];
    const char *path_cs = argv[11];
    const char *path_d  = argv[12];
    int grid_dim_x = std::atoi(argv[13]);
    int grid_dim_y = std::atoi(argv[14]);
    int grid_dim_z = std::atoi(argv[15]);
    const char *path_go = argv[16];
    const char *path_out = argv[17];

    float h2, visc_amp, surf_amp;
    compute_pair_amps(h, sim_scale, mass, &h2, &visc_amp, &surf_amp);

    int n_cells = grid_dim_x * grid_dim_y * grid_dim_z;
    float *h_pos = read_floats_or_die(path_pa, (size_t)n_active * 3);
    float *h_vel = read_floats_or_die(path_va, (size_t)n_active * 3);
    float *h_ss  = read_floats_or_die(path_ss, (size_t)n_static * 3);
    int   *h_cs  = read_ints_or_die  (path_cs, (size_t)(n_cells + 1));
    float *h_d   = read_floats_or_die(path_d,  (size_t)n_active);
    float *h_go  = read_floats_or_die(path_go, 3);

    float3 *d_pos = nullptr, *d_vel = nullptr, *d_ss = nullptr, *d_ea = nullptr;
    int *d_cs = nullptr;
    float *d_d = nullptr;
    CUDA_CHECK(cudaMalloc(&d_pos, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_vel, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_ss,  (size_t)n_static * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_ea,  (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_cs,  (size_t)(n_cells + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_d,   (size_t)n_active * sizeof(float)));

    CUDA_CHECK(cudaMemcpy(d_pos, h_pos, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vel, h_vel, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ss,  h_ss,  (size_t)n_static * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cs,  h_cs,  (size_t)(n_cells + 1) * sizeof(int),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_d,   h_d,   (size_t)n_active * sizeof(float),
                          cudaMemcpyHostToDevice));

    // 32 threads per particle, one block per particle (single warp).
    // shaders.cu hardcodes a 5-round __shfl_down_sync (offsets 16/8/4/2/1)
    // that sums exactly one warp; TPB_PAIR (cuda_common.h) is fixed at 32.
    pair_forces_grid_fwd<<<n_active, TPB_PAIR>>>(
        d_pos, d_vel, d_ss, d_cs, d_d, d_ea,
        h, h2, sim_scale, visc_pair_coef, visc_amp, surf_amp,
        n_active,
        grid_dim_x, grid_dim_y, grid_dim_z,
        h_go[0], h_go[1], h_go[2]);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float *h_ea = (float *)std::malloc((size_t)n_active * 3 * sizeof(float));
    CUDA_CHECK(cudaMemcpy(h_ea, d_ea, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(path_out, h_ea, (size_t)n_active * 3);

    cudaFree(d_pos); cudaFree(d_vel); cudaFree(d_ss);
    cudaFree(d_ea);  cudaFree(d_cs);  cudaFree(d_d);
    std::free(h_pos); std::free(h_vel); std::free(h_ss);
    std::free(h_cs);  std::free(h_d);   std::free(h_go);
    std::free(h_ea);
    return 0;
}

// ── pair_forces_grid_bwd ──────────────────────────────────────────────
// Mirrors run_pair_forces_grid_fwd's CLI convention (grid_origin path on
// argv, not hardcoded /tmp). Two output files: grad_pos and grad_vel.
// Both are pre-zeroed since the kernel does grad_pos[i] += dp_i.
int run_pair_forces_grid_bwd(int argc, char **argv) {
    if (argc != 20) {
        fprintf(stderr,
            "usage: sib_cuda pair_forces_grid_bwd "
            "<n_active> <n_static> <h> <mass> <sim_scale> <visc_pair_coef> "
            "<pos_active.bin> <vel_active.bin> <sorted_static.bin> "
            "<cell_start.bin> <density.bin> "
            "<grid_dim_x> <grid_dim_y> <grid_dim_z> "
            "<grid_origin.bin> <grad_ext_accel.bin> "
            "<grad_pos_out.bin> <grad_vel_out.bin>\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_static = (unsigned int)std::atoi(argv[3]);
    float h              = (float)std::atof(argv[4]);
    float mass           = (float)std::atof(argv[5]);
    float sim_scale      = (float)std::atof(argv[6]);
    float visc_pair_coef = (float)std::atof(argv[7]);
    const char *path_pa = argv[8];
    const char *path_va = argv[9];
    const char *path_ss = argv[10];
    const char *path_cs = argv[11];
    const char *path_d  = argv[12];
    int grid_dim_x = std::atoi(argv[13]);
    int grid_dim_y = std::atoi(argv[14]);
    int grid_dim_z = std::atoi(argv[15]);
    const char *path_go  = argv[16];
    const char *path_gea = argv[17];
    const char *path_gp_out = argv[18];
    const char *path_gv_out = argv[19];

    float h2, visc_amp, surf_amp;
    compute_pair_amps(h, sim_scale, mass, &h2, &visc_amp, &surf_amp);

    int n_cells = grid_dim_x * grid_dim_y * grid_dim_z;
    float *h_pos = read_floats_or_die(path_pa, (size_t)n_active * 3);
    float *h_vel = read_floats_or_die(path_va, (size_t)n_active * 3);
    float *h_ss  = read_floats_or_die(path_ss, (size_t)n_static * 3);
    int   *h_cs  = read_ints_or_die  (path_cs, (size_t)(n_cells + 1));
    float *h_d   = read_floats_or_die(path_d,  (size_t)n_active);
    float *h_go  = read_floats_or_die(path_go, 3);
    float *h_gea = read_floats_or_die(path_gea, (size_t)n_active * 3);

    float3 *d_pos = nullptr, *d_vel = nullptr, *d_ss = nullptr;
    float3 *d_gea = nullptr, *d_gp = nullptr, *d_gv = nullptr;
    int *d_cs = nullptr;
    float *d_d = nullptr;
    CUDA_CHECK(cudaMalloc(&d_pos, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_vel, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_ss,  (size_t)n_static * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gea, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gp,  (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gv,  (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_cs,  (size_t)(n_cells + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_d,   (size_t)n_active * sizeof(float)));

    CUDA_CHECK(cudaMemcpy(d_pos, h_pos, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vel, h_vel, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ss,  h_ss,  (size_t)n_static * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gea, h_gea, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cs,  h_cs,  (size_t)(n_cells + 1) * sizeof(int),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_d,   h_d,   (size_t)n_active * sizeof(float),
                          cudaMemcpyHostToDevice));
    // Kernel does grad_pos[i] += dp_i — pre-zero to match Metal driver.
    CUDA_CHECK(cudaMemset(d_gp, 0, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMemset(d_gv, 0, (size_t)n_active * sizeof(float3)));

    // One thread per particle (Metal: dispatchThreads(n_active, 1, 1) tpg 64).
    const unsigned int TPB = 64;
    unsigned int grid = (n_active + TPB - 1) / TPB;
    pair_forces_grid_bwd<<<grid, TPB>>>(
        d_pos, d_vel, d_ss, d_cs, d_d, d_gea, d_gp, d_gv,
        h, h2, sim_scale, visc_pair_coef, visc_amp, surf_amp,
        n_active,
        grid_dim_x, grid_dim_y, grid_dim_z,
        h_go[0], h_go[1], h_go[2]);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float *h_gp = (float *)std::malloc((size_t)n_active * 3 * sizeof(float));
    float *h_gv = (float *)std::malloc((size_t)n_active * 3 * sizeof(float));
    CUDA_CHECK(cudaMemcpy(h_gp, d_gp, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_gv, d_gv, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(path_gp_out, h_gp, (size_t)n_active * 3);
    write_floats_or_die(path_gv_out, h_gv, (size_t)n_active * 3);

    cudaFree(d_pos); cudaFree(d_vel); cudaFree(d_ss);
    cudaFree(d_gea); cudaFree(d_gp);  cudaFree(d_gv);
    cudaFree(d_cs);  cudaFree(d_d);
    std::free(h_pos); std::free(h_vel); std::free(h_ss);
    std::free(h_cs);  std::free(h_d);   std::free(h_go);
    std::free(h_gea); std::free(h_gp);  std::free(h_gv);
    return 0;
}

// ── apply_ext_accel_backward ──────────────────────────────────────────
int run_apply_ext_accel_bwd(int argc, char **argv) {
    if (argc != 7) {
        fprintf(stderr,
            "usage: sib_cuda apply_ext_accel_bwd "
            "<n> <dt> <grad_v_new.bin> "
            "<grad_v_old_inout.bin> <grad_ext_accel_inout.bin>\n");
        return 1;
    }
    unsigned int n = (unsigned int)std::atoi(argv[2]);
    float dt = (float)std::atof(argv[3]);

    float *h_gn = read_floats_or_die(argv[4], (size_t)n * 3);
    float *h_go = read_floats_or_die(argv[5], (size_t)n * 3);
    float *h_ge = read_floats_or_die(argv[6], (size_t)n * 3);

    float3 *d_gn = nullptr, *d_go = nullptr, *d_ge = nullptr;
    CUDA_CHECK(cudaMalloc(&d_gn, (size_t)n * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_go, (size_t)n * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_ge, (size_t)n * sizeof(float3)));
    CUDA_CHECK(cudaMemcpy(d_gn, h_gn, (size_t)n * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_go, h_go, (size_t)n * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ge, h_ge, (size_t)n * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));

    const unsigned int TPB = 256;
    unsigned int grid = (n + TPB - 1) / TPB;
    apply_ext_accel_backward<<<grid, TPB>>>(d_gn, d_go, d_ge, dt, n);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_go, d_go, (size_t)n * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_ge, d_ge, (size_t)n * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[5], h_go, (size_t)n * 3);
    write_floats_or_die(argv[6], h_ge, (size_t)n * 3);

    cudaFree(d_gn); cudaFree(d_go); cudaFree(d_ge);
    std::free(h_gn); std::free(h_go); std::free(h_ge);
    return 0;
}

// ── visc_K_partial ────────────────────────────────────────────────────
// Outputs per-particle scalar (n_active floats); the caller (xpbd_full_bwd
// in sub-task F) sums to get ∂L/∂(visc_pair_coef). Mirrors the
// run_spring_K_partial convention — driver writes raw per-particle array,
// not the host-summed scalar. We take `mass` on the CLI because the host
// computes visc_amp from it (same formula as pair_forces_grid_fwd). We do
// NOT take visc_pair_coef — by construction the kernel is computing the
// partial w.r.t. it (the multiplier is factored out of the kernel).
int run_visc_K_partial(int argc, char **argv) {
    if (argc != 18) {
        fprintf(stderr,
            "usage: sib_cuda visc_K_partial "
            "<n_active> <n_static> <h> <mass> <sim_scale> "
            "<pos_active.bin> <vel_active.bin> <sorted_static.bin> "
            "<cell_start.bin> <density.bin> "
            "<grid_dim_x> <grid_dim_y> <grid_dim_z> "
            "<grid_origin.bin> <grad_ext_accel.bin> "
            "<per_particle_out.bin>\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_static = (unsigned int)std::atoi(argv[3]);
    float h         = (float)std::atof(argv[4]);
    float mass      = (float)std::atof(argv[5]);
    float sim_scale = (float)std::atof(argv[6]);
    const char *path_pa  = argv[7];
    const char *path_va  = argv[8];
    const char *path_ss  = argv[9];
    const char *path_cs  = argv[10];
    const char *path_d   = argv[11];
    int grid_dim_x = std::atoi(argv[12]);
    int grid_dim_y = std::atoi(argv[13]);
    int grid_dim_z = std::atoi(argv[14]);
    const char *path_go  = argv[15];
    const char *path_gea = argv[16];
    const char *path_out = argv[17];

    // visc_K_partial only needs visc_amp; surf_amp is computed and discarded.
    float h2, visc_amp, surf_amp;
    compute_pair_amps(h, sim_scale, mass, &h2, &visc_amp, &surf_amp);
    (void)surf_amp;

    int n_cells = grid_dim_x * grid_dim_y * grid_dim_z;
    float *h_pos = read_floats_or_die(path_pa, (size_t)n_active * 3);
    float *h_vel = read_floats_or_die(path_va, (size_t)n_active * 3);
    float *h_ss  = read_floats_or_die(path_ss, (size_t)n_static * 3);
    int   *h_cs  = read_ints_or_die  (path_cs, (size_t)(n_cells + 1));
    float *h_d   = read_floats_or_die(path_d,  (size_t)n_active);
    float *h_go  = read_floats_or_die(path_go, 3);
    float *h_gea = read_floats_or_die(path_gea, (size_t)n_active * 3);

    float3 *d_pos = nullptr, *d_vel = nullptr, *d_ss = nullptr;
    float3 *d_gea = nullptr;
    int *d_cs = nullptr;
    float *d_d = nullptr, *d_pp = nullptr;
    CUDA_CHECK(cudaMalloc(&d_pos, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_vel, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_ss,  (size_t)n_static * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gea, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_cs,  (size_t)(n_cells + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_d,   (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_pp,  (size_t)n_active * sizeof(float)));

    CUDA_CHECK(cudaMemcpy(d_pos, h_pos, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vel, h_vel, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ss,  h_ss,  (size_t)n_static * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gea, h_gea, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cs,  h_cs,  (size_t)(n_cells + 1) * sizeof(int),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_d,   h_d,   (size_t)n_active * sizeof(float),
                          cudaMemcpyHostToDevice));
    // Per-particle output is write-only (kernel does `per_particle[i] = ...`),
    // but zero-init for safety in case n_active==0 / partial launches.
    CUDA_CHECK(cudaMemset(d_pp, 0, (size_t)n_active * sizeof(float)));

    // One thread per particle (Metal: dispatchThreads(n_active, 1, 1) tpg 64).
    const unsigned int TPB = 64;
    unsigned int grid = (n_active + TPB - 1) / TPB;
    visc_K_partial<<<grid, TPB>>>(
        d_pos, d_vel, d_ss, d_cs, d_d, d_gea, d_pp,
        h, h2, sim_scale, visc_amp,
        n_active,
        grid_dim_x, grid_dim_y, grid_dim_z,
        h_go[0], h_go[1], h_go[2]);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float *h_pp = (float *)std::malloc((size_t)n_active * sizeof(float));
    CUDA_CHECK(cudaMemcpy(h_pp, d_pp, (size_t)n_active * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(path_out, h_pp, (size_t)n_active);

    cudaFree(d_pos); cudaFree(d_vel); cudaFree(d_ss);
    cudaFree(d_gea); cudaFree(d_cs);  cudaFree(d_d);
    cudaFree(d_pp);
    std::free(h_pos); std::free(h_vel); std::free(h_ss);
    std::free(h_cs);  std::free(h_d);   std::free(h_go);
    std::free(h_gea); std::free(h_pp);
    return 0;
}
