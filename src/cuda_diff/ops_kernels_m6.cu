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

// ops_kernels_m6.cu — M6 op host drivers (wpoly6_inplace, rowsum_density,
// density_constraint_grad, and dist_active_* helpers).
// Kernels in shaders.cu; protos in shaders.h.

#include "cuda_common.h"
#include "shaders.h"
#include "ops.h"


int run_wpoly6_inplace(int argc, char **argv) {
    if (argc != 5 && argc != 6) {
        fprintf(stderr,
            "usage: sib_cuda wpoly6_inplace "
            "<n_total> <h> <inout.bin> [iters]\n");
        return 1;
    }
    unsigned int n_total = (unsigned int)std::atoi(argv[2]);
    float h = (float)std::atof(argv[3]);
    const char *path_inout = argv[4];
    int iters = (argc == 6) ? std::atoi(argv[5]) : 1;
    if (iters < 1) iters = 1;

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));

    FILE *f = std::fopen(path_inout, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path_inout); return 1; }
    float *data = (float *)std::malloc((size_t)n_total * sizeof(float));
    if (std::fread(data, sizeof(float), n_total, f) != n_total) {
        fprintf(stderr, "short read on %s\n", path_inout);
        std::fclose(f);
        std::free(data);
        return 1;
    }
    std::fclose(f);

    float *d_data = nullptr;
    CUDA_CHECK(cudaMalloc(&d_data, (size_t)n_total * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_data, data, (size_t)n_total * sizeof(float),
                          cudaMemcpyHostToDevice));

    const unsigned int TPB = 256;
    unsigned int grid = (n_total + TPB - 1) / TPB;

    // Warm-up + timing strategy matches metal_diff: measure kernel
    // dispatch+sync only. After iter 1 the buffer holds W (mostly zeros)
    // and re-running the kernel gives wrong output but identical compute
    // cost — exactly what we want for steady-state throughput.
    auto t0 = std::chrono::steady_clock::now();
    for (int it = 0; it < iters; it++) {
        wpoly6_inplace<<<grid, TPB>>>(d_data, h2, poly6_const, n_total);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Restore iter-1 result for write-back if iters > 1 (the buffer is
    // currently the iter-N result; rerun once on the input to get iter 1).
    if (iters > 1) {
        CUDA_CHECK(cudaMemcpy(d_data, data, (size_t)n_total * sizeof(float),
                              cudaMemcpyHostToDevice));
        wpoly6_inplace<<<grid, TPB>>>(d_data, h2, poly6_const, n_total);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    CUDA_CHECK(cudaMemcpy(data, d_data, (size_t)n_total * sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_data));

    if (iters > 1) {
        fprintf(stderr, "wpoly6_inplace: %.3f ms/iter (%d iters, %u elems)\n",
                ms / iters, iters, n_total);
    }

    f = std::fopen(path_inout, "wb");
    if (!f) { fprintf(stderr, "cannot open %s for write\n", path_inout);
              std::free(data); return 1; }
    std::fwrite(data, sizeof(float), n_total, f);
    std::fclose(f);
    std::free(data);
    return 0;
}

int run_dist_active_static(int argc, char **argv) {
    if (argc != 7 && argc != 8) {
        fprintf(stderr,
            "usage: sib_cuda dist_active_static "
            "<n_active> <n_static> <active.bin> <static.bin> <out.bin> [iters]\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_static = (unsigned int)std::atoi(argv[3]);
    const char *path_active = argv[4];
    const char *path_static = argv[5];
    const char *path_out    = argv[6];
    int iters = (argc == 8) ? std::atoi(argv[7]) : 1;
    if (iters < 1) iters = 1;

    float *active   = read_floats_or_die(path_active, (size_t)n_active * 3);
    float *static_p = read_floats_or_die(path_static, (size_t)n_static * 3);

    float3 *d_active = nullptr, *d_static = nullptr;
    float *d_dist = nullptr;
    size_t out_bytes = (size_t)n_active * n_static * sizeof(float);
    CUDA_CHECK(cudaMalloc(&d_active, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_static, (size_t)n_static * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_dist, out_bytes));
    CUDA_CHECK(cudaMemcpy(d_active, active,
                          (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_static, static_p,
                          (size_t)n_static * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));

    dim3 threads(16, 16);
    dim3 blocks((n_active + threads.x - 1) / threads.x,
                (n_static + threads.y - 1) / threads.y);

    auto t0 = std::chrono::steady_clock::now();
    for (int it = 0; it < iters; it++) {
        dist_active_static<<<blocks, threads>>>(d_active, d_static, d_dist,
                                                n_active, n_static);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    if (iters > 1) {
        fprintf(stderr, "dist_active_static: %.3f ms/iter (%d iters)\n",
                ms / iters, iters);
    }

    float *out = (float *)std::malloc(out_bytes);
    CUDA_CHECK(cudaMemcpy(out, d_dist, out_bytes, cudaMemcpyDeviceToHost));
    write_floats_or_die(path_out, out, (size_t)n_active * n_static);

    cudaFree(d_active); cudaFree(d_static); cudaFree(d_dist);
    std::free(active); std::free(static_p); std::free(out);
    return 0;
}

int run_dist_active_active(int argc, char **argv) {
    if (argc != 5 && argc != 6) {
        fprintf(stderr,
            "usage: sib_cuda dist_active_active "
            "<n_active> <active.bin> <out.bin> [iters]\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    const char *path_active = argv[3];
    const char *path_out    = argv[4];
    int iters = (argc == 6) ? std::atoi(argv[5]) : 1;
    if (iters < 1) iters = 1;

    float *active = read_floats_or_die(path_active, (size_t)n_active * 3);

    float3 *d_active = nullptr;
    float *d_dist = nullptr;
    size_t out_bytes = (size_t)n_active * n_active * sizeof(float);
    CUDA_CHECK(cudaMalloc(&d_active, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_dist, out_bytes));
    CUDA_CHECK(cudaMemcpy(d_active, active,
                          (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));

    dim3 threads(16, 16);
    dim3 blocks((n_active + threads.x - 1) / threads.x,
                (n_active + threads.y - 1) / threads.y);

    auto t0 = std::chrono::steady_clock::now();
    for (int it = 0; it < iters; it++) {
        dist_active_active<<<blocks, threads>>>(d_active, d_dist, n_active);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    if (iters > 1) {
        fprintf(stderr, "dist_active_active: %.3f ms/iter (%d iters)\n",
                ms / iters, iters);
    }

    float *out = (float *)std::malloc(out_bytes);
    CUDA_CHECK(cudaMemcpy(out, d_dist, out_bytes, cudaMemcpyDeviceToHost));
    write_floats_or_die(path_out, out, (size_t)n_active * n_active);

    cudaFree(d_active); cudaFree(d_dist);
    std::free(active); std::free(out);
    return 0;
}

int run_rowsum_density(int argc, char **argv) {
    if (argc != 7 && argc != 8) {
        fprintf(stderr,
            "usage: sib_cuda rowsum_density "
            "<n_rows> <n_cols> <mass> <W.bin> <density.bin> [iters]\n");
        return 1;
    }
    unsigned int n_rows = (unsigned int)std::atoi(argv[2]);
    unsigned int n_cols = (unsigned int)std::atoi(argv[3]);
    float mass = (float)std::atof(argv[4]);
    const char *path_W   = argv[5];
    const char *path_out = argv[6];
    int iters = (argc == 8) ? std::atoi(argv[7]) : 1;
    if (iters < 1) iters = 1;

    size_t n_W = (size_t)n_rows * n_cols;
    float *W = read_floats_or_die(path_W, n_W);

    float *d_W = nullptr, *d_density = nullptr;
    CUDA_CHECK(cudaMalloc(&d_W, n_W * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_density, (size_t)n_rows * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_W, W, n_W * sizeof(float),
                          cudaMemcpyHostToDevice));

    // rowsum_density declares __shared__ float partials[256] and tree-reduces
    // with stride = T/2; TPB_REDUCE (cuda_common.h) is fixed at 256 to match.
    auto t0 = std::chrono::steady_clock::now();
    for (int it = 0; it < iters; it++) {
        rowsum_density<<<n_rows, TPB_REDUCE>>>(d_W, d_density, mass,
                                               n_cols, n_rows);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    if (iters > 1) {
        fprintf(stderr, "rowsum_density: %.3f ms/iter (%d iters)\n",
                ms / iters, iters);
    }

    float *out = (float *)std::malloc((size_t)n_rows * sizeof(float));
    CUDA_CHECK(cudaMemcpy(out, d_density, (size_t)n_rows * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(path_out, out, (size_t)n_rows);

    cudaFree(d_W); cudaFree(d_density);
    std::free(W); std::free(out);
    return 0;
}

int run_density_constraint_grad(int argc, char **argv) {
    if (argc != 13 && argc != 14) {
        fprintf(stderr,
            "usage: sib_cuda density_constraint_grad "
            "<n_active> <n_static> <h> <mass> <rho_rest> "
            "<active.bin> <static.bin> <r2_aa.bin> <r2_as.bin> "
            "<gradC.bin> <denom_helper.bin> [iters]\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_static = (unsigned int)std::atoi(argv[3]);
    float h        = (float)std::atof(argv[4]);
    float mass     = (float)std::atof(argv[5]);
    float rho_rest = (float)std::atof(argv[6]);
    const char *path_active = argv[7];
    const char *path_static = argv[8];
    const char *path_r2_aa  = argv[9];
    const char *path_r2_as  = argv[10];
    const char *path_gradC  = argv[11];
    const char *path_denomH = argv[12];
    int iters = (argc == 14) ? std::atoi(argv[13]) : 1;
    if (iters < 1) iters = 1;

    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));

    float *active   = read_floats_or_die(path_active, (size_t)n_active * 3);
    float *static_p = read_floats_or_die(path_static, (size_t)n_static * 3);
    float *r2_aa    = read_floats_or_die(path_r2_aa,
                                         (size_t)n_active * n_active);
    float *r2_as    = read_floats_or_die(path_r2_as,
                                         (size_t)n_active * n_static);

    float3 *d_active = nullptr, *d_static = nullptr, *d_grad = nullptr;
    float *d_r2aa = nullptr, *d_r2as = nullptr, *d_denom = nullptr;
    CUDA_CHECK(cudaMalloc(&d_active, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_static, (size_t)n_static * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_r2aa,
                          (size_t)n_active * n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_r2as,
                          (size_t)n_active * n_static * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_grad, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_denom, (size_t)n_active * sizeof(float)));

    CUDA_CHECK(cudaMemcpy(d_active, active,
                          (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_static, static_p,
                          (size_t)n_static * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_r2aa, r2_aa,
                          (size_t)n_active * n_active * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_r2as, r2_as,
                          (size_t)n_active * n_static * sizeof(float),
                          cudaMemcpyHostToDevice));

    // density_constraint_grad declares __shared__ float3 partials_grad[256]
    // + __shared__ float partials_denom[256] and tree-reduces with stride
    // = T/2; TPB_REDUCE (cuda_common.h) is fixed at 256 to match.
    auto t0 = std::chrono::steady_clock::now();
    for (int it = 0; it < iters; it++) {
        density_constraint_grad<<<n_active, TPB_REDUCE>>>(
            d_active, d_static, d_r2aa, d_r2as, d_grad, d_denom,
            h, spiky_const, mass, rho_rest, n_active, n_static);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    if (iters > 1) {
        fprintf(stderr, "density_constraint_grad: %.3f ms/iter (%d iters)\n",
                ms / iters, iters);
    }

    float *grad = (float *)std::malloc((size_t)n_active * 3 * sizeof(float));
    float *denom = (float *)std::malloc((size_t)n_active * sizeof(float));
    CUDA_CHECK(cudaMemcpy(grad, d_grad, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(denom, d_denom, (size_t)n_active * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(path_gradC, grad, (size_t)n_active * 3);
    write_floats_or_die(path_denomH, denom, (size_t)n_active);

    cudaFree(d_active); cudaFree(d_static); cudaFree(d_r2aa);
    cudaFree(d_r2as); cudaFree(d_grad); cudaFree(d_denom);
    std::free(active); std::free(static_p); std::free(r2_aa); std::free(r2_as);
    std::free(grad); std::free(denom);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M7 — XPBD step kernels
// ──────────────────────────────────────────────────────────────────────


int run_dist_active_static_bwd(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr,
            "usage: sib_cuda dist_active_static_bwd "
            "<n_active> <n_static> <active.bin> <static.bin> <grad_r2.bin> "
            "<grad_active_inout.bin>\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_static = (unsigned int)std::atoi(argv[3]);
    float *h_active = read_floats_or_die(argv[4], (size_t)n_active * 3);
    float *h_static = read_floats_or_die(argv[5], (size_t)n_static * 3);
    float *h_gr2    = read_floats_or_die(argv[6], (size_t)n_active * n_static);
    float *h_gact   = read_floats_or_die(argv[7], (size_t)n_active * 3);

    float3 *d_active = nullptr, *d_static = nullptr, *d_gact = nullptr;
    float *d_gr2 = nullptr;
    CUDA_CHECK(cudaMalloc(&d_active, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_static, (size_t)n_static * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gr2, (size_t)n_active * n_static * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_gact, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMemcpy(d_active, h_active, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_static, h_static, (size_t)n_static * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gr2, h_gr2,
                          (size_t)n_active * n_static * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gact, h_gact, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));

    // dist_active_static_bwd uses __shared__ float3 partials[256] + tree
    // reduction with stride = T/2; TPB_REDUCE (cuda_common.h) is fixed at
    // 256 to match the shared-mem size.
    dist_active_static_bwd<<<n_active, TPB_REDUCE>>>(
        d_active, d_static, d_gr2, d_gact, n_active, n_static);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_gact, d_gact, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[7], h_gact, (size_t)n_active * 3);

    cudaFree(d_active); cudaFree(d_static); cudaFree(d_gr2); cudaFree(d_gact);
    std::free(h_active); std::free(h_static); std::free(h_gr2); std::free(h_gact);
    return 0;
}

// M6.3_bwd — backward of dist_active_active.
// Both endpoints are active. Per i:

int run_dist_active_active_bwd(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr,
            "usage: sib_cuda dist_active_active_bwd "
            "<n_active> <active.bin> <grad_r2.bin> <grad_active_inout.bin>\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    float *h_active = read_floats_or_die(argv[3], (size_t)n_active * 3);
    float *h_gr2    = read_floats_or_die(argv[4], (size_t)n_active * n_active);
    float *h_gact   = read_floats_or_die(argv[5], (size_t)n_active * 3);

    float3 *d_active = nullptr, *d_gact = nullptr;
    float *d_gr2 = nullptr;
    CUDA_CHECK(cudaMalloc(&d_active, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gr2, (size_t)n_active * n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_gact, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMemcpy(d_active, h_active, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gr2, h_gr2,
                          (size_t)n_active * n_active * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gact, h_gact, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyHostToDevice));

    const unsigned int TPB = 128;
    unsigned int grid = (n_active + TPB - 1) / TPB;
    dist_active_active_bwd<<<grid, TPB>>>(d_active, d_gr2, d_gact, n_active);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_gact, d_gact, (size_t)n_active * 3 * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[5], h_gact, (size_t)n_active * 3);

    cudaFree(d_active); cudaFree(d_gr2); cudaFree(d_gact);
    std::free(h_active); std::free(h_gr2); std::free(h_gact);
    return 0;
}

// M6.1_bwd — in-place backward of wpoly6_inplace.
// Forward: W = poly6_const · (h² - r²)³ if r² < h² else 0.
//   dW/dr² = -3 · poly6_const · (h² - r²)²
// In-place: input holds ∂L/∂W, output (same buffer) is ∂L/∂r².
int run_wpoly6_inplace_bwd(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr,
            "usage: sib_cuda wpoly6_inplace_bwd "
            "<n_total> <h> <r2_saved.bin> <grad_W_or_r2_inout.bin>\n");
        return 1;
    }
    unsigned int n_total = (unsigned int)std::atoi(argv[2]);
    float h = (float)std::atof(argv[3]);
    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));

    float *h_r2  = read_floats_or_die(argv[4], (size_t)n_total);
    float *h_buf = read_floats_or_die(argv[5], (size_t)n_total);

    float *d_r2 = nullptr, *d_buf = nullptr;
    CUDA_CHECK(cudaMalloc(&d_r2, (size_t)n_total * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_buf, (size_t)n_total * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_r2, h_r2, (size_t)n_total * sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_buf, h_buf, (size_t)n_total * sizeof(float),
                          cudaMemcpyHostToDevice));

    const unsigned int TPB = 256;
    unsigned int grid = (n_total + TPB - 1) / TPB;
    wpoly6_inplace_bwd<<<grid, TPB>>>(d_r2, d_buf, h2, poly6_const, n_total);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_buf, d_buf, (size_t)n_total * sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[5], h_buf, (size_t)n_total);

    cudaFree(d_r2); cudaFree(d_buf);
    std::free(h_r2); std::free(h_buf);
    return 0;
}

// M6.2_bwd — backward of rowsum_density (broadcast).
// Forward: density[i] = mass · Σ_j W[i,j]
//   ∂L/∂W[i,j] = mass · ∂L/∂density[i]   (constant across j)

int run_rowsum_density_bwd(int argc, char **argv) {
    if (argc != 7) {
        fprintf(stderr,
            "usage: sib_cuda rowsum_density_bwd "
            "<n_rows> <n_cols> <mass> <grad_density.bin> <grad_W_out.bin>\n");
        return 1;
    }
    unsigned int n_rows = (unsigned int)std::atoi(argv[2]);
    unsigned int n_cols = (unsigned int)std::atoi(argv[3]);
    float mass = (float)std::atof(argv[4]);
    float *h_gd = read_floats_or_die(argv[5], (size_t)n_rows);

    float *d_gd = nullptr, *d_gW = nullptr;
    size_t W_bytes = (size_t)n_rows * n_cols * sizeof(float);
    CUDA_CHECK(cudaMalloc(&d_gd, (size_t)n_rows * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_gW, W_bytes));
    CUDA_CHECK(cudaMemcpy(d_gd, h_gd, (size_t)n_rows * sizeof(float),
                          cudaMemcpyHostToDevice));

    dim3 threads(16, 16);
    dim3 blocks((n_cols + threads.x - 1) / threads.x,
                (n_rows + threads.y - 1) / threads.y);
    rowsum_density_bwd<<<blocks, threads>>>(d_gd, d_gW, mass, n_rows, n_cols);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float *h_gW = (float *)std::malloc(W_bytes);
    CUDA_CHECK(cudaMemcpy(h_gW, d_gW, W_bytes, cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[6], h_gW, (size_t)n_rows * n_cols);

    cudaFree(d_gd); cudaFree(d_gW);
    std::free(h_gd); std::free(h_gW);
    return 0;
}

// M6.4_bwd — backward of density_constraint_grad.
//
// Forward (per active i):
//   grad_C[i]       = (m/ρ_rest) · Σ_{k≠i} ∇W_spiky(p_i - p_k)
//   denom_helper[i] = Σ_{k≠i} |∇W_spiky(p_i - p_k)|²
// where k iterates over both active (skip self) and static neighbors.
//
// ∇W_spiky(v) = c·(h-r)²·v/r with c=spiky_const, r=|v|.
// Jacobian:
//   J = c·(h-r) · [(h-r)/r · (I - ĝĝᵀ) - 2·ĝĝᵀ]
// For a row vector u:
//   u·J = c·(h-r) · [(h-r)/r · (u - (u·ĝ)·ĝ) - 2·(u·ĝ)·ĝ]
//
// Per pair (i, k), let u_self = (m/ρ)·ω_i + 2·ψ_i·grad_W_ik.
// Contribution to ∂L/∂p_i from row i: +u_self·J.
// If k is active (k=j), the SAME pair (i,k) also appears in row j as
// pair (j,i): grad_W_ji = -grad_W_ij ⟹ u_neigh = (m/ρ)·ω_j - 2·ψ_j·grad_W_ij.
// Then ∂grad_W_ji/∂p_i = -J, so contribution to ∂L/∂p_i: -u_neigh·J.
// Static neighbors contribute only via row i (static positions are frozen).
//

int run_density_constraint_grad_bwd(int argc, char **argv) {
    if (argc != 14) {
        fprintf(stderr,
            "usage: sib_cuda density_constraint_grad_bwd "
            "<n_active> <n_static> <h> <mass> <rho_rest> "
            "<active.bin> <static.bin> <r2_aa.bin> <r2_as.bin> "
            "<grad_grad_C.bin> <grad_denom_h.bin> <grad_active_inout.bin>\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    unsigned int n_static = (unsigned int)std::atoi(argv[3]);
    float h        = (float)std::atof(argv[4]);
    float mass     = (float)std::atof(argv[5]);
    float rho_rest = (float)std::atof(argv[6]);
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));

    float *h_active = read_floats_or_die(argv[7], (size_t)n_active * 3);
    float *h_static = read_floats_or_die(argv[8], (size_t)n_static * 3);
    float *h_r2aa   = read_floats_or_die(argv[9], (size_t)n_active * n_active);
    float *h_r2as   = read_floats_or_die(argv[10], (size_t)n_active * n_static);
    float *h_gg     = read_floats_or_die(argv[11], (size_t)n_active * 3);
    float *h_gph    = read_floats_or_die(argv[12], (size_t)n_active);
    float *h_gact   = read_floats_or_die(argv[13], (size_t)n_active * 3);

    float3 *d_active = nullptr, *d_static = nullptr, *d_gg = nullptr, *d_gact = nullptr;
    float *d_r2aa = nullptr, *d_r2as = nullptr, *d_gph = nullptr;
    CUDA_CHECK(cudaMalloc(&d_active, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_static, (size_t)n_static * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_r2aa, (size_t)n_active * n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_r2as, (size_t)n_active * n_static * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_gg,  (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gph, (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_gact, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMemcpy(d_active, h_active, (size_t)n_active*3*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_static, h_static, (size_t)n_static*3*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_r2aa, h_r2aa,
                          (size_t)n_active*n_active*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_r2as, h_r2as,
                          (size_t)n_active*n_static*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gg,  h_gg, (size_t)n_active*3*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gph, h_gph, (size_t)n_active*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gact, h_gact, (size_t)n_active*3*sizeof(float),
                          cudaMemcpyHostToDevice));

    const unsigned int TPB = 128;
    unsigned int grid = (n_active + TPB - 1) / TPB;
    density_constraint_grad_bwd<<<grid, TPB>>>(
        d_active, d_static, d_r2aa, d_r2as, d_gg, d_gph, d_gact,
        h, spiky_const, mass, rho_rest, n_active, n_static);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_gact, d_gact, (size_t)n_active*3*sizeof(float),
                          cudaMemcpyDeviceToHost));
    write_floats_or_die(argv[13], h_gact, (size_t)n_active * 3);

    cudaFree(d_active); cudaFree(d_static);
    cudaFree(d_r2aa); cudaFree(d_r2as);
    cudaFree(d_gg); cudaFree(d_gph); cudaFree(d_gact);
    std::free(h_active); std::free(h_static);
    std::free(h_r2aa); std::free(h_r2as);
    std::free(h_gg); std::free(h_gph); std::free(h_gact);
    return 0;
}

// M7.1_bwd — backward of solve_density_constraint.
//
// Forward (per active i, if C > 0):
//   C = density/ρ - 1
//   D = |g|²/m + (m/ρ²)·helper + A           (A = α/dt²)
//   Δλ = -(C + A·λ_pre) / D
//   pos_post = pos_pre + g·Δλ/m
//   λ_post  = λ_pre + Δλ
//
// Define chain = (ω·g)/m + ψ. Outputs:
//   ∂L/∂pos_pre   += ω                                  (identity)
//   ∂L/∂λ_pre     = -ω_dot_g · A/(m·D) + ψ·(1 - A/D)
//   ∂L/∂density   = -chain / (D · ρ)
//   ∂L/∂g         = (Δλ/m)·ω - (2Δλ/(m·D))·chain·g
//   ∂L/∂helper    = -Δλ · m / (ρ²·D) · chain
//   ∂L/∂ρ          = chain · [density/(ρ²·D) + 2·Δλ·m·helper/(ρ³·D)]   (per-i)
//   ∂L/∂A          = chain · (-λ_post / D)                              (per-i)
// (∂L/∂ρ and ∂L/∂A are per-particle; host sums to scalar.)
//
// If C ≤ 0: forward skipped (one-sided projection), backward is identity

int run_solve_density_constraint_bwd(int argc, char **argv) {
    if (argc != 19) {
        fprintf(stderr,
            "usage: sib_cuda solve_density_constraint_bwd "
            "<n_active> <rho_rest> <mass> <alpha_inv_dt2> "
            "<density.bin> <grad_C.bin> <denom_helper.bin> <lambda_pre.bin> "
            "<grad_pos_post.bin> <grad_lambda_post.bin> "
            "<grad_pos_pre_inout.bin> <grad_lambda_pre.bin> <grad_density.bin> "
            "<grad_grad_C.bin> <grad_denom_h.bin> <grad_rho_rest.bin> <grad_alpha.bin>\n");
        return 1;
    }
    unsigned int n_active = (unsigned int)std::atoi(argv[2]);
    float rho_rest      = (float)std::atof(argv[3]);
    float mass          = (float)std::atof(argv[4]);
    float alpha_inv_dt2 = (float)std::atof(argv[5]);

    float *h_density = read_floats_or_die(argv[6], (size_t)n_active);
    float *h_gC      = read_floats_or_die(argv[7], (size_t)n_active * 3);
    float *h_dh      = read_floats_or_die(argv[8], (size_t)n_active);
    float *h_lp      = read_floats_or_die(argv[9], (size_t)n_active);
    float *h_gpp     = read_floats_or_die(argv[10], (size_t)n_active * 3);
    float *h_glp     = read_floats_or_die(argv[11], (size_t)n_active);
    float *h_gppre   = read_floats_or_die(argv[12], (size_t)n_active * 3);

    float *d_density = nullptr, *d_dh = nullptr, *d_lp = nullptr;
    float *d_glp = nullptr;
    float3 *d_gC = nullptr, *d_gpp = nullptr, *d_gppre = nullptr;
    float *d_glpre = nullptr, *d_gdens = nullptr, *d_gdh = nullptr;
    float3 *d_ggC = nullptr;
    float *d_grho = nullptr, *d_galpha = nullptr;
    CUDA_CHECK(cudaMalloc(&d_density, (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_gC,  (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_dh,  (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_lp,  (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_gpp, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_glp, (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_gppre, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_glpre, (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_gdens, (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_ggC, (size_t)n_active * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_gdh, (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_grho, (size_t)n_active * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_galpha, (size_t)n_active * sizeof(float)));

    CUDA_CHECK(cudaMemcpy(d_density, h_density, (size_t)n_active*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gC, h_gC, (size_t)n_active*3*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_dh, h_dh, (size_t)n_active*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_lp, h_lp, (size_t)n_active*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gpp, h_gpp, (size_t)n_active*3*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_glp, h_glp, (size_t)n_active*sizeof(float),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gppre, h_gppre, (size_t)n_active*3*sizeof(float),
                          cudaMemcpyHostToDevice));

    const unsigned int TPB = 128;
    unsigned int grid = (n_active + TPB - 1) / TPB;
    solve_density_constraint_bwd<<<grid, TPB>>>(
        d_density, d_gC, d_dh, d_lp, d_gpp, d_glp,
        d_gppre, d_glpre, d_gdens, d_ggC, d_gdh, d_grho, d_galpha,
        rho_rest, mass, alpha_inv_dt2, n_active);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_gppre, d_gppre, (size_t)n_active*3*sizeof(float),
                          cudaMemcpyDeviceToHost));
    float *h_glpre  = (float*)std::malloc((size_t)n_active*sizeof(float));
    float *h_gdens  = (float*)std::malloc((size_t)n_active*sizeof(float));
    float *h_ggC    = (float*)std::malloc((size_t)n_active*3*sizeof(float));
    float *h_gdh    = (float*)std::malloc((size_t)n_active*sizeof(float));
    float *h_grho   = (float*)std::malloc((size_t)n_active*sizeof(float));
    float *h_galpha = (float*)std::malloc((size_t)n_active*sizeof(float));
    CUDA_CHECK(cudaMemcpy(h_glpre,  d_glpre, (size_t)n_active*sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_gdens,  d_gdens, (size_t)n_active*sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_ggC,    d_ggC, (size_t)n_active*3*sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_gdh,    d_gdh, (size_t)n_active*sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_grho,   d_grho, (size_t)n_active*sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_galpha, d_galpha, (size_t)n_active*sizeof(float),
                          cudaMemcpyDeviceToHost));

    write_floats_or_die(argv[12], h_gppre,  (size_t)n_active * 3);
    write_floats_or_die(argv[13], h_glpre,  (size_t)n_active);
    write_floats_or_die(argv[14], h_gdens,  (size_t)n_active);
    write_floats_or_die(argv[15], h_ggC,    (size_t)n_active * 3);
    write_floats_or_die(argv[16], h_gdh,    (size_t)n_active);
    write_floats_or_die(argv[17], h_grho,   (size_t)n_active);
    write_floats_or_die(argv[18], h_galpha, (size_t)n_active);

    cudaFree(d_density); cudaFree(d_gC); cudaFree(d_dh); cudaFree(d_lp);
    cudaFree(d_gpp); cudaFree(d_glp); cudaFree(d_gppre); cudaFree(d_glpre);
    cudaFree(d_gdens); cudaFree(d_ggC); cudaFree(d_gdh);
    cudaFree(d_grho); cudaFree(d_galpha);
    std::free(h_density); std::free(h_gC); std::free(h_dh); std::free(h_lp);
    std::free(h_gpp); std::free(h_glp); std::free(h_gppre); std::free(h_glpre);
    std::free(h_gdens); std::free(h_ggC); std::free(h_gdh);
    std::free(h_grho); std::free(h_galpha);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// xpbd_full_fwd_step / xpbd_full_bwd_step
//
// Single-step differentiable XPBD subset (predict + floor + update only).
// Used as the Phase 4 differentiability proof. Density + bonds backwards
// land in a follow-on once their forward analogs need SGD.
//
// fwd CLI:
//   sib_cuda xpbd_full_fwd_step
//       <n> <dt> <gravity_y> <floor_y> <restitution> <sim_scale>
//       <pos_old.bin> <vel.bin>
//       <pos_pred_save.bin> <pos_post_save.bin> <vel_new.bin>
//
// bwd CLI:
//   sib_cuda xpbd_full_bwd_step
//       <n> <dt> <gravity_y> <floor_y> <restitution> <sim_scale>
//       <pos_pred_save.bin> <pos_old.bin>
//       <grad_vel_new.bin>
//       <grad_pos_old.bin> <grad_vel.bin>
//       <grad_gravity_y.bin> <grad_floor_y.bin> <grad_restitution.bin>
