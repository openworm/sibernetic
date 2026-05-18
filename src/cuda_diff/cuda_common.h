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

// cuda_common.h — shared includes, CUDA error-check macro, and host I/O
// helpers for the differentiable Sibernetic substrate.
//
// Every .cu in this directory includes this. The .cu files compile
// independently with `nvcc -rdc=true` and device-link at the final
// invocation; kernels defined in shaders.cu are launched from drivers
// living in ops_*.cu via the prototypes in shaders.h.
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstdint>

#include <cuda_runtime.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define CUDA_CHECK(call)                                                   \
    do {                                                                   \
        cudaError_t err__ = (call);                                        \
        if (err__ != cudaSuccess) {                                        \
            fprintf(stderr, "CUDA error %s at %s:%d: %s\n",                \
                    #call, __FILE__, __LINE__, cudaGetErrorString(err__)); \
            std::exit(1);                                                  \
        }                                                                  \
    } while (0)

// Host I/O helpers — fp32 binary slurp/dump used by every CLI driver.
// Defined in cuda_common.cu.
float *read_floats_or_die(const char *path, size_t n_floats);
void   write_floats_or_die(const char *path, const float *buf, size_t n_floats);

// Spatial-grid types used by pair_forces_grid and its consumers.
typedef struct {
    int x, y, z;
} GridDim3;

typedef struct {
    float x, y, z;
} GridOrigin3;

// Build a spatial hash over static (boundary) particles. Returns
// sorted_static (3·n_static floats), cell_start ((n_cells+1) ints),
// grid_dim, grid_origin, n_cells via out-params. Caller frees
// sorted_static and cell_start.
//
// If aux_in is non-NULL (3·n_static floats — e.g. boundary normals),
// it gets sorted with the same permutation as positions and returned
// via sorted_aux_out (caller frees). Pass NULL for both to skip.
void build_static_grid(
    const float *pos_static_in, uint32_t n_static, float h,
    float **sorted_static_out, int **cell_start_out,
    GridDim3 *grid_dim_out, GridOrigin3 *grid_origin_out,
    int *n_cells_out,
    const float *aux_in = NULL, float **sorted_aux_out = NULL);

// ── Kernel block-size contracts ──────────────────────────────────────
// These values are baked into kernel implementations in shaders.cu via
// hard-coded shared-memory sizes and warp-shuffle widths. Changing
// either constant here changes every launch site consistently but does
// NOT change the kernel internals — you'd also need to edit shaders.cu.
// Host drivers reference these symbols at launch so the contract is
// visible at every callsite.
//
//   TPB_REDUCE = 256 : kernels with __shared__ partials[256] +
//                      power-of-2 tree reduce
//                      (rowsum_density, density_constraint_grad,
//                      dist_active_static_bwd, density_constraint_grad_bwd).
//   TPB_PAIR   =  32 : kernels with __shfl_down_sync over one warp
//                      (pair_forces_grid_fwd's 5-round shuffle at
//                      offsets 16/8/4/2/1).
constexpr unsigned int TPB_REDUCE = 256;
constexpr unsigned int TPB_PAIR   = 32;
// Compile-time invariants: these constants are baked into shared-memory
// allocations (`__shared__ float partials[256]`) and warp-shuffle widths
// in shaders.cu. Changing them here without updating the kernel internals
// would silently corrupt results — a static_assert is cheaper than that.
static_assert(TPB_REDUCE == 256,
              "shaders.cu reductions hardcode __shared__ partials[256]");
static_assert(TPB_PAIR == 32,
              "pair_forces_grid_fwd warp-shuffle assumes a full 32-lane warp");
