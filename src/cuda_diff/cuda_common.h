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
