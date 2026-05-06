// metal_common.h — shared types and helpers for the sib_metal binary.
//
// Implementation lives in metal_common.mm; every ops_*.mm file includes
// this header to pull in MetalCtx, make_ctx, make_pso, the spatial-grid
// types, and the standard Metal/Foundation imports.
#ifndef METAL_COMMON_H
#define METAL_COMMON_H

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Per-process Metal context. Created once via make_ctx; the device,
// shader library, and command queue are reused across all kernel
// invocations in a single binary launch.
typedef struct {
    id<MTLDevice>       device;
    id<MTLLibrary>      lib;
    id<MTLCommandQueue> queue;
} MetalCtx;

// Locate, read, and cache shaders.metal from disk.
NSString *load_shader_source(void);

// Build a MetalCtx (creates device, compiles shaders, opens queue).
MetalCtx make_ctx(void);

// Materialise a compute pipeline state for a named kernel.
id<MTLComputePipelineState> make_pso(MetalCtx ctx, const char *fn_name);

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
void build_static_grid(
    const float *pos_static_in, uint32_t n_static, float h,
    float **sorted_static_out, int **cell_start_out,
    GridDim3 *grid_dim_out, GridOrigin3 *grid_origin_out,
    int *n_cells_out);

#endif // METAL_COMMON_H
