// metal_common.mm — shared utilities for all sib_metal ops.
//
// Contains:
//   - load_shader_source: locate and cache shaders.metal
//   - MetalCtx + make_ctx + make_pso: per-process device + library + queue
//   - GridDim3 / GridOrigin3 + build_static_grid: spatial-hash for static
//     particles, used by pair_forces_grid

#include "metal_common.h"

#include <libgen.h>      // dirname
#include <unistd.h>      // readlink
#include <mach-o/dyld.h> // _NSGetExecutablePath

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
#include <libgen.h>      // dirname
#include <unistd.h>      // readlink
#include <mach-o/dyld.h> // _NSGetExecutablePath

// Load shader source from src/metal_diff/shaders.metal at runtime.
// Search order:
//   1. Same directory as the running binary (most robust)
//   2. ./shaders.metal (cwd)
//   3. ./src/metal_diff/shaders.metal (when run from repo root)
NSString *load_shader_source(void) {
    static NSString *cached = nil;
    if (cached) return cached;

    // Path 1: alongside the binary
    char binpath[4096];
    uint32_t binpath_len = sizeof(binpath);
    NSMutableArray *candidates = [NSMutableArray array];
    if (_NSGetExecutablePath(binpath, &binpath_len) == 0) {
        char *binpath_dup = strdup(binpath);
        char *bin_dir = dirname(binpath_dup);
        NSString *p1 = [NSString stringWithFormat:@"%s/shaders.metal", bin_dir];
        [candidates addObject:p1];
        free(binpath_dup);
    }
    [candidates addObject:@"shaders.metal"];
    [candidates addObject:@"src/metal_diff/shaders.metal"];

    NSError *err = nil;
    for (NSString *path in candidates) {
        if ([[NSFileManager defaultManager] fileExistsAtPath:path]) {
            cached = [NSString stringWithContentsOfFile:path
                                               encoding:NSUTF8StringEncoding
                                                  error:&err];
            if (cached) return cached;
        }
    }
    fprintf(stderr,
        "shader source not found. Looked in:\n");
    for (NSString *path in candidates)
        fprintf(stderr, "  - %s\n", [path UTF8String]);
    exit(1);
}
// ──────────────────────────────────────────────────────────────────────
// Shared Metal context — created once, reused by all ops in a process.
// Keeps the per-process Metal startup (device init + shader compile,
// ~700ms) amortized across multiple op invocations if a future driver
// calls multiple ops in one binary launch.
// ──────────────────────────────────────────────────────────────────────
MetalCtx make_ctx(void) {
    MetalCtx ctx = {};
    ctx.device = MTLCreateSystemDefaultDevice();
    if (!ctx.device) { fprintf(stderr, "no Metal device\n"); exit(1); }
    NSError *err = nil;
    ctx.lib = [ctx.device
        newLibraryWithSource:load_shader_source()
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

id<MTLComputePipelineState>
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

// PERF — build a uniform spatial hash grid over static particle positions
// (one-time per process invocation). Uses counting sort: O(n_static)
// passes, no comparison-based sort needed.
//
// Outputs:
//   sorted_static[]: permuted positions (cell-grouped contiguous)
//   cell_start[]:    [n_cells + 1] indices; range [cell_start[c],
//                    cell_start[c+1]) is cell c's particles
//   grid_dim, grid_origin, n_cells: returned via out-params
//
// Cell size is the SPH smoothing radius h, so each particle's neighbors
// are guaranteed within the 3×3×3 surrounding cell neighborhood.
void build_static_grid(
    const float *pos_static_in, uint32_t n_static, float h,
    float **sorted_static_out, int **cell_start_out,
    GridDim3 *grid_dim_out, GridOrigin3 *grid_origin_out,
    int *n_cells_out,
    const float *aux_in, float **sorted_aux_out)
{
    if (n_static == 0) {
        *sorted_static_out = NULL;
        *cell_start_out = (int *)calloc(2, sizeof(int));
        grid_dim_out->x = grid_dim_out->y = grid_dim_out->z = 1;
        grid_origin_out->x = grid_origin_out->y = grid_origin_out->z = 0.0f;
        *n_cells_out = 1;
        return;
    }
    // Compute bounding box.
    float bx0 = pos_static_in[0], bx1 = bx0;
    float by0 = pos_static_in[1], by1 = by0;
    float bz0 = pos_static_in[2], bz1 = bz0;
    for (uint32_t i = 1; i < n_static; i++) {
        float x = pos_static_in[i*3+0];
        float y = pos_static_in[i*3+1];
        float z = pos_static_in[i*3+2];
        if (x < bx0) bx0 = x; if (x > bx1) bx1 = x;
        if (y < by0) by0 = y; if (y > by1) by1 = y;
        if (z < bz0) bz0 = z; if (z > bz1) bz1 = z;
    }
    bx0 -= 0.1f; by0 -= 0.1f; bz0 -= 0.1f;
    bx1 += 0.1f; by1 += 0.1f; bz1 += 0.1f;
    int gx = (int)ceilf((bx1 - bx0) / h);
    int gy = (int)ceilf((by1 - by0) / h);
    int gz = (int)ceilf((bz1 - bz0) / h);
    if (gx < 1) gx = 1; if (gy < 1) gy = 1; if (gz < 1) gz = 1;
    int n_cells = gx * gy * gz;
    grid_dim_out->x = gx; grid_dim_out->y = gy; grid_dim_out->z = gz;
    grid_origin_out->x = bx0; grid_origin_out->y = by0; grid_origin_out->z = bz0;
    *n_cells_out = n_cells;

    // Compute cell ID per particle.
    int *cell_ids = (int *)malloc((size_t)n_static * sizeof(int));
    for (uint32_t i = 0; i < n_static; i++) {
        int cx = (int)floorf((pos_static_in[i*3+0] - bx0) / h);
        int cy = (int)floorf((pos_static_in[i*3+1] - by0) / h);
        int cz = (int)floorf((pos_static_in[i*3+2] - bz0) / h);
        if (cx < 0) cx = 0; if (cx >= gx) cx = gx - 1;
        if (cy < 0) cy = 0; if (cy >= gy) cy = gy - 1;
        if (cz < 0) cz = 0; if (cz >= gz) cz = gz - 1;
        cell_ids[i] = cx + cy * gx + cz * gx * gy;
    }
    // Counting sort: cell_start[c] = number of particles with id < c.
    int *cell_start = (int *)calloc((size_t)(n_cells + 1), sizeof(int));
    for (uint32_t i = 0; i < n_static; i++) cell_start[cell_ids[i] + 1]++;
    for (int c = 1; c <= n_cells; c++) cell_start[c] += cell_start[c - 1];
    // Scatter into sorted output. Optionally scatter a parallel aux
    // array (e.g. boundary normals) using the same permutation so it
    // stays index-aligned with sorted_static.
    float *sorted_static = (float *)malloc((size_t)n_static * 3 * sizeof(float));
    float *sorted_aux = NULL;
    if (aux_in != NULL && sorted_aux_out != NULL) {
        sorted_aux = (float *)malloc((size_t)n_static * 3 * sizeof(float));
    }
    int *write_pos = (int *)calloc((size_t)n_cells, sizeof(int));
    for (uint32_t i = 0; i < n_static; i++) {
        int c = cell_ids[i];
        int dst = cell_start[c] + write_pos[c]++;
        sorted_static[dst*3+0] = pos_static_in[i*3+0];
        sorted_static[dst*3+1] = pos_static_in[i*3+1];
        sorted_static[dst*3+2] = pos_static_in[i*3+2];
        if (sorted_aux) {
            sorted_aux[dst*3+0] = aux_in[i*3+0];
            sorted_aux[dst*3+1] = aux_in[i*3+1];
            sorted_aux[dst*3+2] = aux_in[i*3+2];
        }
    }
    free(cell_ids); free(write_pos);
    *sorted_static_out = sorted_static;
    *cell_start_out = cell_start;
    if (sorted_aux_out) *sorted_aux_out = sorted_aux;
}

