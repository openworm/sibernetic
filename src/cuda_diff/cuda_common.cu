// cuda_common.cu — host I/O helpers used by every CLI driver.
#include "cuda_common.h"

float *read_floats_or_die(const char *path, size_t n_floats) {
    FILE *f = std::fopen(path, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); std::exit(1); }
    float *buf = (float *)std::malloc(n_floats * sizeof(float));
    if (std::fread(buf, sizeof(float), n_floats, f) != n_floats) {
        fprintf(stderr, "short read on %s\n", path); std::exit(1);
    }
    std::fclose(f);
    return buf;
}

void write_floats_or_die(const char *path, const float *buf,
                         size_t n_floats) {
    FILE *f = std::fopen(path, "wb");
    if (!f) { fprintf(stderr, "cannot open %s for write\n", path);
              std::exit(1); }
    std::fwrite(buf, sizeof(float), n_floats, f);
    std::fclose(f);
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
