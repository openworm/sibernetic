// ops_worm.cu — CLI drivers for the worm-physics kernel family.
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
