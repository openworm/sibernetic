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
            pair_forces_grid_fwd<<<n_active, 32>>>(
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
        CUDA_CHECK(cudaMemcpyAsync(d_Waa, d_R2aa, r2aa_b,
                                   cudaMemcpyDeviceToDevice));
        wpoly6_inplace<<<(n_aa_total + TPB - 1) / TPB, TPB>>>(
            d_Waa, h2, poly6_const, n_aa_total);
        rowsum_density<<<n_active, 256>>>(d_Waa, d_D_aa, mass,
                                          n_active, n_active);
        if (n_static > 0) {
            {
                dim3 t2(16, 16);
                dim3 b2((n_active + 15) / 16, (n_static + 15) / 16);
                dist_active_static<<<b2, t2>>>(d_Xp, d_Sp, d_R2as,
                                               n_active, n_static);
            }
            unsigned int n_as_total = n_active * n_static;
            CUDA_CHECK(cudaMemcpyAsync(d_Was, d_R2as, r2as_b,
                                       cudaMemcpyDeviceToDevice));
            wpoly6_inplace<<<(n_as_total + TPB - 1) / TPB, TPB>>>(
                d_Was, h2, poly6_const, n_as_total);
            rowsum_density<<<n_active, 256>>>(d_Was, d_D, mass,
                                              n_static, n_active);
            add_inplace<<<gridA, TPB>>>(d_D, d_D_aa, n_active);
        } else {
            // density = density_aa (memcpy).
            CUDA_CHECK(cudaMemcpyAsync(d_D, d_D_aa, s_b,
                                       cudaMemcpyDeviceToDevice));
        }

        // (3) density_constraint_grad.
        density_constraint_grad<<<n_active, 256>>>(
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
