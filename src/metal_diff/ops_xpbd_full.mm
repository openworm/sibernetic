// ops_xpbd_full.mm — differentiable XPBD pipeline.
//
// xpbd_full_fwd: runs N steps, saves all per-step intermediate state
//                needed for backward.
// xpbd_full_bwd: walks the chain in reverse using the saved state, runs
//                each per-step backward kernel, and produces analytic
//                gradients of L wrt initial state and parameters
//                (rho_rest, spring_K, visc_pair_coef, alpha_dens).

#include "metal_common.h"

// ──────────────────────────────────────────────────────────────────────
// Multi-step XPBD forward+backward (xpbd_full_fwd / _bwd).
//
// Pipeline per step k (n_iters=1 for now to keep state simple):
//   predict_positions
//   density chain: dist_aa + dist_as + wpoly6×2 + rowsum×2 + add
//   density_constraint_grad → grad_C, denom_helper
//   solve_density_constraint
//   update_velocities
//
// Saved state per step (all on host, per-step indexing):
//   x_old[k], v_old[k]          ← 6 floats per particle per step
//   density[k]                   ← 1 float per particle per step
//   grad_C[k]                    ← 3 floats per particle per step
//   denom_helper[k]              ← 1 float per particle per step
//
// r²_aa, r²_as are RECOMPUTED in backward (cheap distance ops) to
// avoid memory blowup at demo1 scale (r²_as = 24 MB per step).
//
// No bonds, no floor in this version — to be added as separate
// extensions once base multi-step backward validates.
// ──────────────────────────────────────────────────────────────────────
int run_xpbd_full_fwd(int argc, char **argv) {
    // Optional 15th arg: sim_scale (default 1.0 — toy convention)
    // Optional 16th arg: visc_pair_coef (default 0 — pair forces off)
    if (argc < 15 || argc > 26) {
        fprintf(stderr,
            "usage: sib_metal xpbd_full_fwd "
            "<n_active> <n_static> <K> <h> <mass> <rho_rest> <dt> "
            "<gravity_y> <alpha_density> "
            "<pos0.bin> <vel0.bin> <pos_static.bin> <state_out.bin> "
            "[sim_scale] [visc_pair_coef] [spring_K] [bonds.bin] [floor_y] "
            "[restitution] [n_membranes] [n_elastic] [r0] "
            "[membranes.bin] [pmem_index.bin]\n"
            "  state_out.bin contains the per-step trajectory: \n"
            "    K × [x_old(n*3) + v_old(n*3) + density(n) + grad_C(n*3) + denom_h(n)"
            "    [+ pair_density(n) if visc_pair_coef>0]"
            "    [+ floor_clamped_mask(n) as int32 if floor_y is set]]\n"
            "    Plus K+1 frames of x_post for the trajectory.\n"
            "    Plus final v(n*3).\n"
            "  When spring_K > 0, bonds.bin must be provided (16 bytes per\n"
            "  bond: int32 i, int32 j, float32 rest_len, float32 _pad).\n"
            "  When floor_y is set, elastic floor at y=floor_y is applied\n"
            "  after density solve. restitution ∈ [0, 1] (default 0):\n"
            "    0 = inelastic clamp (legacy), 1 = perfectly elastic bounce.\n"
            "  n_membranes > 0 enables M10 membrane interaction (forward only;\n"
            "    state-file format is unchanged in Phase 3 — xpbd_full_bwd\n"
            "    will not yet attribute gradients through membranes).\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    uint32_t K        = (uint32_t)atoi(argv[4]);
    float h           = (float)atof(argv[5]);
    float mass        = (float)atof(argv[6]);
    float rho_rest    = (float)atof(argv[7]);
    float dt          = (float)atof(argv[8]);
    float g_y         = (float)atof(argv[9]);
    float alpha_dens  = (float)atof(argv[10]);
    float sim_scale   = (argc >= 16) ? (float)atof(argv[15]) : 1.0f;
    float sim_scale_inv = 1.0f / sim_scale;
    float visc_pair_coef = (argc >= 17) ? (float)atof(argv[16]) : 0.0f;
    bool use_pair = visc_pair_coef != 0.0f;
    float spring_K = (argc >= 18) ? (float)atof(argv[17]) : 0.0f;
    const char *bonds_path = (argc >= 19) ? argv[18] : NULL;
    bool use_springs = spring_K != 0.0f;
    bool need_ext_accel = use_pair || use_springs;
    if (use_springs && !bonds_path) {
        fprintf(stderr, "spring_K > 0 requires bonds.bin path as 19th arg\n");
        return 1;
    }
    // Optional 20th arg: floor_y. Default = NaN means "no floor" (legacy).
    bool use_floor = (argc >= 20);
    float floor_y = use_floor ? (float)atof(argv[19]) : 0.0f;
    // Optional 21st arg: restitution coefficient (default 0 = inelastic).
    float restitution = (argc >= 21) ? (float)atof(argv[20]) : 0.0f;
    // M10 membrane interaction (forward only in Phase 3).
    uint32_t n_membranes = (argc >= 22) ? (uint32_t)atoi(argv[21]) : 0u;
    uint32_t n_elastic   = (argc >= 23) ? (uint32_t)atoi(argv[22]) : 0u;
    float r0_mem         = (argc >= 24) ? (float)atof(argv[23]) : 0.0f;
    const char *path_membranes = (argc >= 25) ? argv[24] : NULL;
    const char *path_pmem_idx  = (argc >= 26) ? argv[25] : NULL;
    bool use_membranes = (n_membranes > 0) && (n_elastic > 0)
                         && path_membranes && path_pmem_idx
                         && (r0_mem > 0.0f);
    // Load bonds if springs enabled.
    uint32_t n_bonds = 0;
    int32_t *bond_ij_data = NULL;
    float *bond_rest_data = NULL;
    if (use_springs) {
        FILE *fb = fopen(bonds_path, "rb");
        if (!fb) { fprintf(stderr, "open %s\n", bonds_path); return 1; }
        fseek(fb, 0, SEEK_END);
        long sz = ftell(fb);
        fseek(fb, 0, SEEK_SET);
        n_bonds = (uint32_t)(sz / 16);
        uint8_t *raw = (uint8_t *)malloc((size_t)n_bonds * 16);
        fread(raw, 1, (size_t)n_bonds * 16, fb);
        fclose(fb);
        bond_ij_data = (int32_t *)malloc((size_t)n_bonds * 2 * sizeof(int32_t));
        bond_rest_data = (float *)malloc((size_t)n_bonds * sizeof(float));
        for (uint32_t b = 0; b < n_bonds; b++) {
            memcpy(&bond_ij_data[b * 2], raw + b * 16, 8);
            memcpy(&bond_rest_data[b], raw + b * 16 + 8, 4);
        }
        free(raw);
    }

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));
    float alpha_inv_dt2 = alpha_dens / (dt * dt);
    // Precomputed pair-force amps (in fp64 to dodge fp32 underflow when
    // sim_scale is tiny; cast to fp32 at kernel boundary).
    double h_scaled  = (double)h * (double)sim_scale;
    double h_s6      = pow(h_scaled, 6.0);
    double h_s9      = pow(h_scaled, 9.0);
    double divgradWvisco = 45.0 / (M_PI * h_s6);
    float visc_amp = (float)(1.5 * (double)mass * divgradWvisco
                              * pow((double)sim_scale, 3.0));
    double wpoly6_si = 315.0 / (64.0 * M_PI * h_s9);
    float surf_amp = (float)(-1.7e-9 * (double)mass * wpoly6_si
                              * (double)sim_scale / (double)mass);

    auto rd = ^(const char *p, size_t n_floats) {
        FILE *f = fopen(p, "rb");
        if (!f) { fprintf(stderr, "open %s\n", p); exit(1); }
        float *b = (float *)malloc(n_floats * sizeof(float));
        fread(b, sizeof(float), n_floats, f); fclose(f); return b;
    };
    float *pos0 = rd(argv[11], (size_t)n_active * 3);
    float *vel0 = rd(argv[12], (size_t)n_active * 3);
    float *pos_static = rd(argv[13], (size_t)n_static * 3);

    // Load membrane topology (int32 buffers, sized in bytes — read via fopen).
    int32_t *mem_tris = NULL;
    int32_t *mem_pidx = NULL;
    if (use_membranes) {
        size_t tris_bytes = (size_t)n_membranes * 3 * sizeof(int32_t);
        size_t pidx_bytes = (size_t)n_elastic * 7 * sizeof(int32_t);
        FILE *fa = fopen(path_membranes, "rb");
        if (!fa) { fprintf(stderr, "open %s\n", path_membranes); return 1; }
        mem_tris = (int32_t *)malloc(tris_bytes);
        if (fread(mem_tris, 1, tris_bytes, fa) != tris_bytes) {
            fprintf(stderr, "short read on %s\n", path_membranes); return 1;
        }
        fclose(fa);
        FILE *fb2 = fopen(path_pmem_idx, "rb");
        if (!fb2) { fprintf(stderr, "open %s\n", path_pmem_idx); return 1; }
        mem_pidx = (int32_t *)malloc(pidx_bytes);
        if (fread(mem_pidx, 1, pidx_bytes, fb2) != pidx_bytes) {
            fprintf(stderr, "short read on %s\n", path_pmem_idx); return 1;
        }
        fclose(fb2);
    }

    size_t pos_b = (size_t)n_active * 3 * sizeof(float);
    size_t s_b = (size_t)n_active * sizeof(float);

    // Per-step state arrays (host).
    // Layout per step: [x_old(3n) + v_old(3n) + density(n) + grad_C(3n) + denom_h(n)
    //                   + pair_density(n) if use_pair
    //                   + floor_clamped_mask(n as int32) if use_floor
    //                   + mem_corr(3n) if use_membranes]
    int extra_per_step = (use_pair ? 1 : 0) + (use_floor ? 1 : 0)
                       + (use_membranes ? 3 : 0);
    size_t per_step_floats = (size_t)n_active * (3 + 3 + 1 + 3 + 1 + extra_per_step);
    float *state = (float *)calloc((size_t)K * per_step_floats, sizeof(float));
    // Static spatial grid buffers — declared in outer scope so the
    // free() at function end can see them.
    float *sorted_static_buf = NULL;
    int *cell_start_buf = NULL;
    // Plus K+1 frames of x_post for the trajectory.
    float *traj = (float *)calloc((size_t)(K + 1) * (size_t)n_active * 3, sizeof(float));
    memcpy(traj, pos0, pos_b);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pred  = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_d_aa  = make_pso(ctx, "dist_active_active");
        id<MTLComputePipelineState> pso_d_as  = make_pso(ctx, "dist_active_static");
        id<MTLComputePipelineState> pso_wp    = make_pso(ctx, "wpoly6_inplace");
        id<MTLComputePipelineState> pso_rs    = make_pso(ctx, "rowsum_density");
        id<MTLComputePipelineState> pso_addin = make_pso(ctx, "add_inplace");
        id<MTLComputePipelineState> pso_dgrad = make_pso(ctx, "density_constraint_grad");
        id<MTLComputePipelineState> pso_solv  = make_pso(ctx, "solve_density_constraint");
        id<MTLComputePipelineState> pso_uv    = make_pso(ctx, "update_velocities");
        id<MTLComputePipelineState> pso_pair  = use_pair
            ? make_pso(ctx, "pair_forces_grid") : nil;
        id<MTLComputePipelineState> pso_appext = need_ext_accel
            ? make_pso(ctx, "apply_ext_accel") : nil;
        id<MTLComputePipelineState> pso_spring = use_springs
            ? make_pso(ctx, "spring_bonds_force") : nil;
        id<MTLComputePipelineState> pso_floor = use_floor
            ? make_pso(ctx, "solve_floor_constraint_with_mask") : nil;
        // M10 membrane PSOs.
        id<MTLComputePipelineState> pso_mem_clear = use_membranes
            ? make_pso(ctx, "clear_membrane_correction") : nil;
        id<MTLComputePipelineState> pso_mem_acc   = use_membranes
            ? make_pso(ctx, "accumulate_membrane_correction") : nil;
        id<MTLComputePipelineState> pso_mem_apply = use_membranes
            ? make_pso(ctx, "apply_membrane_correction") : nil;

        size_t r2aa_b = (size_t)n_active * n_active * sizeof(float);
        size_t r2as_b = (size_t)n_active * n_static * sizeof(float);

        // Build static spatial grid if pair_forces is on (one-time cost
        // per process). Identical to xpbd_step's grid setup.
        GridDim3 grid_dim_struct = {0, 0, 0};
        GridOrigin3 grid_origin_struct = {0, 0, 0};
        int n_cells = 0;
        if (use_pair && n_static > 0) {
            build_static_grid(pos_static, n_static, h,
                              &sorted_static_buf, &cell_start_buf,
                              &grid_dim_struct, &grid_origin_struct, &n_cells);
        }

        id<MTLBuffer> bX  = [ctx.device newBufferWithBytes:pos0 length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bV  = [ctx.device newBufferWithBytes:vel0 length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bXp = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSp = [ctx.device newBufferWithBytes:pos_static
            length:(size_t)n_static * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2aa = [ctx.device newBufferWithLength:r2aa_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2as = [ctx.device newBufferWithLength:r2as_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bWaa  = [ctx.device newBufferWithLength:r2aa_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bWas  = [ctx.device newBufferWithLength:r2as_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD_aa = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD    = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGc   = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDh   = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bLam  = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];

        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGy = [ctx.device newBufferWithBytes:&g_y length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH  = [ctx.device newBufferWithBytes:&h length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSk = [ctx.device newBufferWithBytes:&spiky_const length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM  = [ctx.device newBufferWithBytes:&mass length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR  = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bA  = [ctx.device newBufferWithBytes:&alpha_inv_dt2 length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNs = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        uint32_t n_aa_total = n_active * n_active;
        uint32_t n_as_total = n_active * n_static;
        id<MTLBuffer> bNaaTot = [ctx.device newBufferWithBytes:&n_aa_total length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNasTot = [ctx.device newBufferWithBytes:&n_as_total length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        // sim_scale_inv buffer for predict_positions (default toy = 1.0)
        id<MTLBuffer> bSSI = [ctx.device newBufferWithBytes:&sim_scale_inv
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSS  = [ctx.device newBufferWithBytes:&sim_scale
            length:sizeof(float) options:MTLResourceStorageModeShared];
        // Identity unit-scale (used by update_velocities; same as bSS when use_pair)
        float one = 1.0f;
        id<MTLBuffer> bOne = [ctx.device newBufferWithBytes:&one length:sizeof(float) options:MTLResourceStorageModeShared];
        // External-accel scratch buffer (springs and/or pair_forces).
        id<MTLBuffer> bExtA = need_ext_accel
            ? [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared] : nil;
        // Spring bond buffers (when springs on).
        id<MTLBuffer> bBondIJ_sp = use_springs
            ? [ctx.device newBufferWithBytes:bond_ij_data
                  length:(size_t)n_bonds * 2 * sizeof(int32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bBondRest_sp = use_springs
            ? [ctx.device newBufferWithBytes:bond_rest_data
                  length:(size_t)n_bonds * sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bSpringK = use_springs
            ? [ctx.device newBufferWithBytes:&spring_K length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bNbonds_sp = use_springs
            ? [ctx.device newBufferWithBytes:&n_bonds length:sizeof(uint32_t)
                  options:MTLResourceStorageModeShared] : nil;
        // Floor constraint buffers (only when use_floor).
        id<MTLBuffer> bFloorY = use_floor
            ? [ctx.device newBufferWithBytes:&floor_y length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bRestitution = use_floor
            ? [ctx.device newBufferWithBytes:&restitution length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bClamped = use_floor
            ? [ctx.device newBufferWithLength:(size_t)n_active * sizeof(int32_t)
                  options:MTLResourceStorageModeShared] : nil;
        // M10 membrane buffers.
        id<MTLBuffer> bMembranes = use_membranes
            ? [ctx.device newBufferWithBytes:mem_tris
                  length:(size_t)n_membranes * 3 * sizeof(int32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bPmemIdx = use_membranes
            ? [ctx.device newBufferWithBytes:mem_pidx
                  length:(size_t)n_elastic * 7 * sizeof(int32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bMemCorr = use_membranes
            ? [ctx.device newBufferWithLength:pos_b
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bNelastic = use_membranes
            ? [ctx.device newBufferWithBytes:&n_elastic length:sizeof(uint32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bR0       = use_membranes
            ? [ctx.device newBufferWithBytes:&r0_mem length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bSortedS = (use_pair && n_static > 0)
            ? [ctx.device newBufferWithBytes:sorted_static_buf
                  length:(size_t)n_static * 3 * sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bCellStart = (use_pair && n_static > 0)
            ? [ctx.device newBufferWithBytes:cell_start_buf
                  length:(size_t)(n_cells + 1) * sizeof(int)
                  options:MTLResourceStorageModeShared] : nil;
        int grid_dim_packed[3] = {grid_dim_struct.x, grid_dim_struct.y, grid_dim_struct.z};
        float grid_origin_packed[3] = {grid_origin_struct.x, grid_origin_struct.y, grid_origin_struct.z};
        id<MTLBuffer> bGridDim = use_pair
            ? [ctx.device newBufferWithBytes:grid_dim_packed length:3*sizeof(int)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bGridOrigin = use_pair
            ? [ctx.device newBufferWithBytes:grid_origin_packed length:3*sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bViscPair = use_pair
            ? [ctx.device newBufferWithBytes:&visc_pair_coef length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bViscAmp  = use_pair
            ? [ctx.device newBufferWithBytes:&visc_amp length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bSurfAmp  = use_pair
            ? [ctx.device newBufferWithBytes:&surf_amp length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        // Seed bD with rho_rest so the first step's pair_forces has a
        // sensible 1/ρ factor (same convention as xpbd_step).
        if (use_pair) {
            float *dens0 = (float *)[bD contents];
            for (uint32_t i = 0; i < n_active; i++) dens0[i] = rho_rest;
        }

        for (uint32_t k = 0; k < K; k++) {
            // SAVE x_old, v_old (current state at start of step)
            float *step_state = state + (size_t)k * per_step_floats;
            memcpy(step_state, [bX contents], pos_b);
            memcpy(step_state + 3 * n_active, [bV contents], pos_b);

            memset([bLam contents], 0, s_b);

            // SAVE pair_density (the density used by THIS step's
            // pair_forces). For step 0 this is the rho_rest seed; for
            // step k>0 it's the density at the end of step k-1.
            if (use_pair) {
                memcpy(step_state + 11 * n_active, [bD contents], s_b);
            }

            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];

            // (0) External forces (visc + surf + springs) accumulate
            //     into ext_accel, then apply to velocity before predict.
            if (need_ext_accel) {
                memset([bExtA contents], 0, pos_b);  // zero accumulator
            }
            if (use_pair) {
                id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_pair];
                [e setBuffer:bX           offset:0 atIndex:0];
                [e setBuffer:bV           offset:0 atIndex:1];
                [e setBuffer:bSortedS     offset:0 atIndex:2];
                [e setBuffer:bCellStart   offset:0 atIndex:3];
                [e setBuffer:bD           offset:0 atIndex:4];
                [e setBuffer:bExtA        offset:0 atIndex:5];
                [e setBuffer:bH           offset:0 atIndex:6];
                [e setBuffer:bH2          offset:0 atIndex:7];
                [e setBuffer:bM           offset:0 atIndex:8];
                [e setBuffer:bSS          offset:0 atIndex:9];
                [e setBuffer:bViscPair    offset:0 atIndex:10];
                [e setBuffer:bViscAmp     offset:0 atIndex:11];
                [e setBuffer:bSurfAmp     offset:0 atIndex:12];
                [e setBuffer:bNa          offset:0 atIndex:13];
                [e setBuffer:bGridDim     offset:0 atIndex:14];
                [e setBuffer:bGridOrigin  offset:0 atIndex:15];
                [e dispatchThreads:MTLSizeMake(32, n_active, 1)
                    threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
                [e endEncoding];
            }
            if (use_springs) {
                id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_spring];
                [e setBuffer:bX            offset:0 atIndex:0];
                [e setBuffer:bBondIJ_sp    offset:0 atIndex:1];
                [e setBuffer:bBondRest_sp  offset:0 atIndex:2];
                [e setBuffer:bExtA         offset:0 atIndex:3];
                [e setBuffer:bSpringK      offset:0 atIndex:4];
                [e setBuffer:bNbonds_sp    offset:0 atIndex:5];
                [e setBuffer:bNa           offset:0 atIndex:6];
                [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [e endEncoding];
            }
            if (need_ext_accel) {
                id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_appext];
                [e setBuffer:bV    offset:0 atIndex:0];
                [e setBuffer:bExtA offset:0 atIndex:1];
                [e setBuffer:bDt   offset:0 atIndex:2];
                [e setBuffer:bNa   offset:0 atIndex:3];
                [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [e endEncoding];
            }

            // (1) predict
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pred];
              [e setBuffer:bX offset:0 atIndex:0]; [e setBuffer:bV offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bGy offset:0 atIndex:4]; [e setBuffer:bNa offset:0 atIndex:5];
              [e setBuffer:bSSI offset:0 atIndex:6];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (2) dist_aa
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_d_aa];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bR2aa offset:0 atIndex:1];
              [e setBuffer:bNa offset:0 atIndex:2];
              [e dispatchThreads:MTLSizeMake(n_active, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(16,16,1)];
              [e endEncoding]; }
            // (3) dist_as
            if (n_static > 0) {
              id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_d_as];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bSp offset:0 atIndex:1];
              [e setBuffer:bR2as offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
              [e setBuffer:bNs offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n_active, n_static, 1)
                  threadsPerThreadgroup:MTLSizeMake(16,16,1)];
              [e endEncoding];
            }
            // (4) Wpoly6 on aa
            { id<MTLBlitCommandEncoder> b = [cmd blitCommandEncoder];
              [b copyFromBuffer:bR2aa sourceOffset:0 toBuffer:bWaa destinationOffset:0 size:r2aa_b];
              [b endEncoding]; }
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_wp];
              [e setBuffer:bWaa offset:0 atIndex:0]; [e setBuffer:bH2 offset:0 atIndex:1];
              [e setBuffer:bP6 offset:0 atIndex:2]; [e setBuffer:bNaaTot offset:0 atIndex:3];
              [e dispatchThreads:MTLSizeMake(n_aa_total,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
              [e endEncoding]; }
            // (5) Wpoly6 on as
            if (n_static > 0) {
              { id<MTLBlitCommandEncoder> b = [cmd blitCommandEncoder];
                [b copyFromBuffer:bR2as sourceOffset:0 toBuffer:bWas destinationOffset:0 size:r2as_b];
                [b endEncoding]; }
              { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_wp];
                [e setBuffer:bWas offset:0 atIndex:0]; [e setBuffer:bH2 offset:0 atIndex:1];
                [e setBuffer:bP6 offset:0 atIndex:2]; [e setBuffer:bNasTot offset:0 atIndex:3];
                [e dispatchThreads:MTLSizeMake(n_as_total,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
                [e endEncoding]; }
            }
            // (6) rowsum aa → density_aa
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_rs];
              [e setBuffer:bWaa offset:0 atIndex:0]; [e setBuffer:bD_aa offset:0 atIndex:1];
              [e setBuffer:bM offset:0 atIndex:2];
              [e setBuffer:bNa offset:0 atIndex:3]; [e setBuffer:bNa offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(256, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(256,1,1)];
              [e endEncoding]; }
            // (7) rowsum as → density (will accumulate via add)
            if (n_static > 0) {
              id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_rs];
              [e setBuffer:bWas offset:0 atIndex:0]; [e setBuffer:bD offset:0 atIndex:1];
              [e setBuffer:bM offset:0 atIndex:2];
              [e setBuffer:bNs offset:0 atIndex:3]; [e setBuffer:bNa offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(256, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(256,1,1)];
              [e endEncoding];
              // density += density_aa
              { id<MTLComputeCommandEncoder> e2 = [cmd computeCommandEncoder];
                [e2 setComputePipelineState:pso_addin];
                [e2 setBuffer:bD offset:0 atIndex:0]; [e2 setBuffer:bD_aa offset:0 atIndex:1];
                [e2 setBuffer:bNa offset:0 atIndex:2];
                [e2 dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
                [e2 endEncoding]; }
            } else {
              // density = density_aa
              id<MTLBlitCommandEncoder> b = [cmd blitCommandEncoder];
              [b copyFromBuffer:bD_aa sourceOffset:0 toBuffer:bD destinationOffset:0 size:s_b];
              [b endEncoding];
            }
            // (8) density_constraint_grad
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_dgrad];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bSp offset:0 atIndex:1];
              [e setBuffer:bR2aa offset:0 atIndex:2]; [e setBuffer:bR2as offset:0 atIndex:3];
              [e setBuffer:bGc offset:0 atIndex:4]; [e setBuffer:bDh offset:0 atIndex:5];
              [e setBuffer:bH offset:0 atIndex:6]; [e setBuffer:bSk offset:0 atIndex:7];
              [e setBuffer:bM offset:0 atIndex:8]; [e setBuffer:bR offset:0 atIndex:9];
              [e setBuffer:bNa offset:0 atIndex:10]; [e setBuffer:bNs offset:0 atIndex:11];
              [e dispatchThreads:MTLSizeMake(256, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(256,1,1)];
              [e endEncoding]; }
            // (9) solve_density_constraint
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_solv];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bLam offset:0 atIndex:1];
              [e setBuffer:bD offset:0 atIndex:2]; [e setBuffer:bGc offset:0 atIndex:3];
              [e setBuffer:bDh offset:0 atIndex:4]; [e setBuffer:bR offset:0 atIndex:5];
              [e setBuffer:bM offset:0 atIndex:6]; [e setBuffer:bA offset:0 atIndex:7];
              [e setBuffer:bNa offset:0 atIndex:8];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (9b) Optional elastic floor constraint with mask (matches xpbd_step).
            if (use_floor) {
                id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_floor];
                [e setBuffer:bXp           offset:0 atIndex:0];
                [e setBuffer:bClamped      offset:0 atIndex:1];
                [e setBuffer:bFloorY       offset:0 atIndex:2];
                [e setBuffer:bNa           offset:0 atIndex:3];
                [e setBuffer:bRestitution  offset:0 atIndex:4];
                [e dispatchThreads:MTLSizeMake(n_active,1,1)
                    threadsPerThreadgroup:MTLSizeMake(64,1,1)];
                [e endEncoding];
            }
            // (10) update_vel: v_new = (xp - x_old) · sim_scale / dt
            // Runs BEFORE membrane apply — velocity is from post-floor pos
            // only. See xpbd_step.mm for rationale.
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uv];
              [e setBuffer:bV offset:0 atIndex:0]; [e setBuffer:bX offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bNa offset:0 atIndex:4];
              [e setBuffer:bSS offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (10b) M10 membrane interaction — position-only correction,
            // last op of the step (after update_velocities).
            if (use_membranes) {
                { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                  [e setComputePipelineState:pso_mem_clear];
                  [e setBuffer:bMemCorr offset:0 atIndex:0];
                  [e setBuffer:bNa      offset:0 atIndex:1];
                  [e dispatchThreads:MTLSizeMake(n_active,1,1)
                      threadsPerThreadgroup:MTLSizeMake(64,1,1)];
                  [e endEncoding]; }
                { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                  [e setComputePipelineState:pso_mem_acc];
                  [e setBuffer:bXp         offset:0 atIndex:0];
                  [e setBuffer:bMembranes  offset:0 atIndex:1];
                  [e setBuffer:bPmemIdx    offset:0 atIndex:2];
                  [e setBuffer:bMemCorr    offset:0 atIndex:3];
                  [e setBuffer:bNa         offset:0 atIndex:4];
                  [e setBuffer:bNelastic   offset:0 atIndex:5];
                  [e setBuffer:bR0         offset:0 atIndex:6];
                  [e dispatchThreads:MTLSizeMake(n_active,1,1)
                      threadsPerThreadgroup:MTLSizeMake(64,1,1)];
                  [e endEncoding]; }
                { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                  [e setComputePipelineState:pso_mem_apply];
                  [e setBuffer:bXp      offset:0 atIndex:0];
                  [e setBuffer:bMemCorr offset:0 atIndex:1];
                  [e setBuffer:bNa      offset:0 atIndex:2];
                  [e dispatchThreads:MTLSizeMake(n_active,1,1)
                      threadsPerThreadgroup:MTLSizeMake(64,1,1)];
                  [e endEncoding]; }
            }
            [cmd commit]; [cmd waitUntilCompleted];

            // SAVE density, grad_C, denom_helper for backward
            float *p = step_state + 6 * n_active;
            memcpy(p, [bD contents], s_b);                            // density
            memcpy(p + n_active, [bGc contents], pos_b);              // grad_C
            memcpy(p + n_active + 3 * n_active, [bDh contents], s_b); // denom_h
            // Save floor mask (after pair_density slot if use_pair).
            if (use_floor) {
                size_t floor_off = (size_t)(11 + (use_pair ? 1 : 0)) * n_active;
                memcpy(step_state + floor_off, [bClamped contents],
                       (size_t)n_active * sizeof(int32_t));
            }
            // Save mem_corr (after pair_density + floor mask slots).
            if (use_membranes) {
                size_t mc_off = (size_t)(11
                                + (use_pair ? 1 : 0)
                                + (use_floor ? 1 : 0)) * n_active;
                memcpy(step_state + mc_off, [bMemCorr contents], pos_b);
            }
            // Save x_post in trajectory
            memcpy(traj + (size_t)(k + 1) * n_active * 3,
                   [bXp contents], pos_b);
            // Advance: x_old ← x_post
            memcpy([bX contents], [bXp contents], pos_b);
        }

        // Write outputs: state buffer, trajectory, final velocity
        FILE *fs = fopen(argv[14], "wb");
        fwrite(state, sizeof(float), (size_t)K * per_step_floats, fs);
        fwrite(traj, sizeof(float), (size_t)(K + 1) * n_active * 3, fs);
        fwrite([bV contents], 1, pos_b, fs);
        fclose(fs);
    }
    free(pos0); free(vel0); free(pos_static); free(state); free(traj);
    if (sorted_static_buf) free(sorted_static_buf);
    if (cell_start_buf)    free(cell_start_buf);
    if (bond_ij_data) free(bond_ij_data);
    if (bond_rest_data) free(bond_rest_data);
    if (mem_tris) free(mem_tris);
    if (mem_pidx) free(mem_pidx);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// Multi-step backward driver. Walks K steps in reverse, calling all
// the M9 backward kernels per step, accumulating ∂L/∂(rho_rest, x_init,
// v_init).
//
// Inputs:
//   state file (from xpbd_full_fwd)
//   ∂L/∂x_final (input gradient seed)
// Outputs:
//   ∂L/∂x_init, ∂L/∂v_init, ∂L/∂rho_rest (scalar)
// ──────────────────────────────────────────────────────────────────────
int run_xpbd_full_bwd(int argc, char **argv) {
    if (argc < 17 || argc > 31) {
        fprintf(stderr,
            "usage: sib_metal xpbd_full_bwd "
            "<n_active> <n_static> <K> <h> <mass> <rho_rest> <dt> "
            "<gravity_y> <alpha_density> "
            "<state_in.bin> <pos_static.bin> <grad_x_final.bin> "
            "<grad_x_init_out.bin> <grad_v_init_out.bin> <grad_rho_out.bin> "
            "[sim_scale] [visc_pair_coef] [spring_K] [bonds.bin] "
            "[grad_spring_K_out.bin] [grad_visc_K_out.bin] [floor_y] "
            "[grad_alpha_dens_out.bin] [restitution] "
            "[n_membranes] [n_elastic] [r0] [membranes.bin] [pmem_index.bin]\n"
            "       (must match the xpbd_full_fwd args used to produce the "
            "state file)\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    uint32_t K        = (uint32_t)atoi(argv[4]);
    float h           = (float)atof(argv[5]);
    float mass        = (float)atof(argv[6]);
    float rho_rest    = (float)atof(argv[7]);
    float dt          = (float)atof(argv[8]);
    float g_y         = (float)atof(argv[9]);
    float alpha_dens  = (float)atof(argv[10]);
    // Optional 17th arg: sim_scale; 18th: visc_pair_coef. Match xpbd_full_fwd.
    float sim_scale   = (argc >= 18) ? (float)atof(argv[17]) : 1.0f;
    float sim_scale_inv = 1.0f / sim_scale;
    float visc_pair_coef = (argc >= 19) ? (float)atof(argv[18]) : 0.0f;
    bool use_pair = visc_pair_coef != 0.0f;
    float spring_K = (argc >= 20) ? (float)atof(argv[19]) : 0.0f;
    const char *bonds_path = (argc >= 21) ? argv[20] : NULL;
    bool use_springs = spring_K != 0.0f;
    if (use_springs && !bonds_path) {
        fprintf(stderr, "spring_K > 0 requires bonds.bin path as 21st arg\n");
        return 1;
    }
    // Optional floor_y at argv[23] (positions 21, 22 reserved for K/visc grad outs).
    bool use_floor = (argc >= 24);
    float floor_y = use_floor ? (float)atof(argv[23]) : 0.0f;
    // Optional restitution at argv[25] (after grad_alpha_dens output at argv[24]).
    // Default 0 = inelastic (legacy). Pass 0..1 to enable elastic floor.
    float restitution = (argc >= 26) ? (float)atof(argv[25]) : 0.0f;
    // M10 membrane backward (must match xpbd_full_fwd args).
    uint32_t n_membranes = (argc >= 27) ? (uint32_t)atoi(argv[26]) : 0u;
    uint32_t n_elastic   = (argc >= 28) ? (uint32_t)atoi(argv[27]) : 0u;
    float r0_mem         = (argc >= 29) ? (float)atof(argv[28]) : 0.0f;
    const char *path_membranes = (argc >= 30) ? argv[29] : NULL;
    const char *path_pmem_idx  = (argc >= 31) ? argv[30] : NULL;
    bool use_membranes = (n_membranes > 0) && (n_elastic > 0)
                         && path_membranes && path_pmem_idx
                         && (r0_mem > 0.0f);
    // Load bonds for springs (mirrors xpbd_full_fwd loader).
    uint32_t n_bonds = 0;
    int32_t *bond_ij_data = NULL;
    float *bond_rest_data = NULL;
    if (use_springs) {
        FILE *fb = fopen(bonds_path, "rb");
        if (!fb) { fprintf(stderr, "open %s\n", bonds_path); return 1; }
        fseek(fb, 0, SEEK_END);
        long sz = ftell(fb);
        fseek(fb, 0, SEEK_SET);
        n_bonds = (uint32_t)(sz / 16);
        uint8_t *raw = (uint8_t *)malloc((size_t)n_bonds * 16);
        fread(raw, 1, (size_t)n_bonds * 16, fb);
        fclose(fb);
        bond_ij_data = (int32_t *)malloc((size_t)n_bonds * 2 * sizeof(int32_t));
        bond_rest_data = (float *)malloc((size_t)n_bonds * sizeof(float));
        for (uint32_t b = 0; b < n_bonds; b++) {
            memcpy(&bond_ij_data[b * 2], raw + b * 16, 8);
            memcpy(&bond_rest_data[b], raw + b * 16 + 8, 4);
        }
        free(raw);
    }

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));
    float alpha_inv_dt2 = alpha_dens / (dt * dt);
    // Pair-force amps (fp64 for fp32-underflow safety).
    double h_scaled = (double)h * (double)sim_scale;
    double h_s6 = pow(h_scaled, 6.0);
    double h_s9 = pow(h_scaled, 9.0);
    double divgradWvisco = 45.0 / (M_PI * h_s6);
    float visc_amp = (float)(1.5 * (double)mass * divgradWvisco
                              * pow((double)sim_scale, 3.0));
    double wpoly6_si = 315.0 / (64.0 * M_PI * h_s9);
    float surf_amp = (float)(-1.7e-9 * (double)mass * wpoly6_si
                              * (double)sim_scale / (double)mass);

    auto rd = ^(const char *p, size_t n_floats) {
        FILE *f = fopen(p, "rb");
        if (!f) { fprintf(stderr, "open %s\n", p); exit(1); }
        float *b = (float *)malloc(n_floats * sizeof(float));
        fread(b, sizeof(float), n_floats, f); fclose(f); return b;
    };
    // State file layout (from xpbd_full_fwd):
    //   [K × per_step_floats] state
    //   [(K+1) × n*3] traj
    //   [n*3] vel_final
    int extra_per_step = (use_pair ? 1 : 0) + (use_floor ? 1 : 0)
                       + (use_membranes ? 3 : 0);
    size_t per_step_floats = (size_t)n_active * (3 + 3 + 1 + 3 + 1 + extra_per_step);
    size_t state_size = (size_t)K * per_step_floats
                      + (size_t)(K + 1) * n_active * 3
                      + (size_t)n_active * 3;
    float *all = rd(argv[11], state_size);
    float *state = all;
    float *traj = all + (size_t)K * per_step_floats;
    // (vel_final at all + ... + (K+1)*n*3 — not used here)

    float *pos_static = rd(argv[12], (size_t)n_static * 3);
    float *grad_x_fin = rd(argv[13], (size_t)n_active * 3);

    // Load membrane topology if enabled.
    int32_t *mem_tris = NULL;
    int32_t *mem_pidx = NULL;
    if (use_membranes) {
        size_t tris_bytes = (size_t)n_membranes * 3 * sizeof(int32_t);
        size_t pidx_bytes = (size_t)n_elastic * 7 * sizeof(int32_t);
        FILE *fa = fopen(path_membranes, "rb");
        if (!fa) { fprintf(stderr, "open %s\n", path_membranes); return 1; }
        mem_tris = (int32_t *)malloc(tris_bytes);
        if (fread(mem_tris, 1, tris_bytes, fa) != tris_bytes) {
            fprintf(stderr, "short read on %s\n", path_membranes); return 1;
        }
        fclose(fa);
        FILE *fb2 = fopen(path_pmem_idx, "rb");
        if (!fb2) { fprintf(stderr, "open %s\n", path_pmem_idx); return 1; }
        mem_pidx = (int32_t *)malloc(pidx_bytes);
        if (fread(mem_pidx, 1, pidx_bytes, fb2) != pidx_bytes) {
            fprintf(stderr, "short read on %s\n", path_pmem_idx); return 1;
        }
        fclose(fb2);
    }

    size_t pos_b = (size_t)n_active * 3 * sizeof(float);
    size_t s_b = (size_t)n_active * sizeof(float);

    // Static spatial grid buffers — outer scope for cleanup.
    float *sorted_static_buf = NULL;
    int *cell_start_buf = NULL;

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        // Forward kernels (recompute distances and predict x_pre):
        id<MTLComputePipelineState> pso_pred  = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_d_aa  = make_pso(ctx, "dist_active_active");
        id<MTLComputePipelineState> pso_d_as  = make_pso(ctx, "dist_active_static");
        // Backward kernels:
        id<MTLComputePipelineState> pso_uvbw  = make_pso(ctx, "update_velocities_backward");
        id<MTLComputePipelineState> pso_solv_bw = make_pso(ctx, "solve_density_constraint_backward");
        id<MTLComputePipelineState> pso_dgrad_bw = make_pso(ctx, "density_constraint_grad_backward");
        id<MTLComputePipelineState> pso_rs_bw = make_pso(ctx, "rowsum_density_backward");
        id<MTLComputePipelineState> pso_wp_bw = make_pso(ctx, "wpoly6_inplace_backward");
        id<MTLComputePipelineState> pso_d_aa_bw = make_pso(ctx, "dist_active_active_backward");
        id<MTLComputePipelineState> pso_d_as_bw = make_pso(ctx, "dist_active_static_backward");
        id<MTLComputePipelineState> pso_pred_bw = make_pso(ctx, "predict_positions_backward");
        // Pair-force backward (only if use_pair).
        id<MTLComputePipelineState> pso_pair_bw = use_pair
            ? make_pso(ctx, "pair_forces_grid_backward") : nil;
        // Spring-force backward (only if use_springs).
        id<MTLComputePipelineState> pso_spring_bw = use_springs
            ? make_pso(ctx, "spring_bonds_force_backward") : nil;
        // Param-gradient kernels for SGD on spring_K and visc_pair_coef.
        id<MTLComputePipelineState> pso_spring_K_part = use_springs
            ? make_pso(ctx, "spring_K_partial") : nil;
        id<MTLComputePipelineState> pso_visc_K_part   = use_pair
            ? make_pso(ctx, "visc_K_partial") : nil;
        // Floor backward (matches forward solve_floor_constraint_with_mask).
        id<MTLComputePipelineState> pso_floor_bw = use_floor
            ? make_pso(ctx, "solve_floor_constraint_backward") : nil;
        // M10 backward kernels.
        id<MTLComputePipelineState> pso_mem_apply_bw = use_membranes
            ? make_pso(ctx, "apply_membrane_correction_backward") : nil;
        id<MTLComputePipelineState> pso_mem_acc_bw = use_membranes
            ? make_pso(ctx, "accumulate_membrane_correction_backward") : nil;

        // Build static spatial grid for pair_forces backward (matches forward).
        GridDim3 grid_dim_struct = {0, 0, 0};
        GridOrigin3 grid_origin_struct = {0, 0, 0};
        int n_cells = 0;
        if (use_pair && n_static > 0) {
            build_static_grid(pos_static, n_static, h,
                              &sorted_static_buf, &cell_start_buf,
                              &grid_dim_struct, &grid_origin_struct, &n_cells);
        }

        // Persistent buffers.
        id<MTLBuffer> bX  = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bV  = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bXp = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSp = [ctx.device newBufferWithBytes:pos_static
            length:(size_t)n_static * 3 * sizeof(float) options:MTLResourceStorageModeShared];
        size_t r2aa_b = (size_t)n_active * n_active * sizeof(float);
        size_t r2as_b = (size_t)n_active * n_static * sizeof(float);
        id<MTLBuffer> bR2aa = [ctx.device newBufferWithLength:r2aa_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2as = [ctx.device newBufferWithLength:r2as_b options:MTLResourceStorageModeShared];

        // Per-step inputs from state:
        id<MTLBuffer> bD  = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGc = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDh = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bLam = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        memset([bLam contents], 0, s_b);  // λ_pre always 0

        // Running gradients:
        id<MTLBuffer> bGx_running = [ctx.device newBufferWithBytes:grad_x_fin length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGv_running = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        memset([bGv_running contents], 0, pos_b);

        // M10 membrane buffers.
        id<MTLBuffer> bMembranes_bw = use_membranes
            ? [ctx.device newBufferWithBytes:mem_tris
                  length:(size_t)n_membranes * 3 * sizeof(int32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bPmemIdx_bw = use_membranes
            ? [ctx.device newBufferWithBytes:mem_pidx
                  length:(size_t)n_elastic * 7 * sizeof(int32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bGmemCorr = use_membranes
            ? [ctx.device newBufferWithLength:pos_b
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bXpre = use_membranes
            ? [ctx.device newBufferWithLength:pos_b
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bNelastic_bw = use_membranes
            ? [ctx.device newBufferWithBytes:&n_elastic length:sizeof(uint32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bR0_bw = use_membranes
            ? [ctx.device newBufferWithBytes:&r0_mem length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;

        // Per-step working gradients:
        id<MTLBuffer> bGv_in = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGx_pred = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGx_old_new = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGv_old_new = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        // Outputs of solve_density_constraint_backward:
        id<MTLBuffer> bGlam_pre = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGdens = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGgC = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGdh = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGrho_per = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGalpha_per = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];
        // Density chain backward intermediates:
        id<MTLBuffer> bGW_or_Gr2_aa = [ctx.device newBufferWithLength:r2aa_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGW_or_Gr2_as = [ctx.device newBufferWithLength:r2as_b options:MTLResourceStorageModeShared];
        // For predict_backward:
        id<MTLBuffer> bGgy_per = [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared];

        // Constants.
        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGy = [ctx.device newBufferWithBytes:&g_y length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH  = [ctx.device newBufferWithBytes:&h length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSk = [ctx.device newBufferWithBytes:&spiky_const length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM  = [ctx.device newBufferWithBytes:&mass length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR  = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bA  = [ctx.device newBufferWithBytes:&alpha_inv_dt2 length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNs = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        uint32_t n_aa_total = n_active * n_active;
        uint32_t n_as_total = n_active * n_static;
        id<MTLBuffer> bNaaTot = [ctx.device newBufferWithBytes:&n_aa_total length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNasTot = [ctx.device newBufferWithBytes:&n_as_total length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        // sim_scale_inv buffer for predict_positions (default toy = 1.0)
        id<MTLBuffer> bSSI = [ctx.device newBufferWithBytes:&sim_scale_inv
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSS  = [ctx.device newBufferWithBytes:&sim_scale
            length:sizeof(float) options:MTLResourceStorageModeShared];
        // Identity unit-scale (toy convention: pos & vel share unit system).
        float one = 1.0f;
        id<MTLBuffer> bOne = [ctx.device newBufferWithBytes:&one length:sizeof(float) options:MTLResourceStorageModeShared];
        // External-accel gradient scratch (springs and/or pair_forces).
        id<MTLBuffer> bGext = (use_pair || use_springs)
            ? [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared] : nil;
        // Spring bond buffers (when springs on).
        id<MTLBuffer> bBondIJ_sp = use_springs
            ? [ctx.device newBufferWithBytes:bond_ij_data
                  length:(size_t)n_bonds * 2 * sizeof(int32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bBondRest_sp = use_springs
            ? [ctx.device newBufferWithBytes:bond_rest_data
                  length:(size_t)n_bonds * sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bSpringK = use_springs
            ? [ctx.device newBufferWithBytes:&spring_K length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bNbonds_sp = use_springs
            ? [ctx.device newBufferWithBytes:&n_bonds length:sizeof(uint32_t)
                  options:MTLResourceStorageModeShared] : nil;
        // Per-particle scratch for analytic param gradients (SGD on
        // spring_K / visc_pair_coef). Host sums per-particle partials
        // each step into total_grad_spring_K / total_grad_visc_K.
        id<MTLBuffer> bSpringPart = use_springs
            ? [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bViscPart = use_pair
            ? [ctx.device newBufferWithLength:s_b options:MTLResourceStorageModeShared] : nil;
        // Floor mask + scratch grad buffer (only when use_floor).
        id<MTLBuffer> bClamped_bw = use_floor
            ? [ctx.device newBufferWithLength:(size_t)n_active * sizeof(int32_t)
                  options:MTLResourceStorageModeShared] : nil;
        // packed_float3 buffer for grad_pos_pre (ZERO'd then accumulated);
        // sized as n_active * 3 floats (pos_b), NOT s_b which is just n.
        id<MTLBuffer> bGfloorScr = use_floor
            ? [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bRestitution_bw = use_floor
            ? [ctx.device newBufferWithBytes:&restitution length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bSortedS = (use_pair && n_static > 0)
            ? [ctx.device newBufferWithBytes:sorted_static_buf
                  length:(size_t)n_static * 3 * sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bCellStart = (use_pair && n_static > 0)
            ? [ctx.device newBufferWithBytes:cell_start_buf
                  length:(size_t)(n_cells + 1) * sizeof(int)
                  options:MTLResourceStorageModeShared] : nil;
        int grid_dim_packed[3] = {grid_dim_struct.x, grid_dim_struct.y, grid_dim_struct.z};
        float grid_origin_packed[3] = {grid_origin_struct.x, grid_origin_struct.y, grid_origin_struct.z};
        id<MTLBuffer> bGridDim = use_pair
            ? [ctx.device newBufferWithBytes:grid_dim_packed length:3*sizeof(int)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bGridOrigin = use_pair
            ? [ctx.device newBufferWithBytes:grid_origin_packed length:3*sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bViscPair = use_pair
            ? [ctx.device newBufferWithBytes:&visc_pair_coef length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bViscAmp = use_pair
            ? [ctx.device newBufferWithBytes:&visc_amp length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bSurfAmp = use_pair
            ? [ctx.device newBufferWithBytes:&surf_amp length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;

        // Total parameter-gradient accumulators (host-side scalars).
        float total_grad_rho = 0.0f;
        float total_grad_spring_K = 0.0f;  // Σ_steps Σ_i <∂a_i/∂K, ga_i>
        float total_grad_visc_K = 0.0f;    // Σ_steps Σ_i <∂a_i/∂visc_K, ga_i>
        float total_grad_alpha_inv_dt2 = 0.0f;  // Σ_steps Σ_i ∂L/∂A from density solve

        // Truncated BPTT: gradients through chaotic dynamics explode
        // exponentially per step (~3× for our cube-spring system).
        // Chain back only the last `max_bw_steps` steps if BWD_TBPTT is
        // set; default = K (no truncation, may overflow at large K).
        int32_t max_bw_steps = K;
        const char *tbptt_env = getenv("BWD_TBPTT");
        if (tbptt_env) {
            max_bw_steps = atoi(tbptt_env);
            if (max_bw_steps <= 0 || max_bw_steps > (int32_t)K) max_bw_steps = K;
            fprintf(stderr, "[xpbd_full_bwd] truncated BPTT: %d / %u steps\n",
                    max_bw_steps, K);
        }
        int32_t k_stop = (int32_t)K - max_bw_steps;  // walk K-1 down to k_stop

        // Per-step gradient clipping. Each chain step amplifies position
        // gradients by ~sim_scale_inv (≈ 1.35e5) because update_vel
        // backward divides by dt and predict backward multiplies by
        // dt·sim_scale_inv. With Sibernetic's mass=2e-12 / sim_scale=7.4e-6,
        // gradients overflow to NaN within 5–10 chain steps. Clipping the
        // running x and v gradients per step yields a biased-but-bounded
        // gradient that still points in the right direction for SGD.
        // Set BWD_CLIP_NORM=0 (or unset) to disable; default = 1e3.
        float clip_norm_x = 1e3f;
        float clip_norm_v = 1e3f;
        const char *clip_env = getenv("BWD_CLIP_NORM");
        if (clip_env) {
            clip_norm_x = (float)atof(clip_env);
            clip_norm_v = clip_norm_x;
            fprintf(stderr, "[xpbd_full_bwd] per-step gradient clipping: "
                    "|grad_x|, |grad_v| L2 capped at %.3e\n", clip_norm_x);
        }

        // Walk K steps backward (or fewer if truncated BPTT).
        for (int32_t k = (int32_t)K - 1; k >= k_stop; k--) {
            float *step_state = state + (size_t)k * per_step_floats;
            // Load saved state.
            memcpy([bX contents], step_state, pos_b);                  // x_old
            memcpy([bV contents], step_state + 3 * n_active, pos_b);   // v_old
            memcpy([bD contents], step_state + 6 * n_active, s_b);     // density
            memcpy([bGc contents], step_state + 7 * n_active, pos_b);  // grad_C
            memcpy([bDh contents], step_state + 10 * n_active, s_b);   // denom_h

            // Recompute x_pre = x_old + dt²·g + dt·v_old (same as predict).
            id<MTLCommandBuffer> cmd1 = [ctx.queue commandBuffer];
            { id<MTLComputeCommandEncoder> e = [cmd1 computeCommandEncoder];
              [e setComputePipelineState:pso_pred];
              [e setBuffer:bX offset:0 atIndex:0]; [e setBuffer:bV offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bGy offset:0 atIndex:4]; [e setBuffer:bNa offset:0 atIndex:5];
              [e setBuffer:bOne offset:0 atIndex:6];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // Recompute r2_aa, r2_as for the density chain backward.
            { id<MTLComputeCommandEncoder> e = [cmd1 computeCommandEncoder];
              [e setComputePipelineState:pso_d_aa];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bR2aa offset:0 atIndex:1];
              [e setBuffer:bNa offset:0 atIndex:2];
              [e dispatchThreads:MTLSizeMake(n_active, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(16,16,1)];
              [e endEncoding]; }
            if (n_static > 0) {
              id<MTLComputeCommandEncoder> e = [cmd1 computeCommandEncoder];
              [e setComputePipelineState:pso_d_as];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bSp offset:0 atIndex:1];
              [e setBuffer:bR2as offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
              [e setBuffer:bNs offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n_active, n_static, 1)
                  threadsPerThreadgroup:MTLSizeMake(16,16,1)];
              [e endEncoding];
            }
            [cmd1 commit]; [cmd1 waitUntilCompleted];

            memset([bGx_old_new contents], 0, pos_b);
            memset([bGv_old_new contents], 0, pos_b);
            memset([bGx_pred contents], 0, pos_b);
            // Move running v gradient into bGv_in for update_vel_bwd.
            memcpy([bGv_in contents], [bGv_running contents], pos_b);

            // (a-pre) Floor backward — undo forward's last kernel
            //         (solve_floor_constraint_with_mask) before update_v_bw.
            //         Clamping kills y-component of gradient on clamped particles.
            if (use_floor) {
                size_t floor_off = (size_t)(11 + (use_pair ? 1 : 0)) * n_active;
                memcpy([bClamped_bw contents],
                       step_state + floor_off,
                       (size_t)n_active * sizeof(int32_t));
                id<MTLCommandBuffer> cmdF = [ctx.queue commandBuffer];
                id<MTLComputeCommandEncoder> e = [cmdF computeCommandEncoder];
                [e setComputePipelineState:pso_floor_bw];
                // input grad_pos_post = bGx_running; output overwrites it
                // by writing into a temp and copying back. Since the kernel
                // ACCUMULATES into grad_pos_pre, we zero a scratch buffer
                // (bGfloorScr is reused as temp), then copy back.
                memset([bGfloorScr contents], 0, pos_b);
                [e setBuffer:bGx_running   offset:0 atIndex:0];  // grad_pos_post (input)
                // Note: backward kernel signature is grad_pos_post + grad_pos_pre
                // + grad_floor + clamped + n. We need grad_pos_pre as separate
                // accumulator buffer so we use bGx_pred (will be reset later
                // before density-bwd anyway). Actually safer: use a dedicated
                // scratch — overwrite bGx_running with the corrected gradient.
                [e setBuffer:bGfloorScr    offset:0 atIndex:1];  // grad_pos_pre (accum, zero'd)
                [e setBuffer:bGx_pred      offset:0 atIndex:2];  // grad_floor per-particle (we discard)
                [e setBuffer:bClamped_bw   offset:0 atIndex:3];
                [e setBuffer:bNa           offset:0 atIndex:4];
                [e setBuffer:bRestitution_bw offset:0 atIndex:5];
                [e dispatchThreads:MTLSizeMake(n_active,1,1)
                    threadsPerThreadgroup:MTLSizeMake(64,1,1)];
                [e endEncoding];
                [cmdF commit]; [cmdF waitUntilCompleted];
                // Replace bGx_running ← bGfloorScr (the floor-aware gradient).
                memcpy([bGx_running contents], [bGfloorScr contents], pos_b);
                // Reset bGx_pred (we used it as scratch for grad_floor).
                memset([bGx_pred contents], 0, pos_b);
            }

            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            // (a-mem) M10 membrane backward — runs FIRST in reverse since
            //         membrane apply is now the LAST forward op of the
            //         step (after update_velocities). bGx_running enters
            //         the loop body as ∂L/∂x_post_step = ∂L/∂x_post_apply.
            //         apply_bw is identity for x (∂L/∂x_pre_apply ←
            //         ∂L/∂x_post_apply, in place). accumulate_bw chains
            //         ∂L/∂mem_corr through the projection and adds to
            //         ∂L/∂x_pre_acc (which lives in the same bGx_running
            //         buffer). After: bGx_running = ∂L/∂x_post_floor,
            //         which is what the rest of the chain expects.
            if (use_membranes) {
                size_t mc_off = (size_t)(11
                                + (use_pair ? 1 : 0)
                                + (use_floor ? 1 : 0)) * n_active;
                float *mc_step  = step_state + mc_off;
                float *x_post   = traj + (size_t)(k + 1) * n_active * 3;
                float *xpre_buf = (float *)[bXpre contents];
                for (uint32_t i = 0; i < n_active * 3u; i++) {
                    xpre_buf[i] = x_post[i] - mc_step[i];
                }
                memcpy([bGmemCorr contents], [bGx_running contents], pos_b);

                id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_mem_acc_bw];
                [e setBuffer:bXpre          offset:0 atIndex:0];
                [e setBuffer:bMembranes_bw  offset:0 atIndex:1];
                [e setBuffer:bPmemIdx_bw    offset:0 atIndex:2];
                [e setBuffer:bGmemCorr      offset:0 atIndex:3];
                [e setBuffer:bGx_running    offset:0 atIndex:4];   // accumulate
                [e setBuffer:bNa            offset:0 atIndex:5];
                [e setBuffer:bNelastic_bw   offset:0 atIndex:6];
                [e setBuffer:bR0_bw         offset:0 atIndex:7];
                [e dispatchThreads:MTLSizeMake(n_active,1,1)
                    threadsPerThreadgroup:MTLSizeMake(64,1,1)];
                [e endEncoding];
                // Sync after membrane_bw so subsequent kernels see the
                // updated bGx_running.
                [cmd commit]; [cmd waitUntilCompleted];
                cmd = [ctx.queue commandBuffer];
            }
            // (a) update_vel backward: ∂L/∂v_new (=bGv_in) → +bGx_running (∂L/∂x_post_floor), +bGx_old_new (∂L/∂x_old)
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uvbw];
              [e setBuffer:bGv_in offset:0 atIndex:0];
              [e setBuffer:bGx_running offset:0 atIndex:1];
              [e setBuffer:bGx_old_new offset:0 atIndex:2];
              [e setBuffer:bDt offset:0 atIndex:3]; [e setBuffer:bNa offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (b) solve_density_constraint backward
            // ∂L/∂λ_post = 0 (single-iter, no propagation)
            // We use bGdens as a SCRATCH for the zeros.
            // Actually solve bw expects buf at index 5 to be ∂L/∂λ_post — pass bDh re-purposed? No, need a zero buffer.
            // Re-use bGdh for this purpose (we'll overwrite it with output ∂L/∂dh anyway).
            // But we need a SEPARATE zero buffer. Use a small temp.
            // Quick fix: just zero bGdh first and pass it. But bGdh gets overwritten by the kernel.
            // Cleaner: allocate a separate zero buffer for ∂L/∂λ_post.
            // Hack: bLam itself is zeros (we set above) — reuse it.
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_solv_bw];
              [e setBuffer:bD offset:0 atIndex:0];     [e setBuffer:bGc offset:0 atIndex:1];
              [e setBuffer:bDh offset:0 atIndex:2];    [e setBuffer:bLam offset:0 atIndex:3];
              [e setBuffer:bGx_running offset:0 atIndex:4];  // ∂L/∂x_post (input)
              [e setBuffer:bLam offset:0 atIndex:5];   // ∂L/∂λ_post = 0 (reuse bLam as zeros)
              [e setBuffer:bGx_pred offset:0 atIndex:6];  // ∂L/∂x_pre (accum)
              [e setBuffer:bGlam_pre offset:0 atIndex:7];
              [e setBuffer:bGdens offset:0 atIndex:8]; [e setBuffer:bGgC offset:0 atIndex:9];
              [e setBuffer:bGdh offset:0 atIndex:10];  [e setBuffer:bGrho_per offset:0 atIndex:11];
              [e setBuffer:bGalpha_per offset:0 atIndex:12];   // ∂L/∂A per particle (host sums)
              [e setBuffer:bR offset:0 atIndex:13];    [e setBuffer:bM offset:0 atIndex:14];
              [e setBuffer:bA offset:0 atIndex:15];    [e setBuffer:bNa offset:0 atIndex:16];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (c) density_constraint_grad backward: contributes ∂L/∂x via grad_C and denom_helper
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_dgrad_bw];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bSp offset:0 atIndex:1];
              [e setBuffer:bR2aa offset:0 atIndex:2]; [e setBuffer:bR2as offset:0 atIndex:3];
              [e setBuffer:bGgC offset:0 atIndex:4]; [e setBuffer:bGdh offset:0 atIndex:5];
              [e setBuffer:bGx_pred offset:0 atIndex:6];   // accumulate
              [e setBuffer:bH offset:0 atIndex:7]; [e setBuffer:bSk offset:0 atIndex:8];
              [e setBuffer:bM offset:0 atIndex:9]; [e setBuffer:bR offset:0 atIndex:10];
              [e setBuffer:bNa offset:0 atIndex:11]; [e setBuffer:bNs offset:0 atIndex:12];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // (d) density chain backward: ∂L/∂density → ∂L/∂x_pre (additional contribution)
            //     Pipeline: rowsum_bwd → wpoly6_bwd → dist_aa_bwd
            //               rowsum_bwd → wpoly6_bwd → dist_as_bwd
            //     Both contribute to bGx_pred (same buffer, multiple paths add up via dist_*_backward kernels).
            // density_aa branch:
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_rs_bw];
              [e setBuffer:bGdens offset:0 atIndex:0];
              [e setBuffer:bGW_or_Gr2_aa offset:0 atIndex:1];
              [e setBuffer:bM offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
              [e dispatchThreads:MTLSizeMake(n_active, n_active, 1)
                  threadsPerThreadgroup:MTLSizeMake(16,16,1)];
              [e endEncoding]; }
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_wp_bw];
              [e setBuffer:bR2aa offset:0 atIndex:0]; [e setBuffer:bGW_or_Gr2_aa offset:0 atIndex:1];
              [e setBuffer:bH2 offset:0 atIndex:2]; [e setBuffer:bP6 offset:0 atIndex:3];
              [e setBuffer:bNaaTot offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n_aa_total,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
              [e endEncoding]; }
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_d_aa_bw];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bGW_or_Gr2_aa offset:0 atIndex:1];
              [e setBuffer:bGx_pred offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // density_as branch:
            if (n_static > 0) {
              { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_rs_bw];
                [e setBuffer:bGdens offset:0 atIndex:0];
                [e setBuffer:bGW_or_Gr2_as offset:0 atIndex:1];
                [e setBuffer:bM offset:0 atIndex:2]; [e setBuffer:bNs offset:0 atIndex:3];
                [e dispatchThreads:MTLSizeMake(n_static, n_active, 1)
                    threadsPerThreadgroup:MTLSizeMake(16,16,1)];
                [e endEncoding]; }
              { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_wp_bw];
                [e setBuffer:bR2as offset:0 atIndex:0]; [e setBuffer:bGW_or_Gr2_as offset:0 atIndex:1];
                [e setBuffer:bH2 offset:0 atIndex:2]; [e setBuffer:bP6 offset:0 atIndex:3];
                [e setBuffer:bNasTot offset:0 atIndex:4];
                [e dispatchThreads:MTLSizeMake(n_as_total,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
                [e endEncoding]; }
              { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
                [e setComputePipelineState:pso_d_as_bw];
                [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bSp offset:0 atIndex:1];
                [e setBuffer:bGW_or_Gr2_as offset:0 atIndex:2];
                [e setBuffer:bGx_pred offset:0 atIndex:3];
                [e setBuffer:bNa offset:0 atIndex:4]; [e setBuffer:bNs offset:0 atIndex:5];
                [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
                [e endEncoding]; }
            }
            // (e) predict_positions backward: ∂L/∂x_pred → +bGx_old_new, +bGv_old_new, ∂L/∂g_y per particle
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pred_bw];
              [e setBuffer:bGx_pred offset:0 atIndex:0];
              [e setBuffer:bGx_old_new offset:0 atIndex:1];
              [e setBuffer:bGv_old_new offset:0 atIndex:2];
              [e setBuffer:bGgy_per offset:0 atIndex:3];
              [e setBuffer:bDt offset:0 atIndex:4]; [e setBuffer:bNa offset:0 atIndex:5];
              [e setBuffer:bSSI offset:0 atIndex:6];
              [e dispatchThreads:MTLSizeMake(n_active,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // (f) Optional ext_accel backward — pair_forces and/or springs.
            //     Both feed ∂L/∂(ext_accel) = dt · ∂L/∂(v_after_apply).
            //
            //     CRITICAL: must run after the previous cmd buffer
            //     completes — the host-side multiply below reads
            //     bGv_old_new which is written by predict_bw above.
            if (use_pair || use_springs) {
                // host: bGext = dt · bGv_old_new (= ∂L/∂ext_accel)
                float *gnew = (float *)[bGv_old_new contents];
                float *gext = (float *)[bGext contents];
                for (uint32_t i = 0; i < 3 * n_active; i++) {
                    gext[i] = dt * gnew[i];
                }

                id<MTLCommandBuffer> cmd2 = [ctx.queue commandBuffer];

                if (use_pair) {
                    // restore pair_density used by THIS step's forward pair_forces
                    memcpy([bD contents], step_state + 11 * n_active, s_b);

                    id<MTLComputeCommandEncoder> e = [cmd2 computeCommandEncoder];
                    [e setComputePipelineState:pso_pair_bw];
                    [e setBuffer:bX           offset:0 atIndex:0];
                    [e setBuffer:bV           offset:0 atIndex:1];
                    [e setBuffer:bSortedS     offset:0 atIndex:2];
                    [e setBuffer:bCellStart   offset:0 atIndex:3];
                    [e setBuffer:bD           offset:0 atIndex:4];
                    [e setBuffer:bGext        offset:0 atIndex:5];
                    [e setBuffer:bGx_old_new  offset:0 atIndex:6];
                    [e setBuffer:bGv_old_new  offset:0 atIndex:7];
                    [e setBuffer:bH           offset:0 atIndex:8];
                    [e setBuffer:bH2          offset:0 atIndex:9];
                    [e setBuffer:bM           offset:0 atIndex:10];
                    [e setBuffer:bSS          offset:0 atIndex:11];
                    [e setBuffer:bViscPair    offset:0 atIndex:12];
                    [e setBuffer:bViscAmp     offset:0 atIndex:13];
                    [e setBuffer:bSurfAmp     offset:0 atIndex:14];
                    [e setBuffer:bNa          offset:0 atIndex:15];
                    [e setBuffer:bGridDim     offset:0 atIndex:16];
                    [e setBuffer:bGridOrigin  offset:0 atIndex:17];
                    [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [e endEncoding];
                }

                if (use_springs) {
                    id<MTLComputeCommandEncoder> e = [cmd2 computeCommandEncoder];
                    [e setComputePipelineState:pso_spring_bw];
                    [e setBuffer:bX            offset:0 atIndex:0];
                    [e setBuffer:bBondIJ_sp    offset:0 atIndex:1];
                    [e setBuffer:bBondRest_sp  offset:0 atIndex:2];
                    [e setBuffer:bGext         offset:0 atIndex:3];
                    [e setBuffer:bGx_old_new   offset:0 atIndex:4];
                    [e setBuffer:bSpringK      offset:0 atIndex:5];
                    [e setBuffer:bNbonds_sp    offset:0 atIndex:6];
                    [e setBuffer:bNa           offset:0 atIndex:7];
                    [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [e endEncoding];
                }

                [cmd2 commit]; [cmd2 waitUntilCompleted];

                // (g) Param gradients — analytic ∂L/∂(spring_K) and
                //     ∂L/∂(visc_pair_coef) from this step's contribution.
                //     Both kernels read bGext (= dt · ∂L/∂v_after_apply)
                //     and produce per-particle scalar partials; host sums.
                id<MTLCommandBuffer> cmd3 = [ctx.queue commandBuffer];
                if (use_springs) {
                    id<MTLComputeCommandEncoder> e = [cmd3 computeCommandEncoder];
                    [e setComputePipelineState:pso_spring_K_part];
                    [e setBuffer:bX             offset:0 atIndex:0];
                    [e setBuffer:bBondIJ_sp     offset:0 atIndex:1];
                    [e setBuffer:bBondRest_sp   offset:0 atIndex:2];
                    [e setBuffer:bGext          offset:0 atIndex:3];
                    [e setBuffer:bSpringPart    offset:0 atIndex:4];
                    [e setBuffer:bNbonds_sp     offset:0 atIndex:5];
                    [e setBuffer:bNa            offset:0 atIndex:6];
                    [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [e endEncoding];
                }
                if (use_pair) {
                    id<MTLComputeCommandEncoder> e = [cmd3 computeCommandEncoder];
                    [e setComputePipelineState:pso_visc_K_part];
                    [e setBuffer:bX             offset:0 atIndex:0];
                    [e setBuffer:bV             offset:0 atIndex:1];
                    [e setBuffer:bSortedS       offset:0 atIndex:2];
                    [e setBuffer:bCellStart     offset:0 atIndex:3];
                    [e setBuffer:bD             offset:0 atIndex:4];  // pair_density (already restored)
                    [e setBuffer:bGext          offset:0 atIndex:5];
                    [e setBuffer:bViscPart      offset:0 atIndex:6];
                    [e setBuffer:bH             offset:0 atIndex:7];
                    [e setBuffer:bH2            offset:0 atIndex:8];
                    [e setBuffer:bSS            offset:0 atIndex:9];
                    [e setBuffer:bViscAmp       offset:0 atIndex:10];
                    [e setBuffer:bNa            offset:0 atIndex:11];
                    [e setBuffer:bGridDim       offset:0 atIndex:12];
                    [e setBuffer:bGridOrigin    offset:0 atIndex:13];
                    [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [e endEncoding];
                }
                [cmd3 commit]; [cmd3 waitUntilCompleted];

                // Host-side sum of per-particle partials.
                if (use_springs) {
                    float *sp = (float *)[bSpringPart contents];
                    float s = 0.0f;
                    for (uint32_t i = 0; i < n_active; i++) s += sp[i];
                    total_grad_spring_K += s;
                }
                if (use_pair) {
                    float *vp = (float *)[bViscPart contents];
                    float s = 0.0f;
                    for (uint32_t i = 0; i < n_active; i++) s += vp[i];
                    total_grad_visc_K += s;
                }
            }

            // Sum ρ gradient (kernel partial + implicit via grad_C).
            //   implicit = -(grad_grad_C · grad_C) / ρ_rest
            //   Note: bGgC currently holds ∂L/∂grad_C from solve_dens_bwd; bGc holds grad_C (forward saved).
            float kernel_rho = 0.0f;
            float *gr = (float *)[bGrho_per contents];
            for (uint32_t i = 0; i < n_active; i++) kernel_rho += gr[i];
            float implicit_rho = 0.0f;
            float *ggC = (float *)[bGgC contents];
            float *gC_fwd = (float *)[bGc contents];
            for (uint32_t i = 0; i < n_active; i++) {
                float dot = 0.0f;
                for (int ax = 0; ax < 3; ax++)
                    dot += ggC[i * 3 + ax] * gC_fwd[i * 3 + ax];
                implicit_rho -= dot / rho_rest;
            }
            total_grad_rho += kernel_rho + implicit_rho;

            // Sum α gradient (per-particle ∂L/∂A from density solve backward).
            float *ga = (float *)[bGalpha_per contents];
            for (uint32_t i = 0; i < n_active; i++) total_grad_alpha_inv_dt2 += ga[i];

            // Per-step gradient clipping (when enabled via BWD_CLIP_NORM).
            // Apply to bGx_old_new and bGv_old_new BEFORE they become
            // running gradients for the previous chain step. NaN/Inf
            // values are zeroed so they don't propagate forever.
            if (clip_env) {
                float *gx = (float *)[bGx_old_new contents];
                float *gv = (float *)[bGv_old_new contents];
                size_t n_floats = 3 * (size_t)n_active;
                // Zero out NaN/Inf elements first (degenerate-contact
                // particles can produce them even after kernel epsilons).
                int n_bad_x = 0, n_bad_v = 0;
                for (size_t t = 0; t < n_floats; t++) {
                    if (!isfinite(gx[t])) { gx[t] = 0.0f; n_bad_x++; }
                    if (!isfinite(gv[t])) { gv[t] = 0.0f; n_bad_v++; }
                }
                // L2-norm clip.
                double sum_x = 0.0, sum_v = 0.0;
                for (size_t t = 0; t < n_floats; t++) {
                    sum_x += (double)gx[t] * gx[t];
                    sum_v += (double)gv[t] * gv[t];
                }
                float nx = (float)sqrt(sum_x);
                float nv = (float)sqrt(sum_v);
                if (nx > clip_norm_x && nx > 0.0f) {
                    float s = clip_norm_x / nx;
                    for (size_t t = 0; t < n_floats; t++) gx[t] *= s;
                }
                if (nv > clip_norm_v && nv > 0.0f) {
                    float s = clip_norm_v / nv;
                    for (size_t t = 0; t < n_floats; t++) gv[t] *= s;
                }
                if ((n_bad_x || n_bad_v) && k % 50 == 0) {
                    fprintf(stderr, "[xpbd_full_bwd k=%d] zeroed NaN/Inf: "
                            "%d in grad_x, %d in grad_v\n",
                            (int)k, n_bad_x, n_bad_v);
                }
            }

            // Promote per-step grads to running for previous step.
            memcpy([bGx_running contents], [bGx_old_new contents], pos_b);
            memcpy([bGv_running contents], [bGv_old_new contents], pos_b);
        }

        FILE *o1 = fopen(argv[14], "wb"); fwrite([bGx_running contents], 1, pos_b, o1); fclose(o1);
        FILE *o2 = fopen(argv[15], "wb"); fwrite([bGv_running contents], 1, pos_b, o2); fclose(o2);
        FILE *o3 = fopen(argv[16], "wb"); fwrite(&total_grad_rho, 1, sizeof(float), o3); fclose(o3);
        // Optional param-grad outputs (positions 22, 23 in argv).
        if (argc >= 22) {
            FILE *o4 = fopen(argv[21], "wb");
            fwrite(&total_grad_spring_K, 1, sizeof(float), o4); fclose(o4);
        }
        if (argc >= 23) {
            FILE *o5 = fopen(argv[22], "wb");
            fwrite(&total_grad_visc_K, 1, sizeof(float), o5); fclose(o5);
        }
        if (argc >= 25) {
            // ∂L/∂alpha_dens = ∂L/∂(alpha_inv_dt2) · ∂(alpha_inv_dt2)/∂alpha_dens
            //                = total_grad_alpha_inv_dt2 · (1/dt²)
            float total_grad_alpha_dens = total_grad_alpha_inv_dt2 / (dt * dt);
            FILE *o6 = fopen(argv[24], "wb");
            fwrite(&total_grad_alpha_dens, 1, sizeof(float), o6); fclose(o6);
        }
    }
    free(all); free(pos_static); free(grad_x_fin);
    if (mem_tris) free(mem_tris);
    if (mem_pidx) free(mem_pidx);
    if (sorted_static_buf) free(sorted_static_buf);
    if (cell_start_buf)    free(cell_start_buf);
    if (bond_ij_data) free(bond_ij_data);
    if (bond_rest_data) free(bond_rest_data);
    return 0;
}

