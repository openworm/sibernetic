// ops_xpbd_step.mm — M7 imperative XPBD pipeline (run_xpbd_step).
//
// Single-call op that runs N XPBD steps inside one process: predict +
// (density solve, distance constraint or springs, pair forces, floor)
// × n_iters per step + update_velocity. Forward-only; for the
// differentiable pipeline see ops_xpbd_full.mm.

#include "metal_common.h"

// ──────────────────────────────────────────────────────────────────────
// M7 — xpbd_step: one full XPBD timestep.
//
// Pipeline per step:
//   1. predict_positions (apply gravity to vel, integrate to x_pred)
//   2. for iter in 0..n_iters:
//        a. dist_active_active   (recompute distances after each
//        b. dist_active_static    projection — particles moved)
//        c. wpoly6_inplace × 2    (kernel-evaluate r²_aa and r²_as)
//        d. rowsum_density × 2   (compute ρ_aa and ρ_as)
//        e. add_inplace          (ρ = ρ_aa + ρ_as)
//        f. density_constraint_grad
//        g. solve_density_constraint
//        h. solve_floor_constraint
//   3. update_velocities (recover v from position change)
//
// Inputs (binary fp32 files, little-endian):
//   pos_active.bin: [n_active*3] starting positions of dynamic particles
//   vel_active.bin: [n_active*3] starting velocities
//   pos_static.bin: [n_static*3] frozen boundary positions
// Outputs:
//   pos_active_out.bin, vel_active_out.bin
//
// Lagrange multipliers (λ) are reset to zero at the start of each
// step — XPBD's "warm restart" form rather than persistent λ. This
// matches Macklin 2013 PBD-Fluids.
// ──────────────────────────────────────────────────────────────────────
int run_xpbd_step(int argc, char **argv) {
    // Required: op + 14 base + 3 bonds = 18 args minimum
    // Optional: + bench_steps = 19, + sim_scale = 20, + visc_pair_coef = 21
    // Pass n_bonds=0 to skip distance constraints; in that case
    // bonds.bin and alpha_dist are read but ignored.
    // sim_scale is the Sibernetic unit-system bridge: 1 particle unit
    // = sim_scale meters. Default 1.0 means positions and velocity share
    // a unit system (toy-test convention). For Sibernetic configs use
    // sim_scale = 7.4e-6 (≈ Sibernetic's `simulationScale` for the
    // demo1 mass).
    // visc_pair_coef enables the Sibernetic-equivalent viscosity +
    // surface tension pair-force pass. Default 0 disables it; 1e-4 is
    // Sibernetic's main path coefficient (see sphFluid.cl:602-624).
    if (argc < 18 || argc > 42) {
        fprintf(stderr,
            "usage: sib_metal xpbd_step "
            "<n_active> <n_static> <h> <mass> <rho_rest> <dt> <gravity_y> "
            "<floor_y> <alpha_density> <n_iters> "
            "<pos_active.bin> <vel_active.bin> <pos_static.bin> "
            "<n_bonds> <bonds.bin> <alpha_dist> [bench_steps] [sim_scale] "
            "[visc_pair_coef] [spring_K] [restitution] "
            "[n_membranes] [n_elastic] [r0] [membranes.bin] [pmem_index.bin] "
            "[n_anchors] [anchors.bin] [pressure_k]\n"
            "       (outputs written to /tmp/xpbd_{pos,vel}_out.bin)\n"
            "       bonds.bin format: per bond, [int32 i, int32 j, "
            "float32 rest_len, float32 _pad] (16 bytes each)\n"
            "       spring_K > 0 enables Hooke spring bonds AND disables\n"
            "       the rigid XPBD distance constraint (Sibernetic mode).\n"
            "       spring_K = elasticityCoefficient * sim_scale\n"
            "       (Sibernetic default: 3e8 * 7.4e-6 = 2220).\n"
            "       restitution ∈ [0,1]: floor elasticity (default 0).\n"
            "       n_membranes > 0 enables membrane interaction (M10):\n"
            "         membranes.bin: [n_membranes, 3] int32 vertex IDs\n"
            "                        in active-buffer local indexing\n"
            "         pmem_index.bin: [n_elastic*7] int32 incident-tri\n"
            "                         indices per elastic, -1 sentinel\n"
            "         r0: membrane neighbor radius (Sibernetic uses h/2)\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h           = (float)atof(argv[4]);
    float mass        = (float)atof(argv[5]);
    float rho_rest    = (float)atof(argv[6]);
    float dt          = (float)atof(argv[7]);
    float gravity_y   = (float)atof(argv[8]);
    float floor_y     = (float)atof(argv[9]);
    float alpha_dens  = (float)atof(argv[10]);
    uint32_t n_iters  = (uint32_t)atoi(argv[11]);
    const char *path_pos_active  = argv[12];
    const char *path_vel_active  = argv[13];
    const char *path_pos_static  = argv[14];
    uint32_t n_bonds   = (uint32_t)atoi(argv[15]);
    const char *path_bonds       = argv[16];
    float alpha_dist  = (float)atof(argv[17]);
    int bench_steps = (argc >= 19) ? atoi(argv[18]) : 1;
    if (bench_steps < 1) bench_steps = 1;
    float sim_scale = (argc >= 20) ? (float)atof(argv[19]) : 1.0f;
    float sim_scale_inv = 1.0f / sim_scale;
    float visc_pair_coef = (argc >= 21) ? (float)atof(argv[20]) : 0.0f;
    bool use_pair_forces = visc_pair_coef != 0.0f;
    float spring_K       = (argc >= 22) ? (float)atof(argv[21]) : 0.0f;
    bool use_springs     = spring_K != 0.0f;
    float restitution    = (argc >= 23) ? (float)atof(argv[22]) : 0.0f;
    // M10 membrane interaction (off when n_membranes==0).
    uint32_t n_membranes = (argc >= 24) ? (uint32_t)atoi(argv[23]) : 0u;
    uint32_t n_elastic   = (argc >= 25) ? (uint32_t)atoi(argv[24]) : 0u;
    float r0_mem         = (argc >= 26) ? (float)atof(argv[25]) : 0.0f;
    const char *path_membranes = (argc >= 27) ? argv[26] : NULL;
    const char *path_pmem_idx  = (argc >= 28) ? argv[27] : NULL;
    bool use_membranes = (n_membranes > 0) && (n_elastic > 0)
                         && path_membranes && path_pmem_idx
                         && (r0_mem > 0.0f);
    // M11 boundary anchor springs (off when n_anchors == 0).
    uint32_t n_anchors = (argc >= 29) ? (uint32_t)atoi(argv[28]) : 0u;
    const char *path_anchors = (argc >= 30) ? argv[29] : NULL;
    bool use_anchors = (n_anchors > 0) && path_anchors;
    // M12 SPH pressure force (counter to surface-tension cohesion).
    float pressure_k = (argc >= 31) ? (float)atof(argv[30]) : 0.0f;
    bool use_pressure = pressure_k > 0.0f;
    // Anchor stiffness — separate from bond stiffness so sheets can flex
    // at edges while internal bonds stay rigid. 0 means use spring_K.
    float anchor_k_in = (argc >= 32) ? (float)atof(argv[31]) : 0.0f;
    float anchor_k = (anchor_k_in > 0.0f) ? anchor_k_in : spring_K;
    // Velocity clamp magnitude (m/s) — caps |v| to prevent CFL-violating
    // teleport past boundary walls. 0 = disabled.
    float vel_clamp = (argc >= 33) ? (float)atof(argv[32]) : 0.0f;
    bool use_vel_clamp = vel_clamp > 0.0f;
    // Box clamp: 6 args xmin xmax ymin ymax zmin zmax. Activates when
    // xmax > xmin (sentinel = all zeros = disabled).
    float box_min[3] = {0,0,0}, box_max[3] = {0,0,0};
    bool use_box_clamp = false;
    if (argc >= 39) {
        box_min[0] = (float)atof(argv[33]);
        box_max[0] = (float)atof(argv[34]);
        box_min[1] = (float)atof(argv[35]);
        box_max[1] = (float)atof(argv[36]);
        box_min[2] = (float)atof(argv[37]);
        box_max[2] = (float)atof(argv[38]);
        use_box_clamp = box_max[0] > box_min[0];
    }
    // Active muscle drive: anchor rest length sinusoidal modulation.
    // muscle_amp = 0 (default) = pure passive Hooke spring.
    float muscle_freq = (argc >= 40) ? (float)atof(argv[39]) : 0.0f;
    float muscle_amp  = (argc >= 41) ? (float)atof(argv[40]) : 0.0f;
    // Time offset (seconds) — added to step*dt so muscle phase persists
    // across chunked invocations. dump_metal_trajectory.py tracks
    // cumulative steps and passes the offset.
    float time_offset_s = (argc >= 42) ? (float)atof(argv[41]) : 0.0f;
    // Anchors are only meaningful when springs are on (uses spring_K).
    bool need_ext_accel_with_anchors = (use_anchors && (spring_K != 0.0f));
    // Springs replace the rigid XPBD distance constraint when on.
    bool use_xpbd_bonds  = (n_bonds > 0) && !use_springs;
    // Springs feed into the same ext_accel buffer as pair_forces, so
    // we need the apply_ext_accel scaffolding even if pair_forces is off.
    bool need_ext_accel  = use_pair_forces || use_springs || need_ext_accel_with_anchors || use_pressure;
    // Sibernetic-equivalent precomputed amps (see owPhysicsConstant.h).
    // h_scaled = h * sim_scale; divgradWviscoCoeff = 45/(π·h_s^6);
    // surfTensCoeff = mass·Wpoly6Coef·sim_scale; Wpoly6Coef = 315/(64π·h_s^9).
    //
    // Computed in fp64 because intermediates underflow fp32 — e.g.
    // h_s^9 ≈ 7.5e-42 is subnormal in fp32, making 1/h_s^9 inf.
    // The final surf_amp itself fits fp32 fine (~5e27).
    double h_scaled  = (double)h * (double)sim_scale;
    double h_s6      = pow(h_scaled, 6.0);
    double h_s9      = pow(h_scaled, 9.0);
    double divgradWvisco = 45.0 / (M_PI * h_s6);
    // Sibernetic divides viscosity by density-in-kg/m³, but our density
    // is in kg/unit³ — a factor of sim_scale^3 smaller. Compensate
    // by pre-multiplying the amp so the in-kernel division by ρ
    // (in our units) gives the right Sibernetic-equivalent acceleration.
    double visc_amp_d = 1.5 * (double)mass * divgradWvisco
                              * pow((double)sim_scale, 3.0);
    double wpoly6_si  = 315.0 / (64.0 * M_PI * h_s9);
    double surfTensCoef_si = (double)mass * wpoly6_si * (double)sim_scale;
    double surf_amp_d = -1.7e-9 * surfTensCoef_si / (double)mass;
    float visc_amp = (float)visc_amp_d;
    float surf_amp = (float)surf_amp_d;

    const char *out_pos_path = "/tmp/xpbd_pos_out.bin";
    const char *out_vel_path = "/tmp/xpbd_vel_out.bin";

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));
    float alpha_inv_dt2 = alpha_dens / (dt * dt);
    float alpha_dist_inv_dt2 = alpha_dist / (dt * dt);
    float mass_inv = 1.0f / mass;

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        if (fread(buf, sizeof(float), n_floats, f) != n_floats) {
            fprintf(stderr, "short read on %s\n", path); exit(1);
        }
        fclose(f);
        return buf;
    };

    float *pos_active_init = read_floats(path_pos_active, (size_t)n_active * 3);
    float *vel_active_init = read_floats(path_vel_active, (size_t)n_active * 3);
    float *pos_static      = read_floats(path_pos_static, (size_t)n_static * 3);

    // Build static spatial grid (one-time per process invocation).
    float *sorted_static = NULL;
    int *cell_start = NULL;
    GridDim3 grid_dim_struct;
    GridOrigin3 grid_origin_struct;
    int n_cells = 0;
    build_static_grid(pos_static, n_static, h,
                      &sorted_static, &cell_start,
                      &grid_dim_struct, &grid_origin_struct, &n_cells);

    // Load bonds. Format per bond: [i:int32, j:int32, rest_len:float32, pad:float32]
    // Total 16 bytes per bond. Read raw and unpack.
    void *bonds_raw = NULL;
    if (n_bonds > 0) {
        FILE *f = fopen(path_bonds, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path_bonds); exit(1); }
        bonds_raw = malloc((size_t)n_bonds * 16);
        if (fread(bonds_raw, 16, n_bonds, f) != n_bonds) {
            fprintf(stderr, "short read on %s\n", path_bonds); exit(1);
        }
        fclose(f);
    }

    // Load anchor data if enabled (8 floats per anchor, see shaders.metal).
    float *anchors_data = NULL;
    if (use_anchors) {
        size_t anchor_bytes = (size_t)n_anchors * 8 * sizeof(float);
        FILE *fa = fopen(path_anchors, "rb");
        if (!fa) { fprintf(stderr, "cannot open %s\n", path_anchors); exit(1); }
        anchors_data = (float *)malloc(anchor_bytes);
        if (fread(anchors_data, 1, anchor_bytes, fa) != anchor_bytes) {
            fprintf(stderr, "short read on %s\n", path_anchors); exit(1);
        }
        fclose(fa);
    }

    // Load membrane topology if enabled.
    int32_t *mem_tris = NULL;
    int32_t *mem_pidx = NULL;
    if (use_membranes) {
        size_t tris_bytes = (size_t)n_membranes * 3 * sizeof(int32_t);
        size_t pidx_bytes = (size_t)n_elastic * 7 * sizeof(int32_t);
        FILE *f1 = fopen(path_membranes, "rb");
        if (!f1) { fprintf(stderr, "cannot open %s\n", path_membranes); exit(1); }
        mem_tris = (int32_t *)malloc(tris_bytes);
        if (fread(mem_tris, 1, tris_bytes, f1) != tris_bytes) {
            fprintf(stderr, "short read on %s\n", path_membranes); exit(1);
        }
        fclose(f1);
        FILE *f2 = fopen(path_pmem_idx, "rb");
        if (!f2) { fprintf(stderr, "cannot open %s\n", path_pmem_idx); exit(1); }
        mem_pidx = (int32_t *)malloc(pidx_bytes);
        if (fread(mem_pidx, 1, pidx_bytes, f2) != pidx_bytes) {
            fprintf(stderr, "short read on %s\n", path_pmem_idx); exit(1);
        }
        fclose(f2);
    }

    @autoreleasepool {
        MetalCtx ctx = make_ctx();

        // Compile all the PSOs we need (this is the expensive Metal
        // setup; per-step only command-buffer creation/dispatch).
        id<MTLComputePipelineState> pso_dist_aa  = make_pso(ctx, "dist_active_active");
        id<MTLComputePipelineState> pso_dist_as  = make_pso(ctx, "dist_active_static");
        id<MTLComputePipelineState> pso_wpoly6   = make_pso(ctx, "wpoly6_inplace");
        id<MTLComputePipelineState> pso_density  = make_pso(ctx, "rowsum_density");
        id<MTLComputePipelineState> pso_wp_rs_fused = make_pso(ctx, "wpoly6_rowsum_density_fused");
        id<MTLComputePipelineState> pso_d_grad_combined = make_pso(ctx, "density_grad_combined");
        id<MTLComputePipelineState> pso_d_grad_mega = make_pso(ctx, "density_grad_mega_fused");
        id<MTLComputePipelineState> pso_d_grad_grid = make_pso(ctx, "density_grad_mega_grid");
        (void)pso_wpoly6; (void)pso_density; (void)pso_wp_rs_fused;
        (void)pso_d_grad_combined; (void)pso_dist_aa; (void)pso_dist_as;
        (void)pso_d_grad_mega;
        // pso_d_grad_grid is THE inner-loop kernel: spatial-grid neighbor
        // search for static, plus dense for active-active.
        id<MTLComputePipelineState> pso_addin    = make_pso(ctx, "add_inplace");
        id<MTLComputePipelineState> pso_grad_C   = make_pso(ctx, "density_constraint_grad");
        id<MTLComputePipelineState> pso_predict  = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_solve_d  = make_pso(ctx, "solve_density_constraint");
        id<MTLComputePipelineState> pso_solve_f  = make_pso(ctx, "solve_floor_constraint");
        id<MTLComputePipelineState> pso_updvel   = make_pso(ctx, "update_velocities");
        id<MTLComputePipelineState> pso_solve_b  = use_xpbd_bonds
            ? make_pso(ctx, "solve_distance_constraints_seq") : nil;
        id<MTLComputePipelineState> pso_pair_forces = use_pair_forces
            ? make_pso(ctx, "pair_forces_grid") : nil;
        id<MTLComputePipelineState> pso_apply_ext   = need_ext_accel
            ? make_pso(ctx, "apply_ext_accel") : nil;
        id<MTLComputePipelineState> pso_spring      = use_springs
            ? make_pso(ctx, "spring_bonds_force") : nil;
        id<MTLComputePipelineState> pso_anchor      = (use_anchors && use_springs)
            ? make_pso(ctx, "spring_anchor_force") : nil;
        id<MTLComputePipelineState> pso_pressure    = use_pressure
            ? make_pso(ctx, "pressure_force_grid") : nil;
        id<MTLComputePipelineState> pso_clampv      = use_vel_clamp
            ? make_pso(ctx, "clamp_velocity") : nil;
        id<MTLComputePipelineState> pso_clampbox    = use_box_clamp
            ? make_pso(ctx, "clamp_to_box") : nil;
        // M10 membrane PSOs (compiled only when use_membranes).
        id<MTLComputePipelineState> pso_mem_clear   = use_membranes
            ? make_pso(ctx, "clear_membrane_correction") : nil;
        id<MTLComputePipelineState> pso_mem_acc     = use_membranes
            ? make_pso(ctx, "accumulate_membrane_correction") : nil;
        id<MTLComputePipelineState> pso_mem_apply   = use_membranes
            ? make_pso(ctx, "apply_membrane_correction") : nil;

        size_t pos_a_bytes = (size_t)n_active * 3 * sizeof(float);
        size_t pos_s_bytes = (size_t)n_static * 3 * sizeof(float);
        size_t r2_aa_bytes = (size_t)n_active * n_active * sizeof(float);
        size_t r2_as_bytes = (size_t)n_active * n_static * sizeof(float);
        size_t dens_bytes  = (size_t)n_active * sizeof(float);

        // Persistent buffers across all steps + iters.
        id<MTLBuffer> bufPosOld   = [ctx.device newBufferWithBytes:pos_active_init
            length:pos_a_bytes options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufPosPred  = [ctx.device newBufferWithLength:pos_a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufVel      = [ctx.device newBufferWithBytes:vel_active_init
            length:pos_a_bytes options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufPosStat  = [ctx.device newBufferWithBytes:pos_static
            length:pos_s_bytes options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufSortedStatic = (n_static > 0)
            ? [ctx.device newBufferWithBytes:sorted_static
                length:pos_s_bytes options:MTLResourceStorageModeShared]
            : [ctx.device newBufferWithLength:16 options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufCellStart = [ctx.device newBufferWithBytes:cell_start
            length:(size_t)(n_cells + 1) * sizeof(int)
            options:MTLResourceStorageModeShared];
        // Pack grid_dim as int3 (4 ints in metal alignment) and grid_origin as packed_float3.
        int grid_dim_packed[4] = {grid_dim_struct.x, grid_dim_struct.y, grid_dim_struct.z, 0};
        float grid_origin_packed[3] = {grid_origin_struct.x, grid_origin_struct.y, grid_origin_struct.z};
        id<MTLBuffer> bufGridDim = [ctx.device newBufferWithBytes:grid_dim_packed
            length:sizeof(int) * 4 options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufGridOrigin = [ctx.device newBufferWithBytes:grid_origin_packed
            length:sizeof(float) * 3 options:MTLResourceStorageModeShared];
        (void)bufPosStat; // kept for backward compat with other kernels
        id<MTLBuffer> bufR2aa     = [ctx.device newBufferWithLength:r2_aa_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufR2as     = [ctx.device newBufferWithLength:r2_as_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufWaa      = [ctx.device newBufferWithLength:r2_aa_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufWas      = [ctx.device newBufferWithLength:r2_as_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufDensAa   = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufDens     = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufGradC    = [ctx.device newBufferWithLength:pos_a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufDenomH   = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufLambda   = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];
        // Per-particle external acceleration (m/s²) from viscosity +
        // surface tension. Recomputed each step in pair_forces_grid.
        id<MTLBuffer> bufExtAccel = need_ext_accel
            ? [ctx.device newBufferWithLength:pos_a_bytes
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufSpringK = use_springs
            ? [ctx.device newBufferWithBytes:&spring_K length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        // Anchor buffers (only when use_anchors AND use_springs).
        id<MTLBuffer> bufAnchors = (use_anchors && use_springs)
            ? [ctx.device newBufferWithBytes:anchors_data
                  length:(size_t)n_anchors * 8 * sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufNanchors = (use_anchors && use_springs)
            ? [ctx.device newBufferWithBytes:&n_anchors length:sizeof(uint32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufPressureK = use_pressure
            ? [ctx.device newBufferWithBytes:&pressure_k length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufAnchorK = (use_anchors && use_springs)
            ? [ctx.device newBufferWithBytes:&anchor_k length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufVelClamp = use_vel_clamp
            ? [ctx.device newBufferWithBytes:&vel_clamp length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufBoxMin = use_box_clamp
            ? [ctx.device newBufferWithBytes:box_min length:3*sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufBoxMax = use_box_clamp
            ? [ctx.device newBufferWithBytes:box_max length:3*sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        // Muscle drive: scalars + per-step time_t buffer (updated each step).
        id<MTLBuffer> bufMuscleFreq = (use_anchors && use_springs)
            ? [ctx.device newBufferWithBytes:&muscle_freq length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufMuscleAmp = (use_anchors && use_springs)
            ? [ctx.device newBufferWithBytes:&muscle_amp length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        float time_t_init = 0.0f;
        id<MTLBuffer> bufTimeT = (use_anchors && use_springs)
            ? [ctx.device newBufferWithBytes:&time_t_init length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;
        // M10 membrane buffers.
        id<MTLBuffer> bufMembranes = use_membranes
            ? [ctx.device newBufferWithBytes:mem_tris
                  length:(size_t)n_membranes * 3 * sizeof(int32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufPmemIdx   = use_membranes
            ? [ctx.device newBufferWithBytes:mem_pidx
                  length:(size_t)n_elastic * 7 * sizeof(int32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufMemCorr   = use_membranes
            ? [ctx.device newBufferWithLength:pos_a_bytes
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufNelastic  = use_membranes
            ? [ctx.device newBufferWithBytes:&n_elastic
                  length:sizeof(uint32_t)
                  options:MTLResourceStorageModeShared] : nil;
        id<MTLBuffer> bufR0        = use_membranes
            ? [ctx.device newBufferWithBytes:&r0_mem length:sizeof(float)
                  options:MTLResourceStorageModeShared] : nil;

        // Bond buffers: layout in memory is [int32 i, int32 j, float32 rest, float32 pad]
        // We split into bond_ij (int2 per bond) and rest_len (float per bond)
        // for clean kernel access. Done by repacking bonds_raw on host.
        id<MTLBuffer> bufBondIJ   = nil;
        id<MTLBuffer> bufBondRest = nil;
        id<MTLBuffer> bufBondLam  = nil;
        id<MTLBuffer> bufNbonds   = [ctx.device newBufferWithBytes:&n_bonds
            length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufAdistDt2 = [ctx.device newBufferWithBytes:&alpha_dist_inv_dt2
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufMassInv  = [ctx.device newBufferWithBytes:&mass_inv
            length:sizeof(float) options:MTLResourceStorageModeShared];
        if (n_bonds > 0) {
            int32_t *bond_ij   = (int32_t *)malloc((size_t)n_bonds * 2 * sizeof(int32_t));
            float   *bond_rest = (float *)malloc((size_t)n_bonds * sizeof(float));
            uint8_t *raw       = (uint8_t *)bonds_raw;
            for (uint32_t b = 0; b < n_bonds; b++) {
                memcpy(&bond_ij[b * 2], raw + b * 16, 8);  // i, j
                memcpy(&bond_rest[b],   raw + b * 16 + 8, 4); // rest_len
                // raw + b*16 + 12 is padding, ignored
            }
            bufBondIJ   = [ctx.device newBufferWithBytes:bond_ij
                length:(size_t)n_bonds * 2 * sizeof(int32_t)
                options:MTLResourceStorageModeShared];
            bufBondRest = [ctx.device newBufferWithBytes:bond_rest
                length:(size_t)n_bonds * sizeof(float)
                options:MTLResourceStorageModeShared];
            bufBondLam  = [ctx.device newBufferWithLength:(size_t)n_bonds * sizeof(float)
                options:MTLResourceStorageModeShared];
            free(bond_ij); free(bond_rest);
        }

        // Constants — wrapped in shared buffers so kernels can read them.
        id<MTLBuffer> bufNa  = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNs  = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufH   = [ctx.device newBufferWithBytes:&h       length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufH2  = [ctx.device newBufferWithBytes:&h2      length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufPoly6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufSpiky = [ctx.device newBufferWithBytes:&spiky_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufMass  = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufRho   = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufDt    = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufG     = [ctx.device newBufferWithBytes:&gravity_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufRestitution = [ctx.device newBufferWithBytes:&restitution length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufFloor = [ctx.device newBufferWithBytes:&floor_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufAdt2  = [ctx.device newBufferWithBytes:&alpha_inv_dt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufSimScaleInv = [ctx.device newBufferWithBytes:&sim_scale_inv
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufSimScale    = [ctx.device newBufferWithBytes:&sim_scale
            length:sizeof(float) options:MTLResourceStorageModeShared];
        // Pair-force scalars (only used when use_pair_forces is true).
        id<MTLBuffer> bufViscPair    = [ctx.device newBufferWithBytes:&visc_pair_coef
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufViscAmp     = [ctx.device newBufferWithBytes:&visc_amp
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufSurfAmp     = [ctx.device newBufferWithBytes:&surf_amp
            length:sizeof(float) options:MTLResourceStorageModeShared];
        // n_total for wpoly6_inplace (different per call).
        uint32_t n_aa_total = n_active * n_active;
        uint32_t n_as_total = n_active * n_static;
        id<MTLBuffer> bufNaaTot = [ctx.device newBufferWithBytes:&n_aa_total length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNasTot = [ctx.device newBufferWithBytes:&n_as_total length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        // Pair-force step needs a density value at start of step. For
        // step 0 we have nothing yet, so seed bufDens with rho_rest as
        // a one-time approximation. Subsequent steps inherit density
        // from the previous step's last inner-XPBD iteration.
        if (use_pair_forces) {
            float *dens0 = (float *)[bufDens contents];
            for (uint32_t i = 0; i < n_active; i++) dens0[i] = rho_rest;
        }

        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int step = 0; step < bench_steps; step++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];

            // ── 0. (optional) External forces (visc + surf + springs)
            //    accumulate into ext_accel, then apply to velocity ──
            if (need_ext_accel) {
                // Zero the ext_accel buffer at start of step (kernels accumulate).
                memset([bufExtAccel contents], 0, pos_a_bytes);
            }
            if (use_pair_forces) {
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_pair_forces];
                [enc setBuffer:bufPosOld       offset:0 atIndex:0];
                [enc setBuffer:bufVel          offset:0 atIndex:1];
                [enc setBuffer:bufSortedStatic offset:0 atIndex:2];
                [enc setBuffer:bufCellStart    offset:0 atIndex:3];
                [enc setBuffer:bufDens         offset:0 atIndex:4];
                [enc setBuffer:bufExtAccel     offset:0 atIndex:5];
                [enc setBuffer:bufH            offset:0 atIndex:6];
                [enc setBuffer:bufH2           offset:0 atIndex:7];
                [enc setBuffer:bufMass         offset:0 atIndex:8];
                [enc setBuffer:bufSimScale     offset:0 atIndex:9];
                [enc setBuffer:bufViscPair     offset:0 atIndex:10];
                [enc setBuffer:bufViscAmp      offset:0 atIndex:11];
                [enc setBuffer:bufSurfAmp      offset:0 atIndex:12];
                [enc setBuffer:bufNa           offset:0 atIndex:13];
                [enc setBuffer:bufGridDim      offset:0 atIndex:14];
                [enc setBuffer:bufGridOrigin   offset:0 atIndex:15];
                [enc dispatchThreads:MTLSizeMake(32, n_active, 1)
                    threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
                [enc endEncoding];
            }
            if (use_springs) {
                // Hooke spring bond forces, accumulate into ext_accel.
                // (Replaces XPBD distance constraint for Sibernetic-equivalent
                // visco-elastic behavior.)
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_spring];
                [enc setBuffer:bufPosOld   offset:0 atIndex:0];
                [enc setBuffer:bufBondIJ   offset:0 atIndex:1];
                [enc setBuffer:bufBondRest offset:0 atIndex:2];
                [enc setBuffer:bufExtAccel offset:0 atIndex:3];
                [enc setBuffer:bufSpringK  offset:0 atIndex:4];
                {
                    uint32_t nb = n_bonds;
                    id<MTLBuffer> bNb = [ctx.device newBufferWithBytes:&nb
                        length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
                    [enc setBuffer:bNb offset:0 atIndex:5];
                }
                [enc setBuffer:bufNa       offset:0 atIndex:6];
                [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [enc endEncoding];
            }
            if (use_anchors && use_springs) {
                // Sheet anchors: elastic→boundary bonds. Without these,
                // sheets fall under gravity and collapse to floor.
                // Optional muscle drive: anchor rest length oscillates
                // sinusoidally, producing visible periodic contraction.
                float current_t = time_offset_s + (float)step * dt;
                memcpy([bufTimeT contents], &current_t, sizeof(float));
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_anchor];
                [enc setBuffer:bufPosOld     offset:0 atIndex:0];
                [enc setBuffer:bufAnchors    offset:0 atIndex:1];
                [enc setBuffer:bufExtAccel   offset:0 atIndex:2];
                [enc setBuffer:bufAnchorK    offset:0 atIndex:3];
                [enc setBuffer:bufNanchors   offset:0 atIndex:4];
                [enc setBuffer:bufNa         offset:0 atIndex:5];
                [enc setBuffer:bufMuscleFreq offset:0 atIndex:6];
                [enc setBuffer:bufMuscleAmp  offset:0 atIndex:7];
                [enc setBuffer:bufTimeT      offset:0 atIndex:8];
                [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [enc endEncoding];
            }
            if (use_pressure) {
                // SPH pressure force counters surface-tension cohesion
                // (M12). Uses density from PREVIOUS step (one-step lag,
                // acceptable since density changes slowly per step).
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_pressure];
                [enc setBuffer:bufPosOld    offset:0 atIndex:0];
                [enc setBuffer:bufDens      offset:0 atIndex:1];
                [enc setBuffer:bufExtAccel  offset:0 atIndex:2];
                [enc setBuffer:bufH         offset:0 atIndex:3];
                [enc setBuffer:bufH2        offset:0 atIndex:4];
                [enc setBuffer:bufMass      offset:0 atIndex:5];
                [enc setBuffer:bufSimScale  offset:0 atIndex:6];
                [enc setBuffer:bufRho       offset:0 atIndex:7];
                [enc setBuffer:bufPressureK offset:0 atIndex:8];
                [enc setBuffer:bufNa        offset:0 atIndex:9];
                [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [enc endEncoding];
            }
            if (need_ext_accel) {
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_apply_ext];
                [enc setBuffer:bufVel       offset:0 atIndex:0];
                [enc setBuffer:bufExtAccel  offset:0 atIndex:1];
                [enc setBuffer:bufDt        offset:0 atIndex:2];
                [enc setBuffer:bufNa        offset:0 atIndex:3];
                [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [enc endEncoding];
            }

            // ── 1. predict_positions ──
            {
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_predict];
                [enc setBuffer:bufPosOld       offset:0 atIndex:0];
                [enc setBuffer:bufVel          offset:0 atIndex:1];
                [enc setBuffer:bufPosPred      offset:0 atIndex:2];
                [enc setBuffer:bufDt           offset:0 atIndex:3];
                [enc setBuffer:bufG            offset:0 atIndex:4];
                [enc setBuffer:bufNa           offset:0 atIndex:5];
                [enc setBuffer:bufSimScaleInv  offset:0 atIndex:6];
                [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [enc endEncoding];
            }

            // Reset Lagrange multipliers to zero at start of step.
            // (The memset is safe here because the previous step's
            // [cmd waitUntilCompleted] already finished.)
            memset([bufLambda contents], 0, dens_bytes);
            if (n_bonds > 0) {
                memset([bufBondLam contents], 0,
                       (size_t)n_bonds * sizeof(float));
            }

            // PERF MEGA-FUSED — distance computation moved INTO the
            // density_grad_combined inner-loop kernel. We no longer
            // materialize r²_aa / r²_as. Each inner iter computes
            // distances inline from current positions (more accurate
            // than the prior distance-reuse trick) and produces
            // density + grad_C + denom_helper in one kernel pass.

            // ── 2. inner XPBD iterations ──
            for (uint32_t it = 0; it < n_iters; it++) {
                // PERF MEGA-GRID kernel: spatial grid lookup for static
                // neighbors (skip 99% non-neighbor pairs at demo1 scale)
                // + dense N² for active-active.
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_d_grad_grid];
                    [enc setBuffer:bufPosPred     offset:0 atIndex:0];
                    [enc setBuffer:bufSortedStatic offset:0 atIndex:1];
                    [enc setBuffer:bufCellStart   offset:0 atIndex:2];
                    [enc setBuffer:bufDens        offset:0 atIndex:3];
                    [enc setBuffer:bufGradC       offset:0 atIndex:4];
                    [enc setBuffer:bufDenomH      offset:0 atIndex:5];
                    [enc setBuffer:bufH           offset:0 atIndex:6];
                    [enc setBuffer:bufH2          offset:0 atIndex:7];
                    [enc setBuffer:bufPoly6       offset:0 atIndex:8];
                    [enc setBuffer:bufSpiky       offset:0 atIndex:9];
                    [enc setBuffer:bufMass        offset:0 atIndex:10];
                    [enc setBuffer:bufRho         offset:0 atIndex:11];
                    [enc setBuffer:bufNa          offset:0 atIndex:12];
                    [enc setBuffer:bufGridDim     offset:0 atIndex:13];
                    [enc setBuffer:bufGridOrigin  offset:0 atIndex:14];
                    // PERF: with grid lookup each row has only ~30
                    // neighbors of work — 32 threads (1 simdgroup) per
                    // row is enough; smaller TG = more rows in flight
                    // simultaneously, better GPU saturation.
                    [enc dispatchThreads:MTLSizeMake(32, n_active, 1)
                        threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
                    [enc endEncoding];
                }
                // i. solve_density_constraint (uses denom_helper)
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_solve_d];
                    [enc setBuffer:bufPosPred offset:0 atIndex:0];
                    [enc setBuffer:bufLambda  offset:0 atIndex:1];
                    [enc setBuffer:bufDens    offset:0 atIndex:2];
                    [enc setBuffer:bufGradC   offset:0 atIndex:3];
                    [enc setBuffer:bufDenomH  offset:0 atIndex:4];
                    [enc setBuffer:bufRho     offset:0 atIndex:5];
                    [enc setBuffer:bufMass    offset:0 atIndex:6];
                    [enc setBuffer:bufAdt2    offset:0 atIndex:7];
                    [enc setBuffer:bufNa      offset:0 atIndex:8];
                    [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [enc endEncoding];
                }
                // j. solve_distance_constraints (elastic bonds, sequential)
                //    Skipped when springs are on — Hooke forces handle bonds.
                if (use_xpbd_bonds) {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_solve_b];
                    [enc setBuffer:bufPosPred offset:0 atIndex:0];
                    [enc setBuffer:bufBondLam offset:0 atIndex:1];
                    [enc setBuffer:bufBondIJ  offset:0 atIndex:2];
                    [enc setBuffer:bufBondRest offset:0 atIndex:3];
                    [enc setBuffer:bufAdistDt2 offset:0 atIndex:4];
                    [enc setBuffer:bufMassInv  offset:0 atIndex:5];
                    [enc setBuffer:bufNbonds   offset:0 atIndex:6];
                    [enc dispatchThreads:MTLSizeMake(1, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(1, 1, 1)];
                    [enc endEncoding];
                }
                // k. solve_floor_constraint (elastic with restitution)
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_solve_f];
                    [enc setBuffer:bufPosPred offset:0 atIndex:0];
                    [enc setBuffer:bufFloor   offset:0 atIndex:1];
                    [enc setBuffer:bufNa      offset:0 atIndex:2];
                    [enc setBuffer:bufRestitution offset:0 atIndex:3];
                    [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [enc endEncoding];
                }
            }

            // ── 3. update_velocities ──
            // Note: must run BEFORE membrane apply so velocity is computed
            // from the post-floor (pre-membrane) position. Otherwise the
            // membrane's position bump (~r0/2 per step) gets converted into
            // a huge velocity injection (delta/dt) that rockets liquid into
            // space. Membrane is a position-only soft constraint here; it
            // doesn't change momentum.
            {
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_updvel];
                [enc setBuffer:bufVel        offset:0 atIndex:0];
                [enc setBuffer:bufPosOld     offset:0 atIndex:1];
                [enc setBuffer:bufPosPred    offset:0 atIndex:2];
                [enc setBuffer:bufDt         offset:0 atIndex:3];
                [enc setBuffer:bufNa         offset:0 atIndex:4];
                [enc setBuffer:bufSimScale   offset:0 atIndex:5];
                [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [enc endEncoding];
            }
            // ── 3.1 Velocity clamp ── prevents CFL teleport past walls.
            if (use_vel_clamp) {
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_clampv];
                [enc setBuffer:bufVel       offset:0 atIndex:0];
                [enc setBuffer:bufVelClamp  offset:0 atIndex:1];
                [enc setBuffer:bufNa        offset:0 atIndex:2];
                [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [enc endEncoding];
            }
            // ── 3.2 Box clamp ── hard wall — particles outside sim box
            // get pushed back to the boundary face with zeroed normal velocity.
            if (use_box_clamp) {
                id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                [enc setComputePipelineState:pso_clampbox];
                [enc setBuffer:bufPosPred offset:0 atIndex:0];
                [enc setBuffer:bufVel     offset:0 atIndex:1];
                [enc setBuffer:bufBoxMin  offset:0 atIndex:2];
                [enc setBuffer:bufBoxMax  offset:0 atIndex:3];
                [enc setBuffer:bufNa      offset:0 atIndex:4];
                [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                [enc endEncoding];
            }

            // ── 3.5 Membrane interaction (M10) ──
            // Runs LAST in the step (after update_velocities) so it acts
            // as a position-only correction. The pos buffer at this point
            // is x_post_floor; we accumulate mem_corr based on those
            // positions and then add it. update_velocities already saw
            // x_post_floor so velocity is unaffected by the bump.
            if (use_membranes) {
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_mem_clear];
                    [enc setBuffer:bufMemCorr offset:0 atIndex:0];
                    [enc setBuffer:bufNa      offset:0 atIndex:1];
                    [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [enc endEncoding];
                }
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_mem_acc];
                    [enc setBuffer:bufPosPred   offset:0 atIndex:0];
                    [enc setBuffer:bufMembranes offset:0 atIndex:1];
                    [enc setBuffer:bufPmemIdx   offset:0 atIndex:2];
                    [enc setBuffer:bufMemCorr   offset:0 atIndex:3];
                    [enc setBuffer:bufNa        offset:0 atIndex:4];
                    [enc setBuffer:bufNelastic  offset:0 atIndex:5];
                    [enc setBuffer:bufR0        offset:0 atIndex:6];
                    [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [enc endEncoding];
                }
                {
                    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
                    [enc setComputePipelineState:pso_mem_apply];
                    [enc setBuffer:bufPosPred offset:0 atIndex:0];
                    [enc setBuffer:bufMemCorr offset:0 atIndex:1];
                    [enc setBuffer:bufNa      offset:0 atIndex:2];
                    [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
                        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
                    [enc endEncoding];
                }
            }

            [cmd commit];
            [cmd waitUntilCompleted];

            // For multi-step bench: x_new becomes x_old for next step.
            memcpy([bufPosOld contents], [bufPosPred contents], pos_a_bytes);
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_step_ms = (t1 - t0) * 1000.0 / bench_steps;
        fprintf(stderr, "xpbd_step: %d steps × %u inner iters, "
                        "%.3f ms/step (steady state)\n",
                bench_steps, n_iters, per_step_ms);

        // Write outputs (final state).
        FILE *fp = fopen(out_pos_path, "wb");
        if (!fp) { fprintf(stderr, "cannot open %s\n", out_pos_path); return 1; }
        fwrite([bufPosPred contents], 1, pos_a_bytes, fp);
        fclose(fp);
        FILE *fv = fopen(out_vel_path, "wb");
        if (!fv) { fprintf(stderr, "cannot open %s\n", out_vel_path); return 1; }
        fwrite([bufVel contents], 1, pos_a_bytes, fv);
        fclose(fv);
    }
    free(pos_active_init); free(vel_active_init); free(pos_static);
    free(sorted_static); free(cell_start);
    if (bonds_raw) free(bonds_raw);
    if (mem_tris) free(mem_tris);
    if (mem_pidx) free(mem_pidx);
    if (anchors_data) free(anchors_data);
    return 0;
}

