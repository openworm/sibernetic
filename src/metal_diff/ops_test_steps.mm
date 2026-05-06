// ops_test_steps.mm — kernel-chain test ops.
//
// Each pair (step_*_fwd / step_*_bwd) exercises a small slice of the
// XPBD pipeline against an FD-validatable Python reference:
//   step_simple   = predict + update_vel (no constraints)
//   step_floor    = predict + floor_constraint + update_vel
//   step_bond     = predict + distance_constraint + update_vel

#include "metal_common.h"

// ──────────────────────────────────────────────────────────────────────
// ────────────────────────────────────────────────────────────────────
// M7.C floor — backward of solve_floor_constraint, with mask-emitting
// forward variant so the backward knows which particles were clamped.
// ────────────────────────────────────────────────────────────────────
// Forward variant of solve_floor_constraint that ALSO writes a per-
// particle "was-clamped" mask buffer. The original solve_floor_constraint
// stays unchanged so xpbd_step (M7.A+B forward-only path) doesn't need
// to allocate the mask buffer.
//
// (Defined inside the embedded shader source — see KERNELS block above.)
//
// Host runner that uses the with-mask forward + backward lives further
// down (run_step_floor_fwd / run_step_floor_bwd).
// ────────────────────────────────────────────────────────────────────

// M7.C — step_simple_fwd: ONE step of predict + update_velocities only.
// No constraints. Inputs: x_old, v_old. Outputs: x_pred, v_new.
//
// Used for validating the M7.C backward kernels via numerical-gradient
// comparison. The simplified pipeline has known closed-form derivatives
// so we can check the hand-derived backward matches finite-difference.
// ──────────────────────────────────────────────────────────────────────
int run_step_simple_fwd(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr,
            "usage: sib_metal step_simple_fwd "
            "<n> <dt> <gravity_y> <pos_old.bin> <vel_old.bin> "
            "<pos_pred_out.bin> <vel_new_out.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    float dt = (float)atof(argv[3]);
    float g_y = (float)atof(argv[4]);

    FILE *fp = fopen(argv[5], "rb");
    FILE *fv = fopen(argv[6], "rb");
    if (!fp || !fv) { fprintf(stderr, "cannot open input\n"); return 1; }
    float *x_old = (float *)malloc((size_t)n * 3 * sizeof(float));
    float *v_old = (float *)malloc((size_t)n * 3 * sizeof(float));
    fread(x_old, sizeof(float), n * 3, fp);
    fread(v_old, sizeof(float), n * 3, fv);
    fclose(fp); fclose(fv);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_predict = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_updvel  = make_pso(ctx, "update_velocities");

        size_t nb = (size_t)n * 3 * sizeof(float);
        id<MTLBuffer> bX = [ctx.device newBufferWithBytes:x_old length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bV = [ctx.device newBufferWithBytes:v_old length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bXp = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bVn = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDt  = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bG   = [ctx.device newBufferWithBytes:&g_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN   = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        // Identity unit-scale (toy convention: pos & vel share unit system).
        float one = 1.0f;
        id<MTLBuffer> bOne = [ctx.device newBufferWithBytes:&one length:sizeof(float)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // predict
        {
            id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
            [e setComputePipelineState:pso_predict];
            [e setBuffer:bX   offset:0 atIndex:0];
            [e setBuffer:bV   offset:0 atIndex:1];
            [e setBuffer:bXp  offset:0 atIndex:2];
            [e setBuffer:bDt  offset:0 atIndex:3];
            [e setBuffer:bG   offset:0 atIndex:4];
            [e setBuffer:bN   offset:0 atIndex:5];
            [e setBuffer:bOne offset:0 atIndex:6];
            [e dispatchThreads:MTLSizeMake(n, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
            [e endEncoding];
        }
        // update_vel
        {
            id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
            [e setComputePipelineState:pso_updvel];
            [e setBuffer:bVn  offset:0 atIndex:0];
            [e setBuffer:bX   offset:0 atIndex:1];
            [e setBuffer:bXp  offset:0 atIndex:2];
            [e setBuffer:bDt  offset:0 atIndex:3];
            [e setBuffer:bN   offset:0 atIndex:4];
            [e setBuffer:bOne offset:0 atIndex:5];
            [e dispatchThreads:MTLSizeMake(n, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
            [e endEncoding];
        }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fxp = fopen(argv[7], "wb");
        fwrite([bXp contents], 1, nb, fxp); fclose(fxp);
        FILE *fvn = fopen(argv[8], "wb");
        fwrite([bVn contents], 1, nb, fvn); fclose(fvn);
    }
    free(x_old); free(v_old);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M7.C — step_simple_bwd: backward of the simplified step.
// Inputs: grad_x_pred, grad_v_new (both [n,3] fp32 arrays)
// Outputs: grad_x_old, grad_v_old (both [n,3] fp32), grad_g_y_per (n fp32)
// Host can sum grad_g_y_per to get scalar grad_g_y.
// ──────────────────────────────────────────────────────────────────────
int run_step_simple_bwd(int argc, char **argv) {
    // op + 7 args = 9 total (argv[0]=binary, argv[1]=op,
    //                        argv[2..8] = n, dt, 2 inputs, 3 outputs).
    if (argc != 9) {
        fprintf(stderr,
            "usage: sib_metal step_simple_bwd "
            "<n> <dt> <grad_x_pred.bin> <grad_v_new.bin> "
            "<grad_x_old_out.bin> <grad_v_old_out.bin> "
            "<grad_g_y_out.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    float dt = (float)atof(argv[3]);

    FILE *fxp = fopen(argv[4], "rb");
    FILE *fvn = fopen(argv[5], "rb");
    float *g_xp = (float *)malloc((size_t)n * 3 * sizeof(float));
    float *g_vn = (float *)malloc((size_t)n * 3 * sizeof(float));
    fread(g_xp, sizeof(float), n * 3, fxp);
    fread(g_vn, sizeof(float), n * 3, fvn);
    fclose(fxp); fclose(fvn);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pbw  = make_pso(ctx, "predict_positions_backward");
        id<MTLComputePipelineState> pso_uvbw = make_pso(ctx, "update_velocities_backward");

        size_t nb = (size_t)n * 3 * sizeof(float);
        size_t ns = (size_t)n * sizeof(float);
        // Output gradients start at zero (kernels accumulate via +=).
        id<MTLBuffer> bGxp = [ctx.device newBufferWithBytes:g_xp length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGvn = [ctx.device newBufferWithBytes:g_vn length:nb
            options:MTLResourceStorageModeShared];
        // Initialize grad outputs to zero (alloc with newBufferWithLength
        // is undefined per Metal docs — explicitly memset for safety).
        id<MTLBuffer> bGxo = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGvo = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGgy = [ctx.device newBufferWithLength:ns
            options:MTLResourceStorageModeShared];
        memset([bGxo contents], 0, nb);
        memset([bGvo contents], 0, nb);
        memset([bGgy contents], 0, ns);

        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN  = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        // Identity unit-scale (toy convention).
        float one_pbw = 1.0f;
        id<MTLBuffer> bOne_pbw = [ctx.device newBufferWithBytes:&one_pbw length:sizeof(float)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // Reverse order: backward through update_velocities first, then predict.
        {
            id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
            [e setComputePipelineState:pso_uvbw];
            [e setBuffer:bGvn offset:0 atIndex:0];
            [e setBuffer:bGxp offset:0 atIndex:1];   // accumulate into existing grad_x_pred
            [e setBuffer:bGxo offset:0 atIndex:2];   // accumulate into grad_x_old
            [e setBuffer:bDt  offset:0 atIndex:3];
            [e setBuffer:bN   offset:0 atIndex:4];
            [e dispatchThreads:MTLSizeMake(n, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
            [e endEncoding];
        }
        {
            id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
            [e setComputePipelineState:pso_pbw];
            [e setBuffer:bGxp offset:0 atIndex:0];
            [e setBuffer:bGxo offset:0 atIndex:1];
            [e setBuffer:bGvo offset:0 atIndex:2];
            [e setBuffer:bGgy offset:0 atIndex:3];
            [e setBuffer:bDt  offset:0 atIndex:4];
            [e setBuffer:bN   offset:0 atIndex:5];
            [e setBuffer:bOne_pbw offset:0 atIndex:6];
            [e dispatchThreads:MTLSizeMake(n, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
            [e endEncoding];
        }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *o1 = fopen(argv[6], "wb");
        fwrite([bGxo contents], 1, nb, o1); fclose(o1);
        FILE *o2 = fopen(argv[7], "wb");
        fwrite([bGvo contents], 1, nb, o2); fclose(o2);
        FILE *o3 = fopen(argv[8], "wb");
        fwrite([bGgy contents], 1, ns, o3); fclose(o3);
    }
    free(g_xp); free(g_vn);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M7.C — step_floor_fwd: predict + floor_with_mask + update_velocities
// over K steps. Saves per-step trajectory of positions and clamp masks
// so the backward can replay through them.
//
// Argv layout:
//   [2]=n  [3]=K  [4]=dt  [5]=gravity_y  [6]=floor_y
//   [7]=pos0.bin  [8]=vel0.bin
//   [9]=traj_pos.bin (out, [K+1, n, 3] fp32 — start-of-step positions)
//   [10]=traj_clamp.bin (out, [K, n] int32 — per-step masks)
//   [11]=vel_final.bin (out, [n, 3] fp32 — velocity at end of K steps)
// argc = 12 (op + 10 positional args).
// ──────────────────────────────────────────────────────────────────────
int run_step_floor_fwd(int argc, char **argv) {
    if (argc != 12) {
        fprintf(stderr,
            "usage: sib_metal step_floor_fwd "
            "<n> <K> <dt> <gravity_y> <floor_y> "
            "<pos0.bin> <vel0.bin> "
            "<traj_pos.bin> <traj_clamp.bin> <vel_final.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    uint32_t K = (uint32_t)atoi(argv[3]);
    float dt   = (float)atof(argv[4]);
    float g_y  = (float)atof(argv[5]);
    float fl_y = (float)atof(argv[6]);
    const char *path_x0    = argv[7];
    const char *path_v0    = argv[8];
    const char *path_tpos  = argv[9];
    const char *path_tclmp = argv[10];
    const char *path_vfinal= argv[11];

    size_t nb = (size_t)n * 3 * sizeof(float);
    size_t nm = (size_t)n * sizeof(int32_t);
    FILE *fp = fopen(path_x0, "rb");
    FILE *fv = fopen(path_v0, "rb");
    if (!fp || !fv) { fprintf(stderr, "cannot open input\n"); return 1; }
    float *x0 = (float *)malloc(nb);
    float *v0 = (float *)malloc(nb);
    fread(x0, 1, nb, fp); fread(v0, 1, nb, fv);
    fclose(fp); fclose(fv);

    float *traj_pos = (float *)malloc((size_t)(K + 1) * nb);
    int32_t *traj_clamp = (int32_t *)malloc((size_t)K * nm);
    memcpy(traj_pos, x0, nb);  // frame 0 = initial positions

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pred  = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_floor = make_pso(ctx, "solve_floor_constraint_with_mask");
        id<MTLComputePipelineState> pso_uv    = make_pso(ctx, "update_velocities");

        id<MTLBuffer> bX  = [ctx.device newBufferWithBytes:x0 length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bV  = [ctx.device newBufferWithBytes:v0 length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bXp = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bClmp = [ctx.device newBufferWithLength:nm
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGy = [ctx.device newBufferWithBytes:&g_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bFy = [ctx.device newBufferWithBytes:&fl_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN  = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        // Identity unit-scale (toy convention: pos & vel share unit system).
        float one = 1.0f;
        id<MTLBuffer> bOne = [ctx.device newBufferWithBytes:&one length:sizeof(float)
            options:MTLResourceStorageModeShared];
        // Legacy step_floor_fwd test op uses inelastic floor (e=0).
        float restitution_zero = 0.0f;
        id<MTLBuffer> bRest = [ctx.device newBufferWithBytes:&restitution_zero
            length:sizeof(float) options:MTLResourceStorageModeShared];

        for (uint32_t k = 0; k < K; k++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            // predict
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pred];
              [e setBuffer:bX offset:0 atIndex:0]; [e setBuffer:bV offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bGy offset:0 atIndex:4]; [e setBuffer:bN offset:0 atIndex:5];
              [e setBuffer:bOne offset:0 atIndex:6];
              [e dispatchThreads:MTLSizeMake(n,1,1)
                  threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // floor with mask (legacy: e=0, inelastic)
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_floor];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bClmp offset:0 atIndex:1];
              [e setBuffer:bFy offset:0 atIndex:2]; [e setBuffer:bN offset:0 atIndex:3];
              [e setBuffer:bRest offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n,1,1)
                  threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // update_vel
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uv];
              [e setBuffer:bV offset:0 atIndex:0]; [e setBuffer:bX offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bN offset:0 atIndex:4];
              [e setBuffer:bOne offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n,1,1)
                  threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // Save this step's mask + post-floor positions, advance state.
            memcpy(&traj_clamp[k * n], [bClmp contents], nm);
            memcpy((char *)traj_pos + (size_t)(k + 1) * nb, [bXp contents], nb);
            memcpy([bX contents], [bXp contents], nb);
        }

        FILE *fvf = fopen(path_vfinal, "wb");
        if (!fvf) { fprintf(stderr, "cannot open vel_final out\n"); return 1; }
        fwrite([bV contents], 1, nb, fvf); fclose(fvf);
    }
    FILE *ftp = fopen(path_tpos, "wb");
    fwrite(traj_pos, 1, (size_t)(K + 1) * nb, ftp); fclose(ftp);
    FILE *ftc = fopen(path_tclmp, "wb");
    fwrite(traj_clamp, 1, (size_t)K * nm, ftc); fclose(ftc);

    free(x0); free(v0); free(traj_pos); free(traj_clamp);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M7.C — step_floor_bwd: backward through K steps of (predict + floor +
// update_velocities), given trajectory and masks from forward.
//
// Inputs:  ∂L/∂x_final (gradient on final positions, [n, 3] fp32)
// Outputs: ∂L/∂x_init, ∂L/∂v_init, ∂L/∂g_y (scalar), ∂L/∂floor_y (scalar)
//
// Argv layout:
//   [2]=n  [3]=K  [4]=dt
//   [5]=traj_pos.bin  [6]=traj_clamp.bin
//   [7]=grad_x_final.bin
//   [8]=grad_x_init.bin (out)  [9]=grad_v_init.bin (out)
//   [10]=grad_g_y.bin (out, single fp32)  [11]=grad_floor_y.bin (out, single fp32)
// argc = 12.
// ──────────────────────────────────────────────────────────────────────
int run_step_floor_bwd(int argc, char **argv) {
    if (argc != 12) {
        fprintf(stderr,
            "usage: sib_metal step_floor_bwd "
            "<n> <K> <dt> <traj_pos.bin> <traj_clamp.bin> "
            "<grad_x_final.bin> "
            "<grad_x_init.bin> <grad_v_init.bin> "
            "<grad_g_y.bin> <grad_floor_y.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    uint32_t K = (uint32_t)atoi(argv[3]);
    float dt = (float)atof(argv[4]);

    size_t nb = (size_t)n * 3 * sizeof(float);
    size_t nm = (size_t)n * sizeof(int32_t);
    size_t ns = (size_t)n * sizeof(float);

    FILE *ftp = fopen(argv[5], "rb");
    FILE *ftc = fopen(argv[6], "rb");
    FILE *fgx = fopen(argv[7], "rb");
    if (!ftp || !ftc || !fgx) { fprintf(stderr, "cannot open input\n"); return 1; }
    float   *traj_pos   = (float *)malloc((size_t)(K + 1) * nb);
    int32_t *traj_clamp = (int32_t *)malloc((size_t)K * nm);
    float   *grad_x_fin = (float *)malloc(nb);
    fread(traj_pos,   1, (size_t)(K + 1) * nb, ftp);
    fread(traj_clamp, 1, (size_t)K * nm, ftc);
    fread(grad_x_fin, 1, nb, fgx);
    fclose(ftp); fclose(ftc); fclose(fgx);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pbw   = make_pso(ctx, "predict_positions_backward");
        id<MTLComputePipelineState> pso_uvbw  = make_pso(ctx, "update_velocities_backward");
        id<MTLComputePipelineState> pso_flbw  = make_pso(ctx, "solve_floor_constraint_backward");

        // Persistent gradient buffers (accumulate across timesteps).
        id<MTLBuffer> bGx_old = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];  // ∂L/∂x_old (running)
        id<MTLBuffer> bGv_old = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];  // ∂L/∂v_old (running)
        // Per-step working buffers.
        id<MTLBuffer> bGx_pred = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGv_new  = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGx_pre_floor = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        // Per-particle parameter gradients (summed host-side).
        id<MTLBuffer> bGgy_per = [ctx.device newBufferWithLength:ns
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGfy_per = [ctx.device newBufferWithLength:ns
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bClmp = [ctx.device newBufferWithLength:nm
            options:MTLResourceStorageModeShared];

        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN  = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        // Initialize: at the END of the simulation, gradient on the final
        // post-floor position is grad_x_fin. v_new gradient = 0 (loss only
        // depends on final positions in this demo). x_old grad = 0 too.
        memcpy([bGx_old contents], grad_x_fin, nb);  // initial x_old grad = grad_x_fin
        memset([bGv_old contents], 0, nb);
        // Accumulators for parameter gradients.
        float total_g_y = 0.0f;
        float total_fy = 0.0f;

        // Walk backwards through K steps.
        for (int32_t k = (int32_t)K - 1; k >= 0; k--) {
            // At this step's BACKWARD entry: bGx_old holds ∂L/∂x_post_floor
            // for step k (i.e., ∂L/∂position at end of step k = start of step k+1).
            // We need to chain through: update_vel ← floor ← predict.

            // Reset per-step working grads.
            memset([bGx_pred contents], 0, nb);
            memset([bGv_new contents], 0, nb);
            memcpy([bClmp contents], &traj_clamp[k * n], nm);

            // Note: in the forward, after this step's update_vel:
            //   v_new[k] is the "vel" at end of step (input to next step's predict)
            //   But our loss only depends on FINAL positions (step K). So
            //   ∂L/∂v_new[K-1] = 0. For interior steps, the gradient flows
            //   from next step's ∂L/∂v_old (which we'd accumulate into here).
            //   In our simplified driver: we maintain bGv_old as the running
            //   "what was ∂L/∂v at next-step-start" i.e. ∂L/∂v_new[k].
            // So: ∂L/∂v_new[k] = bGv_old (which is the running v gradient from
            //                              prior backward steps, or 0 at start).

            // bGx_old at this point = ∂L/∂x_post_floor[k] from prior steps.
            // We need to also incorporate any ∂L/∂x_old that came from this
            // step's update_velocities forward (which used x_old[k] as input).
            // → handled by update_velocities_backward.

            // ── Backward through update_velocities ──
            // Forward: v_new = (x_post_floor - x_old) / dt
            // Inputs to grad: ∂L/∂v_new (in bGv_old), need to flow into
            //   ∂L/∂x_post_floor (add to bGx_old) and ∂L/∂x_old (output).
            // Note: bGx_old accumulates the "x_post_floor" gradient from
            // both the next step's x_old grad AND this step's update_vel
            // contribution. That's exactly what we want.
            // After update_vel backward, we will use bGx_old as
            // ∂L/∂x_post_floor and accumulate ∂L/∂x_old separately into
            // a fresh buffer (then promote it to bGx_old at end of step).
            id<MTLBuffer> bGx_old_new = [ctx.device newBufferWithLength:nb
                options:MTLResourceStorageModeShared];
            memset([bGx_old_new contents], 0, nb);

            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uvbw];
              [e setBuffer:bGv_old   offset:0 atIndex:0];   // ∂L/∂v_new
              [e setBuffer:bGx_old   offset:0 atIndex:1];   // accumulate into ∂L/∂x_post_floor
              [e setBuffer:bGx_old_new offset:0 atIndex:2]; // accumulate into ∂L/∂x_old (this step)
              [e setBuffer:bDt offset:0 atIndex:3]; [e setBuffer:bN offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }

            // ── Backward through floor ──
            // Inputs: ∂L/∂x_post_floor (= bGx_old). Outputs: ∂L/∂x_pre_floor
            // (which we'll call bGx_pred since that's pos before floor).
            // Legacy step_floor_bwd uses inelastic floor (e=0).
            float restitution_zero_b = 0.0f;
            id<MTLBuffer> bRestZero = [ctx.device newBufferWithBytes:&restitution_zero_b
                length:sizeof(float) options:MTLResourceStorageModeShared];
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_flbw];
              [e setBuffer:bGx_old offset:0 atIndex:0];   // ∂L/∂x_post_floor
              [e setBuffer:bGx_pred offset:0 atIndex:1];  // accumulate into ∂L/∂x_pre_floor (= ∂L/∂x_pred)
              [e setBuffer:bGfy_per offset:0 atIndex:2];
              [e setBuffer:bClmp offset:0 atIndex:3]; [e setBuffer:bN offset:0 atIndex:4];
              [e setBuffer:bRestZero offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }

            // ── Backward through predict ──
            // Inputs: ∂L/∂x_pred (bGx_pred). Outputs: ∂L/∂x_old (accumulate
            // into bGx_old_new), ∂L/∂v (accumulate into a fresh bGv_old_new),
            // ∂L/∂g_y (per-particle into bGgy_per).
            id<MTLBuffer> bGv_old_new = [ctx.device newBufferWithLength:nb
                options:MTLResourceStorageModeShared];
            memset([bGv_old_new contents], 0, nb);
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pbw];
              [e setBuffer:bGx_pred offset:0 atIndex:0];
              [e setBuffer:bGx_old_new offset:0 atIndex:1];   // accumulate
              [e setBuffer:bGv_old_new offset:0 atIndex:2];   // accumulate
              [e setBuffer:bGgy_per offset:0 atIndex:3];
              [e setBuffer:bDt offset:0 atIndex:4]; [e setBuffer:bN offset:0 atIndex:5];
              float one_pf = 1.0f;
              id<MTLBuffer> bOne_pf = [ctx.device newBufferWithBytes:&one_pf
                  length:sizeof(float) options:MTLResourceStorageModeShared];
              [e setBuffer:bOne_pf offset:0 atIndex:6];
              [e dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }

            [cmd commit]; [cmd waitUntilCompleted];

            // Sum per-particle parameter grads into running totals.
            float *gy = (float *)[bGgy_per contents];
            float *fy = (float *)[bGfy_per contents];
            for (uint32_t i = 0; i < n; i++) { total_g_y += gy[i]; total_fy += fy[i]; }

            // Promote per-step grads to running grads for the next (earlier) step.
            memcpy([bGx_old contents], [bGx_old_new contents], nb);
            memcpy([bGv_old contents], [bGv_old_new contents], nb);
        }

        // Write outputs.
        FILE *o1 = fopen(argv[8], "wb");
        fwrite([bGx_old contents], 1, nb, o1); fclose(o1);
        FILE *o2 = fopen(argv[9], "wb");
        fwrite([bGv_old contents], 1, nb, o2); fclose(o2);
        FILE *o3 = fopen(argv[10], "wb");
        fwrite(&total_g_y, 1, sizeof(float), o3); fclose(o3);
        FILE *o4 = fopen(argv[11], "wb");
        fwrite(&total_fy, 1, sizeof(float), o4); fclose(o4);
    }
    free(traj_pos); free(traj_clamp); free(grad_x_fin);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M8 — step_bond_fwd: K-step forward of (predict + bonds_with_save +
// update_velocities). NO floor for simplicity — focus on isolating the
// bond physics so backward correctness is unambiguous. n_iters fixed at 1.
//
// Saves:
//   traj_pos.bin    — [K+1, n, 3] fp32 positions at start of each step
//   traj_state.bin  — [K, n_bonds*7] fp32 per-bond pre-state per step
//   traj_lambda.bin — [K, n_bonds] fp32 lambda before each step's projection
// ──────────────────────────────────────────────────────────────────────
int run_step_bond_fwd(int argc, char **argv) {
    if (argc != 14) {
        fprintf(stderr,
            "usage: sib_metal step_bond_fwd "
            "<n> <K> <n_bonds> <dt> <gravity_y> <mass> <alpha_dist> "
            "<pos0.bin> <vel0.bin> <bonds.bin> "
            "<traj_pos.bin> <traj_state.bin> <vel_final.bin>\n");
        return 1;
    }
    uint32_t n        = (uint32_t)atoi(argv[2]);
    uint32_t K        = (uint32_t)atoi(argv[3]);
    uint32_t n_bonds  = (uint32_t)atoi(argv[4]);
    float dt          = (float)atof(argv[5]);
    float g_y         = (float)atof(argv[6]);
    float mass        = (float)atof(argv[7]);
    float alpha_dist  = (float)atof(argv[8]);
    const char *p_x0    = argv[9];
    const char *p_v0    = argv[10];
    const char *p_bonds = argv[11];
    const char *p_tpos  = argv[12];
    const char *p_tstate = argv[13];
    const char *p_vfin  = "/tmp/sibm_step_bond_vfinal.bin";  // fixed for now

    float alpha_inv_dt2 = alpha_dist / (dt * dt);
    float mass_inv = 1.0f / mass;

    size_t nb = (size_t)n * 3 * sizeof(float);
    size_t state_bytes_per_step = (size_t)n_bonds * 7 * sizeof(float);

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        if (fread(buf, sizeof(float), n_floats, f) != n_floats) {
            fprintf(stderr, "short read on %s\n", path); exit(1);
        }
        fclose(f); return buf;
    };
    float *x0 = read_floats(p_x0, (size_t)n * 3);
    float *v0 = read_floats(p_v0, (size_t)n * 3);

    // Read bonds and unpack to (i,j) and rest_len arrays.
    FILE *fb = fopen(p_bonds, "rb");
    if (!fb) { fprintf(stderr, "cannot open %s\n", p_bonds); return 1; }
    void *bonds_raw = malloc((size_t)n_bonds * 16);
    fread(bonds_raw, 16, n_bonds, fb); fclose(fb);
    int32_t *bond_ij = (int32_t *)malloc((size_t)n_bonds * 2 * sizeof(int32_t));
    float *bond_rest = (float *)malloc((size_t)n_bonds * sizeof(float));
    for (uint32_t b = 0; b < n_bonds; b++) {
        memcpy(&bond_ij[b * 2], (uint8_t *)bonds_raw + b * 16, 8);
        memcpy(&bond_rest[b],   (uint8_t *)bonds_raw + b * 16 + 8, 4);
    }
    free(bonds_raw);

    float *traj_pos = (float *)malloc((size_t)(K + 1) * nb);
    float *traj_state = (float *)malloc((size_t)K * state_bytes_per_step);
    memcpy(traj_pos, x0, nb);  // frame 0

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pred = make_pso(ctx, "predict_positions");
        id<MTLComputePipelineState> pso_bond = make_pso(ctx, "solve_distance_constraints_seq_with_save");
        id<MTLComputePipelineState> pso_uv   = make_pso(ctx, "update_velocities");

        id<MTLBuffer> bX = [ctx.device newBufferWithBytes:x0 length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bV = [ctx.device newBufferWithBytes:v0 length:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bXp = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bBondIJ = [ctx.device newBufferWithBytes:bond_ij
            length:(size_t)n_bonds * 2 * sizeof(int32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bBondRest = [ctx.device newBufferWithBytes:bond_rest
            length:(size_t)n_bonds * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bLambda = [ctx.device newBufferWithLength:(size_t)n_bonds * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bState = [ctx.device newBufferWithLength:state_bytes_per_step
            options:MTLResourceStorageModeShared];

        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGy = [ctx.device newBufferWithBytes:&g_y length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bAdt2 = [ctx.device newBufferWithBytes:&alpha_inv_dt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bMinv = [ctx.device newBufferWithBytes:&mass_inv length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNb = [ctx.device newBufferWithBytes:&n_bonds length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        // Identity unit-scale (toy convention: pos & vel share unit system).
        float one = 1.0f;
        id<MTLBuffer> bOne = [ctx.device newBufferWithBytes:&one length:sizeof(float)
            options:MTLResourceStorageModeShared];

        for (uint32_t k = 0; k < K; k++) {
            // Reset lambdas to zero at start of each step.
            memset([bLambda contents], 0, (size_t)n_bonds * sizeof(float));

            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            // predict
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pred];
              [e setBuffer:bX  offset:0 atIndex:0]; [e setBuffer:bV  offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bGy offset:0 atIndex:4]; [e setBuffer:bN  offset:0 atIndex:5];
              [e setBuffer:bOne offset:0 atIndex:6];
              [e dispatchThreads:MTLSizeMake(n,1,1)
                  threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            // bonds (n_iters=1)
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_bond];
              [e setBuffer:bXp offset:0 atIndex:0]; [e setBuffer:bLambda offset:0 atIndex:1];
              [e setBuffer:bBondIJ offset:0 atIndex:2]; [e setBuffer:bBondRest offset:0 atIndex:3];
              [e setBuffer:bState offset:0 atIndex:4]; [e setBuffer:bAdt2 offset:0 atIndex:5];
              [e setBuffer:bMinv offset:0 atIndex:6]; [e setBuffer:bNb offset:0 atIndex:7];
              [e dispatchThreads:MTLSizeMake(1,1,1)
                  threadsPerThreadgroup:MTLSizeMake(1,1,1)];
              [e endEncoding]; }
            // update_vel
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uv];
              [e setBuffer:bV offset:0 atIndex:0]; [e setBuffer:bX offset:0 atIndex:1];
              [e setBuffer:bXp offset:0 atIndex:2]; [e setBuffer:bDt offset:0 atIndex:3];
              [e setBuffer:bN offset:0 atIndex:4];
              [e setBuffer:bOne offset:0 atIndex:5];
              [e dispatchThreads:MTLSizeMake(n,1,1)
                  threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // Save state for this step, advance positions.
            memcpy((char *)traj_state + (size_t)k * state_bytes_per_step,
                   [bState contents], state_bytes_per_step);
            memcpy((char *)traj_pos + (size_t)(k + 1) * nb,
                   [bXp contents], nb);
            memcpy([bX contents], [bXp contents], nb);
        }

        FILE *fvf = fopen(p_vfin, "wb");
        fwrite([bV contents], 1, nb, fvf); fclose(fvf);
    }
    FILE *ftp = fopen(p_tpos, "wb");
    fwrite(traj_pos, 1, (size_t)(K + 1) * nb, ftp); fclose(ftp);
    FILE *fts = fopen(p_tstate, "wb");
    fwrite(traj_state, 1, (size_t)K * state_bytes_per_step, fts); fclose(fts);

    free(x0); free(v0); free(bond_ij); free(bond_rest);
    free(traj_pos); free(traj_state);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M8 — step_bond_bwd: backward through K-step bond simulation. Walks
// trajectory in REVERSE, calling per-step backward kernels in order:
// update_vel_backward → bond_backward → predict_backward.
//
// Inputs:
//   grad_x_final.bin — [n, 3] fp32 ∂L/∂(positions at end)
//   plus traj_pos, traj_state, bonds from forward
//
// Outputs:
//   grad_x_init, grad_v_init, grad_alpha (single-element fp32)
// ──────────────────────────────────────────────────────────────────────
int run_step_bond_bwd(int argc, char **argv) {
    if (argc != 16) {
        fprintf(stderr,
            "usage: sib_metal step_bond_bwd "
            "<n> <K> <n_bonds> <dt> <mass> <alpha_dist> "
            "<bonds.bin> <traj_pos.bin> <traj_state.bin> <grad_x_final.bin> "
            "<grad_x_init.bin> <grad_v_init.bin> <grad_alpha.bin>\n");
        return 1;
    }
    uint32_t n       = (uint32_t)atoi(argv[2]);
    uint32_t K       = (uint32_t)atoi(argv[3]);
    uint32_t n_bonds = (uint32_t)atoi(argv[4]);
    float dt         = (float)atof(argv[5]);
    float mass       = (float)atof(argv[6]);
    float alpha_dist = (float)atof(argv[7]);
    const char *p_bonds  = argv[8];
    const char *p_tpos   = argv[9];
    const char *p_tstate = argv[10];
    const char *p_gxf    = argv[11];
    const char *p_gxi    = argv[12];
    const char *p_gvi    = argv[13];
    const char *p_galpha = argv[14];
    // argv[15] currently unused — placeholder for grad_L_rest if added.

    float alpha_inv_dt2 = alpha_dist / (dt * dt);
    float dt2 = dt * dt;
    float mass_inv = 1.0f / mass;

    size_t nb = (size_t)n * 3 * sizeof(float);
    size_t state_bytes_per_step = (size_t)n_bonds * 7 * sizeof(float);

    // Read bonds → (ij, rest)
    FILE *fb = fopen(p_bonds, "rb");
    void *bonds_raw = malloc((size_t)n_bonds * 16);
    fread(bonds_raw, 16, n_bonds, fb); fclose(fb);
    int32_t *bond_ij = (int32_t *)malloc((size_t)n_bonds * 2 * sizeof(int32_t));
    float *bond_rest = (float *)malloc((size_t)n_bonds * sizeof(float));
    for (uint32_t b = 0; b < n_bonds; b++) {
        memcpy(&bond_ij[b * 2], (uint8_t *)bonds_raw + b * 16, 8);
        memcpy(&bond_rest[b],   (uint8_t *)bonds_raw + b * 16 + 8, 4);
    }
    free(bonds_raw);

    float *traj_state = (float *)malloc((size_t)K * state_bytes_per_step);
    float *grad_x_fin = (float *)malloc(nb);
    FILE *fts = fopen(p_tstate, "rb");
    FILE *fgx = fopen(p_gxf, "rb");
    fread(traj_state, 1, (size_t)K * state_bytes_per_step, fts);
    fread(grad_x_fin, 1, nb, fgx);
    fclose(fts); fclose(fgx);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_pbw  = make_pso(ctx, "predict_positions_backward");
        id<MTLComputePipelineState> pso_uvbw = make_pso(ctx, "update_velocities_backward");
        id<MTLComputePipelineState> pso_bbw  = make_pso(ctx, "solve_distance_constraints_seq_backward");

        id<MTLBuffer> bGx_old = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGv_old = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGx_pred = [ctx.device newBufferWithLength:nb
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGlam = [ctx.device newBufferWithLength:(size_t)n_bonds * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGalpha = [ctx.device newBufferWithLength:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGgy_per = [ctx.device newBufferWithLength:(size_t)n * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bState = [ctx.device newBufferWithLength:state_bytes_per_step
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bBondIJ = [ctx.device newBufferWithBytes:bond_ij
            length:(size_t)n_bonds * 2 * sizeof(int32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bBondRest = [ctx.device newBufferWithBytes:bond_rest
            length:(size_t)n_bonds * sizeof(float)
            options:MTLResourceStorageModeShared];

        id<MTLBuffer> bDt = [ctx.device newBufferWithBytes:&dt length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bDt2 = [ctx.device newBufferWithBytes:&dt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bAdt2 = [ctx.device newBufferWithBytes:&alpha_inv_dt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bAlpha = [ctx.device newBufferWithBytes:&alpha_dist length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bMinv = [ctx.device newBufferWithBytes:&mass_inv length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNb = [ctx.device newBufferWithBytes:&n_bonds length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        // Initial: ∂L/∂x_final → ∂L/∂x_old (running across reverse steps).
        memcpy([bGx_old contents], grad_x_fin, nb);
        memset([bGv_old contents], 0, nb);
        memset([bGalpha contents], 0, sizeof(float));

        // Walk backward through K steps.
        for (int32_t k = (int32_t)K - 1; k >= 0; k--) {
            // Per-step lambda gradient resets at start (forward resets λ to 0).
            memset([bGlam contents], 0, (size_t)n_bonds * sizeof(float));
            // Per-step working buffers.
            memset([bGx_pred contents], 0, nb);

            // Load this step's saved bond state.
            memcpy([bState contents],
                   (char *)traj_state + (size_t)k * state_bytes_per_step,
                   state_bytes_per_step);

            // bGx_old_new accumulates ∂L/∂x_old for THIS step (separate
            // buffer so we don't mix with the running bGx_old that holds
            // the post-step gradient).
            id<MTLBuffer> bGx_old_new = [ctx.device newBufferWithLength:nb
                options:MTLResourceStorageModeShared];
            id<MTLBuffer> bGv_old_new = [ctx.device newBufferWithLength:nb
                options:MTLResourceStorageModeShared];
            memset([bGx_old_new contents], 0, nb);
            memset([bGv_old_new contents], 0, nb);

            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            // (1) update_vel backward: ∂L/∂v_new (=bGv_old running) flows into
            //     ∂L/∂x_post_bonds (accumulate into bGx_old) and ∂L/∂x_old
            //     (accumulate into bGx_old_new).
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_uvbw];
              [e setBuffer:bGv_old offset:0 atIndex:0];
              [e setBuffer:bGx_old offset:0 atIndex:1];
              [e setBuffer:bGx_old_new offset:0 atIndex:2];
              [e setBuffer:bDt offset:0 atIndex:3]; [e setBuffer:bN offset:0 atIndex:4];
              [e dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // bGx_old now = ∂L/∂x_post_bonds for this step. The bond_backward
            // overwrites it in place (read pos_grad[i,j], write pos_grad[i,j]).
            cmd = [ctx.queue commandBuffer];
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_bbw];
              [e setBuffer:bGx_old offset:0 atIndex:0];   // pos_grad in/out
              [e setBuffer:bGlam offset:0 atIndex:1];     // lambda_grad in/out
              [e setBuffer:bGalpha offset:0 atIndex:2];   // alpha grad accum
              [e setBuffer:bBondIJ offset:0 atIndex:3];
              [e setBuffer:bBondRest offset:0 atIndex:4];
              [e setBuffer:bState offset:0 atIndex:5];
              [e setBuffer:bAdt2 offset:0 atIndex:6];
              [e setBuffer:bAlpha offset:0 atIndex:7];
              [e setBuffer:bDt2 offset:0 atIndex:8];
              [e setBuffer:bMinv offset:0 atIndex:9];
              [e setBuffer:bNb offset:0 atIndex:10];
              [e dispatchThreads:MTLSizeMake(1,1,1) threadsPerThreadgroup:MTLSizeMake(1,1,1)];
              [e endEncoding]; }
            // bGx_old now = ∂L/∂x_pre_bonds = ∂L/∂x_pred (post-predict).
            // Move into bGx_pred for predict_backward.
            // (We could just pass bGx_old to predict_backward directly,
            // but keep separate for clarity.)
            { id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
              [blit copyFromBuffer:bGx_old sourceOffset:0
                          toBuffer:bGx_pred destinationOffset:0
                              size:nb];
              [blit endEncoding]; }
            { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
              [e setComputePipelineState:pso_pbw];
              [e setBuffer:bGx_pred offset:0 atIndex:0];
              [e setBuffer:bGx_old_new offset:0 atIndex:1];   // accumulate
              [e setBuffer:bGv_old_new offset:0 atIndex:2];   // accumulate
              [e setBuffer:bGgy_per offset:0 atIndex:3];
              [e setBuffer:bDt offset:0 atIndex:4]; [e setBuffer:bN offset:0 atIndex:5];
              float one_pb = 1.0f;
              id<MTLBuffer> bOne_pb = [ctx.device newBufferWithBytes:&one_pb
                  length:sizeof(float) options:MTLResourceStorageModeShared];
              [e setBuffer:bOne_pb offset:0 atIndex:6];
              [e dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
              [e endEncoding]; }
            [cmd commit]; [cmd waitUntilCompleted];

            // Promote per-step gradients to running gradients for previous step.
            memcpy([bGx_old contents], [bGx_old_new contents], nb);
            memcpy([bGv_old contents], [bGv_old_new contents], nb);
        }

        FILE *o1 = fopen(p_gxi, "wb");
        fwrite([bGx_old contents], 1, nb, o1); fclose(o1);
        FILE *o2 = fopen(p_gvi, "wb");
        fwrite([bGv_old contents], 1, nb, o2); fclose(o2);
        FILE *o3 = fopen(p_galpha, "wb");
        fwrite([bGalpha contents], 1, sizeof(float), o3); fclose(o3);
    }
    free(bond_ij); free(bond_rest); free(traj_state); free(grad_x_fin);
    return 0;
}

