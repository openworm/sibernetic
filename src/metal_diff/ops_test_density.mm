// ops_test_density.mm — density-pipeline test ops.
//
// Standalone forward and backward ops for the chain that turns particle
// positions into density into a constraint correction:
//   density_as / density_aa: per-pair W contributions to density
//   density_constraint_grad_bwd: backward of the M6.4 forward op
//   solve_density_constraint_fwd / bwd: the actual XPBD projection

#include "metal_common.h"

// ──────────────────────────────────────────────────────────────────────
// M9.A — density_as_fwd: compute density at each active particle from
// active×static neighbor SPH kernel evaluation. No active-active term;
// kept minimal to validate the density backward chain in isolation.
//
// Pipeline: dist_active_static → wpoly6_inplace → rowsum_density.
// Saves r² (post-distance, pre-Wpoly6) for backward.
// ──────────────────────────────────────────────────────────────────────
int run_density_as_fwd(int argc, char **argv) {
    if (argc != 10) {
        fprintf(stderr,
            "usage: sib_metal density_as_fwd "
            "<n_active> <n_static> <h> <mass> "
            "<active.bin> <static.bin> <density_out.bin> <r2_saved_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h    = (float)atof(argv[4]);
    float mass = (float)atof(argv[5]);

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    uint32_t n_r2 = n_active * n_static;

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        if (fread(buf, sizeof(float), n_floats, f) != n_floats) {
            fprintf(stderr, "short read on %s\n", path); exit(1);
        }
        fclose(f); return buf;
    };
    float *active   = read_floats(argv[6], (size_t)n_active * 3);
    float *static_p = read_floats(argv[7], (size_t)n_static * 3);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_dist = make_pso(ctx, "dist_active_static");
        id<MTLComputePipelineState> pso_wp   = make_pso(ctx, "wpoly6_inplace");
        id<MTLComputePipelineState> pso_rs   = make_pso(ctx, "rowsum_density");

        size_t r2_bytes = (size_t)n_r2 * sizeof(float);
        size_t dens_bytes = (size_t)n_active * sizeof(float);

        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:active
            length:(size_t)n_active * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bS = [ctx.device newBufferWithBytes:static_p
            length:(size_t)n_static * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2 = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bW = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];

        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNs = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNr2 = [ctx.device newBufferWithBytes:&n_r2 length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // distance
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_dist];
          [e setBuffer:bA offset:0 atIndex:0]; [e setBuffer:bS offset:0 atIndex:1];
          [e setBuffer:bR2 offset:0 atIndex:2];
          [e setBuffer:bNa offset:0 atIndex:3]; [e setBuffer:bNs offset:0 atIndex:4];
          [e dispatchThreads:MTLSizeMake(n_active, n_static, 1)
              threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
          [e endEncoding]; }
        // GPU-side blit r2 → W (so we can later transform W in-place
        // while keeping r2 saved for backward).
        { id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
          [blit copyFromBuffer:bR2 sourceOffset:0
                      toBuffer:bW destinationOffset:0
                          size:r2_bytes];
          [blit endEncoding]; }
        // wpoly6 in-place on W
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_wp];
          [e setBuffer:bW offset:0 atIndex:0]; [e setBuffer:bH2 offset:0 atIndex:1];
          [e setBuffer:bP6 offset:0 atIndex:2]; [e setBuffer:bNr2 offset:0 atIndex:3];
          [e dispatchThreads:MTLSizeMake(n_r2, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        // rowsum density
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_rs];
          [e setBuffer:bW offset:0 atIndex:0]; [e setBuffer:bD offset:0 atIndex:1];
          [e setBuffer:bM offset:0 atIndex:2];
          [e setBuffer:bNs offset:0 atIndex:3]; [e setBuffer:bNa offset:0 atIndex:4];
          [e dispatchThreads:MTLSizeMake(256, n_active, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fd = fopen(argv[8], "wb");
        fwrite([bD contents], 1, dens_bytes, fd); fclose(fd);
        FILE *fr = fopen(argv[9], "wb");
        fwrite([bR2 contents], 1, r2_bytes, fr); fclose(fr);
    }
    free(active); free(static_p);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.A — density_as_bwd: backward through density_as_fwd.
//
// Inputs: ∂L/∂density [n_active], saved r², positions
// Outputs: ∂L/∂active [n_active × 3]  (static positions are frozen)
//
// Reverse pipeline:
//   ∂L/∂W   = rowsum_density_backward(∂L/∂density)         (broadcast)
//   ∂L/∂r²  = wpoly6_inplace_backward(∂L/∂W; saved r²)     (in-place)
//   ∂L/∂act = dist_active_static_backward(∂L/∂r²)
// ──────────────────────────────────────────────────────────────────────
int run_density_as_bwd(int argc, char **argv) {
    if (argc != 11) {
        fprintf(stderr,
            "usage: sib_metal density_as_bwd "
            "<n_active> <n_static> <h> <mass> "
            "<active.bin> <static.bin> <r2_saved.bin> "
            "<grad_density.bin> <grad_active_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h    = (float)atof(argv[4]);
    float mass = (float)atof(argv[5]);

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    uint32_t n_r2 = n_active * n_static;

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        fread(buf, sizeof(float), n_floats, f); fclose(f); return buf;
    };
    float *active   = read_floats(argv[6], (size_t)n_active * 3);
    float *static_p = read_floats(argv[7], (size_t)n_static * 3);
    float *r2_saved = read_floats(argv[8], (size_t)n_r2);
    float *grad_d   = read_floats(argv[9], (size_t)n_active);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_rs_bw = make_pso(ctx, "rowsum_density_backward");
        id<MTLComputePipelineState> pso_wp_bw = make_pso(ctx, "wpoly6_inplace_backward");
        id<MTLComputePipelineState> pso_d_bw  = make_pso(ctx, "dist_active_static_backward");

        size_t r2_bytes = (size_t)n_r2 * sizeof(float);
        size_t a_bytes  = (size_t)n_active * 3 * sizeof(float);

        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:active length:a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bS = [ctx.device newBufferWithBytes:static_p
            length:(size_t)n_static * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2 = [ctx.device newBufferWithBytes:r2_saved length:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGd = [ctx.device newBufferWithBytes:grad_d
            length:(size_t)n_active * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGW_or_Gr2 = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGa = [ctx.device newBufferWithLength:a_bytes
            options:MTLResourceStorageModeShared];
        memset([bGa contents], 0, a_bytes);

        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNs = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNr2 = [ctx.device newBufferWithBytes:&n_r2 length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // (1) rowsum_bwd: broadcast grad_density → grad_W (write)
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_rs_bw];
          [e setBuffer:bGd offset:0 atIndex:0];
          [e setBuffer:bGW_or_Gr2 offset:0 atIndex:1];
          [e setBuffer:bM offset:0 atIndex:2];
          [e setBuffer:bNs offset:0 atIndex:3];
          [e dispatchThreads:MTLSizeMake(n_static, n_active, 1)
              threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
          [e endEncoding]; }
        // (2) wpoly6_bwd in-place: grad_W → grad_r²
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_wp_bw];
          [e setBuffer:bR2 offset:0 atIndex:0];
          [e setBuffer:bGW_or_Gr2 offset:0 atIndex:1];
          [e setBuffer:bH2 offset:0 atIndex:2];
          [e setBuffer:bP6 offset:0 atIndex:3];
          [e setBuffer:bNr2 offset:0 atIndex:4];
          [e dispatchThreads:MTLSizeMake(n_r2, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        // (3) dist_bwd: grad_r² → grad_active
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_d_bw];
          [e setBuffer:bA offset:0 atIndex:0];
          [e setBuffer:bS offset:0 atIndex:1];
          [e setBuffer:bGW_or_Gr2 offset:0 atIndex:2];
          [e setBuffer:bGa offset:0 atIndex:3];
          [e setBuffer:bNa offset:0 atIndex:4];
          [e setBuffer:bNs offset:0 atIndex:5];
          [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
          [e endEncoding]; }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fo = fopen("/tmp/density_as_grad_active.bin", "wb");
        fwrite([bGa contents], 1, a_bytes, fo); fclose(fo);
        // Allow caller to override output path via argv[9] if it's
        // different — for simplicity fixed path.
    }
    free(active); free(static_p); free(r2_saved); free(grad_d);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.B — density_aa_fwd: density from active×active SPH neighbors only.
// (No static/boundary contribution; for testing dist_active_active_bwd
// in isolation. Self-contribution at r=0 is included.)
// ──────────────────────────────────────────────────────────────────────
int run_density_aa_fwd(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr,
            "usage: sib_metal density_aa_fwd "
            "<n_active> <h> <mass> "
            "<active.bin> <density_out.bin> <r2_saved_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    float h    = (float)atof(argv[3]);
    float mass = (float)atof(argv[4]);
    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    uint32_t n_r2 = n_active * n_active;

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        fread(buf, sizeof(float), n_floats, f); fclose(f); return buf;
    };
    float *active = read_floats(argv[5], (size_t)n_active * 3);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_dist = make_pso(ctx, "dist_active_active");
        id<MTLComputePipelineState> pso_wp   = make_pso(ctx, "wpoly6_inplace");
        id<MTLComputePipelineState> pso_rs   = make_pso(ctx, "rowsum_density");

        size_t r2_bytes = (size_t)n_r2 * sizeof(float);
        size_t dens_bytes = (size_t)n_active * sizeof(float);
        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:active
            length:(size_t)n_active * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2 = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bW = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD = [ctx.device newBufferWithLength:dens_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNr2 = [ctx.device newBufferWithBytes:&n_r2 length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // dist_active_active
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_dist];
          [e setBuffer:bA offset:0 atIndex:0]; [e setBuffer:bR2 offset:0 atIndex:1];
          [e setBuffer:bNa offset:0 atIndex:2];
          [e dispatchThreads:MTLSizeMake(n_active, n_active, 1)
              threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
          [e endEncoding]; }
        // blit r2 → W (preserve r2 for backward)
        { id<MTLBlitCommandEncoder> b = [cmd blitCommandEncoder];
          [b copyFromBuffer:bR2 sourceOffset:0 toBuffer:bW destinationOffset:0
                       size:r2_bytes];
          [b endEncoding]; }
        // wpoly6 in-place on W
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_wp];
          [e setBuffer:bW offset:0 atIndex:0]; [e setBuffer:bH2 offset:0 atIndex:1];
          [e setBuffer:bP6 offset:0 atIndex:2]; [e setBuffer:bNr2 offset:0 atIndex:3];
          [e dispatchThreads:MTLSizeMake(n_r2, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        // rowsum density (n_cols = n_rows = n_active)
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_rs];
          [e setBuffer:bW offset:0 atIndex:0]; [e setBuffer:bD offset:0 atIndex:1];
          [e setBuffer:bM offset:0 atIndex:2];
          [e setBuffer:bNa offset:0 atIndex:3]; [e setBuffer:bNa offset:0 atIndex:4];
          [e dispatchThreads:MTLSizeMake(256, n_active, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fd = fopen(argv[6], "wb"); fwrite([bD contents], 1, dens_bytes, fd); fclose(fd);
        FILE *fr = fopen(argv[7], "wb"); fwrite([bR2 contents], 1, r2_bytes, fr); fclose(fr);
    }
    free(active);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.B — density_aa_bwd: backward through density_aa_fwd.
// Reverse: rowsum_bwd → wpoly6_bwd → dist_active_active_bwd.
// ──────────────────────────────────────────────────────────────────────
int run_density_aa_bwd(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr,
            "usage: sib_metal density_aa_bwd "
            "<n_active> <h> <mass> "
            "<active.bin> <r2_saved.bin> <grad_density.bin> <grad_active_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    float h    = (float)atof(argv[3]);
    float mass = (float)atof(argv[4]);
    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));
    uint32_t n_r2 = n_active * n_active;

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        fread(buf, sizeof(float), n_floats, f); fclose(f); return buf;
    };
    float *active   = read_floats(argv[5], (size_t)n_active * 3);
    float *r2_saved = read_floats(argv[6], (size_t)n_r2);
    float *grad_d   = read_floats(argv[7], (size_t)n_active);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso_rs_bw = make_pso(ctx, "rowsum_density_backward");
        id<MTLComputePipelineState> pso_wp_bw = make_pso(ctx, "wpoly6_inplace_backward");
        id<MTLComputePipelineState> pso_d_bw  = make_pso(ctx, "dist_active_active_backward");

        size_t r2_bytes = (size_t)n_r2 * sizeof(float);
        size_t a_bytes  = (size_t)n_active * 3 * sizeof(float);
        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:active length:a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2 = [ctx.device newBufferWithBytes:r2_saved length:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGd = [ctx.device newBufferWithBytes:grad_d
            length:(size_t)n_active * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGW_or_Gr2 = [ctx.device newBufferWithLength:r2_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGa = [ctx.device newBufferWithLength:a_bytes
            options:MTLResourceStorageModeShared];
        memset([bGa contents], 0, a_bytes);
        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bP6 = [ctx.device newBufferWithBytes:&poly6_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNr2 = [ctx.device newBufferWithBytes:&n_r2 length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        // (1) rowsum_bwd: broadcast grad_density → grad_W
        // n_cols=n_active here.
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_rs_bw];
          [e setBuffer:bGd offset:0 atIndex:0];
          [e setBuffer:bGW_or_Gr2 offset:0 atIndex:1];
          [e setBuffer:bM offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
          [e dispatchThreads:MTLSizeMake(n_active, n_active, 1)
              threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
          [e endEncoding]; }
        // (2) wpoly6_bwd in-place: grad_W → grad_r²
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_wp_bw];
          [e setBuffer:bR2 offset:0 atIndex:0]; [e setBuffer:bGW_or_Gr2 offset:0 atIndex:1];
          [e setBuffer:bH2 offset:0 atIndex:2]; [e setBuffer:bP6 offset:0 atIndex:3];
          [e setBuffer:bNr2 offset:0 atIndex:4];
          [e dispatchThreads:MTLSizeMake(n_r2, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [e endEncoding]; }
        // (3) dist_active_active_bwd: grad_r² → grad_active
        { id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
          [e setComputePipelineState:pso_d_bw];
          [e setBuffer:bA offset:0 atIndex:0]; [e setBuffer:bGW_or_Gr2 offset:0 atIndex:1];
          [e setBuffer:bGa offset:0 atIndex:2]; [e setBuffer:bNa offset:0 atIndex:3];
          [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
          [e endEncoding]; }
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fo = fopen(argv[8], "wb");
        fwrite([bGa contents], 1, a_bytes, fo); fclose(fo);
    }
    free(active); free(r2_saved); free(grad_d);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.C — density_constraint_grad_backward host op.
// CLI:
//   density_constraint_grad_bwd <n_active> <n_static> <h> <mass> <rho_rest>
//     <active.bin> <static.bin> <r2_aa.bin> <r2_as.bin>
//     <grad_grad_C.bin> <grad_denom_helper.bin> <grad_active_out.bin>
// ──────────────────────────────────────────────────────────────────────
int run_density_constraint_grad_bwd(int argc, char **argv) {
    if (argc != 14) {
        fprintf(stderr,
            "usage: sib_metal density_constraint_grad_bwd "
            "<n_active> <n_static> <h> <mass> <rho_rest> "
            "<active.bin> <static.bin> <r2_aa.bin> <r2_as.bin> "
            "<grad_grad_C.bin> <grad_denom_helper.bin> <grad_active_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h        = (float)atof(argv[4]);
    float mass     = (float)atof(argv[5]);
    float rho_rest = (float)atof(argv[6]);
    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));

    auto read_floats = ^(const char *path, size_t n_floats) {
        FILE *f = fopen(path, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
        float *buf = (float *)malloc(n_floats * sizeof(float));
        fread(buf, sizeof(float), n_floats, f); fclose(f); return buf;
    };
    float *active   = read_floats(argv[7], (size_t)n_active * 3);
    float *static_p = read_floats(argv[8], (size_t)n_static * 3);
    float *r2_aa    = read_floats(argv[9], (size_t)n_active * n_active);
    float *r2_as    = read_floats(argv[10], (size_t)n_active * n_static);
    float *g_gC     = read_floats(argv[11], (size_t)n_active * 3);
    float *g_dh     = read_floats(argv[12], (size_t)n_active);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso = make_pso(ctx, "density_constraint_grad_backward");

        size_t a_bytes = (size_t)n_active * 3 * sizeof(float);
        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:active length:a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bS = [ctx.device newBufferWithBytes:static_p
            length:(size_t)n_static * 3 * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2aa = [ctx.device newBufferWithBytes:r2_aa
            length:(size_t)n_active * n_active * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR2as = [ctx.device newBufferWithBytes:r2_as
            length:(size_t)n_active * n_static * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGgC = [ctx.device newBufferWithBytes:g_gC length:a_bytes
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGdh = [ctx.device newBufferWithBytes:g_dh
            length:(size_t)n_active * sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGa = [ctx.device newBufferWithLength:a_bytes
            options:MTLResourceStorageModeShared];
        memset([bGa contents], 0, a_bytes);

        id<MTLBuffer> bH = [ctx.device newBufferWithBytes:&h length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSp = [ctx.device newBufferWithBytes:&spiky_const length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNs = [ctx.device newBufferWithBytes:&n_static length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
        [e setComputePipelineState:pso];
        [e setBuffer:bA offset:0 atIndex:0];     [e setBuffer:bS offset:0 atIndex:1];
        [e setBuffer:bR2aa offset:0 atIndex:2];  [e setBuffer:bR2as offset:0 atIndex:3];
        [e setBuffer:bGgC offset:0 atIndex:4];   [e setBuffer:bGdh offset:0 atIndex:5];
        [e setBuffer:bGa offset:0 atIndex:6];
        [e setBuffer:bH offset:0 atIndex:7];     [e setBuffer:bSp offset:0 atIndex:8];
        [e setBuffer:bM offset:0 atIndex:9];     [e setBuffer:bR offset:0 atIndex:10];
        [e setBuffer:bNa offset:0 atIndex:11];   [e setBuffer:bNs offset:0 atIndex:12];
        [e dispatchThreads:MTLSizeMake(n_active, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
        [e endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fo = fopen(argv[13], "wb");
        fwrite([bGa contents], 1, a_bytes, fo); fclose(fo);
    }
    free(active); free(static_p); free(r2_aa); free(r2_as);
    free(g_gC); free(g_dh);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M9.D — solve_density_constraint isolated forward + backward ops.
//
// solve_density_constraint_fwd:
//   inputs: pos_pre, lambda_pre, density, grad_C, denom_helper,
//           constants (rho_rest, mass, alpha_inv_dt2)
//   outputs: pos_post, lambda_post (both modified in place)
//   This is the same kernel xpbd_step uses internally, exposed for
//   isolated testing.
//
// solve_density_constraint_bwd:
//   inputs: density, grad_C, denom_helper, lambda_pre (FROM FORWARD),
//           grad_pos_post, grad_lambda_post (gradient seeds)
//   outputs: grad_pos_pre, grad_lambda_pre, grad_density, grad_grad_C,
//            grad_denom_h, grad_rho_rest (per-particle, host sums)
// ──────────────────────────────────────────────────────────────────────
int run_solve_density_constraint_fwd(int argc, char **argv) {
    if (argc != 13) {
        fprintf(stderr,
            "usage: sib_metal solve_density_constraint_fwd "
            "<n_active> <rho_rest> <mass> <alpha_inv_dt2> "
            "<pos_pre.bin> <lambda_pre.bin> <density.bin> "
            "<grad_C.bin> <denom_helper.bin> "
            "<pos_post_out.bin> <lambda_post_out.bin>\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    float rho_rest = (float)atof(argv[3]);
    float mass     = (float)atof(argv[4]);
    float adt2     = (float)atof(argv[5]);

    auto rd = ^(const char *p, size_t n_floats) {
        FILE *f = fopen(p, "rb");
        if (!f) { fprintf(stderr, "open %s\n", p); exit(1); }
        float *b = (float *)malloc(n_floats * sizeof(float));
        fread(b, sizeof(float), n_floats, f); fclose(f); return b;
    };
    float *pos = rd(argv[6], (size_t)n * 3);
    float *lam = rd(argv[7], (size_t)n);
    float *dens = rd(argv[8], (size_t)n);
    float *gC = rd(argv[9], (size_t)n * 3);
    float *dh = rd(argv[10], (size_t)n);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso = make_pso(ctx, "solve_density_constraint");

        size_t pos_b = (size_t)n * 3 * sizeof(float);
        size_t s_b = (size_t)n * sizeof(float);
        id<MTLBuffer> bP = [ctx.device newBufferWithBytes:pos length:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bL = [ctx.device newBufferWithBytes:lam length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD = [ctx.device newBufferWithBytes:dens length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bG = [ctx.device newBufferWithBytes:gC length:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH = [ctx.device newBufferWithBytes:dh length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:&adt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
        [e setComputePipelineState:pso];
        [e setBuffer:bP offset:0 atIndex:0]; [e setBuffer:bL offset:0 atIndex:1];
        [e setBuffer:bD offset:0 atIndex:2]; [e setBuffer:bG offset:0 atIndex:3];
        [e setBuffer:bH offset:0 atIndex:4]; [e setBuffer:bR offset:0 atIndex:5];
        [e setBuffer:bM offset:0 atIndex:6]; [e setBuffer:bA offset:0 atIndex:7];
        [e setBuffer:bN offset:0 atIndex:8];
        [e dispatchThreads:MTLSizeMake(n,1,1)
            threadsPerThreadgroup:MTLSizeMake(64,1,1)];
        [e endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *o1 = fopen(argv[11], "wb"); fwrite([bP contents], 1, pos_b, o1); fclose(o1);
        FILE *o2 = fopen(argv[12], "wb"); fwrite([bL contents], 1, s_b, o2); fclose(o2);
    }
    free(pos); free(lam); free(dens); free(gC); free(dh);
    return 0;
}

int run_solve_density_constraint_bwd(int argc, char **argv) {
    if (argc != 18 && argc != 19) {
        fprintf(stderr,
            "usage: sib_metal solve_density_constraint_bwd "
            "<n_active> <rho_rest> <mass> <alpha_inv_dt2> "
            "<density.bin> <grad_C.bin> <denom_helper.bin> <lambda_pre.bin> "
            "<grad_pos_post.bin> <grad_lambda_post.bin> "
            "<grad_pos_pre_out.bin> <grad_lambda_pre_out.bin> "
            "<grad_density_out.bin> <grad_grad_C_out.bin> "
            "<grad_denom_helper_out.bin> <grad_rho_rest_per_particle.bin> "
            "[grad_alpha_inv_dt2_per_particle.bin]\n");
        return 1;
    }
    uint32_t n = (uint32_t)atoi(argv[2]);
    float rho_rest = (float)atof(argv[3]);
    float mass     = (float)atof(argv[4]);
    float adt2     = (float)atof(argv[5]);

    auto rd = ^(const char *p, size_t n_floats) {
        FILE *f = fopen(p, "rb");
        if (!f) { fprintf(stderr, "open %s\n", p); exit(1); }
        float *b = (float *)malloc(n_floats * sizeof(float));
        fread(b, sizeof(float), n_floats, f); fclose(f); return b;
    };
    float *dens = rd(argv[6], (size_t)n);
    float *gC   = rd(argv[7], (size_t)n * 3);
    float *dh   = rd(argv[8], (size_t)n);
    float *lam  = rd(argv[9], (size_t)n);
    float *gP   = rd(argv[10], (size_t)n * 3);
    float *gL   = rd(argv[11], (size_t)n);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso = make_pso(ctx, "solve_density_constraint_backward");

        size_t pos_b = (size_t)n * 3 * sizeof(float);
        size_t s_b = (size_t)n * sizeof(float);
        id<MTLBuffer> bD = [ctx.device newBufferWithBytes:dens length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGc = [ctx.device newBufferWithBytes:gC length:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH = [ctx.device newBufferWithBytes:dh length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bL = [ctx.device newBufferWithBytes:lam length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGp = [ctx.device newBufferWithBytes:gP length:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGl = [ctx.device newBufferWithBytes:gL length:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGpout = [ctx.device newBufferWithLength:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGlout = [ctx.device newBufferWithLength:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGdout = [ctx.device newBufferWithLength:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGgcout = [ctx.device newBufferWithLength:pos_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGhout = [ctx.device newBufferWithLength:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGrout = [ctx.device newBufferWithLength:s_b
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGaout = [ctx.device newBufferWithLength:s_b
            options:MTLResourceStorageModeShared];
        memset([bGpout contents], 0, pos_b);
        // Other outputs are written (not accumulated) so no need to zero.

        id<MTLBuffer> bR = [ctx.device newBufferWithBytes:&rho_rest length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:&mass length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bA = [ctx.device newBufferWithBytes:&adt2 length:sizeof(float)
            options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN = [ctx.device newBufferWithBytes:&n length:sizeof(uint32_t)
            options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> e = [cmd computeCommandEncoder];
        [e setComputePipelineState:pso];
        [e setBuffer:bD offset:0 atIndex:0];     [e setBuffer:bGc offset:0 atIndex:1];
        [e setBuffer:bH offset:0 atIndex:2];     [e setBuffer:bL offset:0 atIndex:3];
        [e setBuffer:bGp offset:0 atIndex:4];    [e setBuffer:bGl offset:0 atIndex:5];
        [e setBuffer:bGpout offset:0 atIndex:6]; [e setBuffer:bGlout offset:0 atIndex:7];
        [e setBuffer:bGdout offset:0 atIndex:8]; [e setBuffer:bGgcout offset:0 atIndex:9];
        [e setBuffer:bGhout offset:0 atIndex:10]; [e setBuffer:bGrout offset:0 atIndex:11];
        [e setBuffer:bGaout offset:0 atIndex:12];
        [e setBuffer:bR offset:0 atIndex:13];    [e setBuffer:bM offset:0 atIndex:14];
        [e setBuffer:bA offset:0 atIndex:15];    [e setBuffer:bN offset:0 atIndex:16];
        [e dispatchThreads:MTLSizeMake(n,1,1)
            threadsPerThreadgroup:MTLSizeMake(64,1,1)];
        [e endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *o1 = fopen(argv[12], "wb"); fwrite([bGpout contents], 1, pos_b, o1); fclose(o1);
        FILE *o2 = fopen(argv[13], "wb"); fwrite([bGlout contents], 1, s_b, o2); fclose(o2);
        FILE *o3 = fopen(argv[14], "wb"); fwrite([bGdout contents], 1, s_b, o3); fclose(o3);
        FILE *o4 = fopen(argv[15], "wb"); fwrite([bGgcout contents], 1, pos_b, o4); fclose(o4);
        FILE *o5 = fopen(argv[16], "wb"); fwrite([bGhout contents], 1, s_b, o5); fclose(o5);
        FILE *o6 = fopen(argv[17], "wb"); fwrite([bGrout contents], 1, s_b, o6); fclose(o6);
        if (argc == 19) {
            FILE *o7 = fopen(argv[18], "wb"); fwrite([bGaout contents], 1, s_b, o7); fclose(o7);
        }
    }
    free(dens); free(gC); free(dh); free(lam); free(gP); free(gL);
    return 0;
}

