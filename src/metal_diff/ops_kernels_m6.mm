// ops_kernels_m6.mm — atomic kernel ops (M6 substrate).
//
// Each op runs a single Metal kernel against pre-computed input buffers
// and writes its output. Useful as building blocks for FD-validation
// tests (Numpy reference vs. Metal output).
//
//   dist_active_static, dist_active_active: pairwise r² matrices
//   wpoly6_inplace:                         W = poly6_const · (h²-r²)³
//   rowsum_density:                         density[i] = mass · ΣW[i,j]
//   density_constraint_grad:                ∇C and Σ|∇W|² for solver

#include "metal_common.h"



int run_dist_active_static(int argc, char **argv) {
    if (argc != 7 && argc != 8) {
        fprintf(stderr,
            "usage: sib_metal dist_active_static "
            "<n_active> <n_static> <active.bin> <static.bin> <out.bin> [iters]\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    const char *path_active = argv[4];
    const char *path_static = argv[5];
    const char *path_out    = argv[6];
    int iters = (argc == 8) ? atoi(argv[7]) : 1;
    if (iters < 1) iters = 1;

    FILE *fa = fopen(path_active, "rb");
    FILE *fs = fopen(path_static, "rb");
    if (!fa || !fs) {
        fprintf(stderr, "cannot open input files\n");
        return 1;
    }
    float *active = (float *)malloc((size_t)n_active * 3 * sizeof(float));
    float *static_p = (float *)malloc((size_t)n_static * 3 * sizeof(float));
    if (fread(active, sizeof(float), n_active * 3, fa) != n_active * 3 ||
        fread(static_p, sizeof(float), n_static * 3, fs) != n_static * 3) {
        fprintf(stderr, "short read on input file\n");
        return 1;
    }
    fclose(fa);
    fclose(fs);

    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            fprintf(stderr, "no Metal device\n");
            return 1;
        }
        NSError *err = nil;
        id<MTLLibrary> lib = [device
            newLibraryWithSource:load_shader_source()
                         options:nil
                           error:&err];
        if (!lib) {
            fprintf(stderr, "shader compile: %s\n",
                [[err localizedDescription] UTF8String]);
            return 1;
        }
        id<MTLFunction> fn = [lib newFunctionWithName:@"dist_active_static"];
        id<MTLComputePipelineState> pso =
            [device newComputePipelineStateWithFunction:fn error:&err];
        if (!pso) {
            fprintf(stderr, "pipeline: %s\n",
                [[err localizedDescription] UTF8String]);
            return 1;
        }
        id<MTLCommandQueue> queue = [device newCommandQueue];

        id<MTLBuffer> bufActive = [device
            newBufferWithBytes:active
                        length:(size_t)n_active * 3 * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufStatic = [device
            newBufferWithBytes:static_p
                        length:(size_t)n_static * 3 * sizeof(float)
                       options:MTLResourceStorageModeShared];
        size_t out_bytes = (size_t)n_active * n_static * sizeof(float);
        id<MTLBuffer> bufDist = [device
            newBufferWithLength:out_bytes
                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNa = [device
            newBufferWithBytes:&n_active
                        length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNs = [device
            newBufferWithBytes:&n_static
                        length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];

        // 16x16 threadgroup: a reasonable default for 2D grids on Apple GPUs.
        // Total threads = 256 per group, fits Apple's 1024 max comfortably.
        MTLSize threadgroup = MTLSizeMake(16, 16, 1);
        MTLSize grid = MTLSizeMake(n_active, n_static, 1);

        // Steady-state benchmark: run `iters` times, time only the
        // dispatch+wait, skip device/library/buffer setup.
        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int it = 0; it < iters; it++) {
            id<MTLCommandBuffer> cmd = [queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufActive offset:0 atIndex:0];
            [enc setBuffer:bufStatic offset:0 atIndex:1];
            [enc setBuffer:bufDist  offset:0 atIndex:2];
            [enc setBuffer:bufNa    offset:0 atIndex:3];
            [enc setBuffer:bufNs    offset:0 atIndex:4];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_iter_ms = (t1 - t0) * 1000.0 / iters;
        fprintf(stderr, "kernel: %d iters, %.3f ms/iter (steady state)\n",
                iters, per_iter_ms);

        FILE *fo = fopen(path_out, "wb");
        if (!fo) { fprintf(stderr, "cannot open output\n"); return 1; }
        fwrite([bufDist contents], 1, out_bytes, fo);
        fclose(fo);
    }
    free(active);
    free(static_p);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M6.1 — wpoly6_inplace
// ──────────────────────────────────────────────────────────────────────
int run_wpoly6_inplace(int argc, char **argv) {
    if (argc != 5 && argc != 6) {
        fprintf(stderr,
            "usage: sib_metal wpoly6_inplace "
            "<n_total> <h> <inout.bin> [iters]\n");
        return 1;
    }
    uint32_t n_total = (uint32_t)atoi(argv[2]);
    float h = (float)atof(argv[3]);
    const char *path_inout = argv[4];
    int iters = (argc == 6) ? atoi(argv[5]) : 1;
    if (iters < 1) iters = 1;

    float h2 = h * h;
    float poly6_const = 315.0f / (64.0f * (float)M_PI * powf(h, 9.0f));

    FILE *f = fopen(path_inout, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path_inout); return 1; }
    float *data = (float *)malloc((size_t)n_total * sizeof(float));
    if (fread(data, sizeof(float), n_total, f) != n_total) {
        fprintf(stderr, "short read on %s\n", path_inout); return 1;
    }
    fclose(f);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso = make_pso(ctx, "wpoly6_inplace");

        id<MTLBuffer> bufData = [ctx.device
            newBufferWithBytes:data
                        length:(size_t)n_total * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufH2 = [ctx.device
            newBufferWithBytes:&h2 length:sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufC = [ctx.device
            newBufferWithBytes:&poly6_const length:sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufN = [ctx.device
            newBufferWithBytes:&n_total length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];

        MTLSize threadgroup = MTLSizeMake(256, 1, 1);
        MTLSize grid = MTLSizeMake(n_total, 1, 1);

        // For steady-state timing of an in-place kernel, the per-iter
        // memcpy needed to restore the input would inflate the number.
        // We measure ONLY kernel dispatch+wait. After iter 1 the buffer
        // holds the output (W values, mostly zeros) — running the
        // kernel again gives wrong outputs but the same compute cost,
        // which is what we want to measure. Final on-disk output is
        // restored to the iter-1 result before write-back.
        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int it = 0; it < iters; it++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufData offset:0 atIndex:0];
            [enc setBuffer:bufH2   offset:0 atIndex:1];
            [enc setBuffer:bufC    offset:0 atIndex:2];
            [enc setBuffer:bufN    offset:0 atIndex:3];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_iter_ms = (t1 - t0) * 1000.0 / iters;
        fprintf(stderr, "kernel: %d iters, %.3f ms/iter (steady state)\n",
                iters, per_iter_ms);

        // Restore the iter-1 result for write-back: re-upload input,
        // run kernel once more.
        if (iters > 1) {
            memcpy([bufData contents], data, (size_t)n_total * sizeof(float));
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufData offset:0 atIndex:0];
            [enc setBuffer:bufH2   offset:0 atIndex:1];
            [enc setBuffer:bufC    offset:0 atIndex:2];
            [enc setBuffer:bufN    offset:0 atIndex:3];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }

        FILE *fo = fopen(path_inout, "wb");
        if (!fo) { fprintf(stderr, "cannot open output\n"); return 1; }
        fwrite([bufData contents], 1, (size_t)n_total * sizeof(float), fo);
        fclose(fo);
    }
    free(data);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M6.2 — rowsum_density
// ──────────────────────────────────────────────────────────────────────
int run_rowsum_density(int argc, char **argv) {
    if (argc != 7 && argc != 8) {
        fprintf(stderr,
            "usage: sib_metal rowsum_density "
            "<n_rows> <n_cols> <mass> <W.bin> <density.bin> [iters]\n");
        return 1;
    }
    uint32_t n_rows = (uint32_t)atoi(argv[2]);
    uint32_t n_cols = (uint32_t)atoi(argv[3]);
    float    mass   = (float)atof(argv[4]);
    const char *path_W       = argv[5];
    const char *path_density = argv[6];
    int iters = (argc == 8) ? atoi(argv[7]) : 1;
    if (iters < 1) iters = 1;

    size_t w_count = (size_t)n_rows * n_cols;
    FILE *f = fopen(path_W, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path_W); return 1; }
    float *W = (float *)malloc(w_count * sizeof(float));
    if (fread(W, sizeof(float), w_count, f) != w_count) {
        fprintf(stderr, "short read on %s\n", path_W); return 1;
    }
    fclose(f);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso = make_pso(ctx, "rowsum_density");

        id<MTLBuffer> bufW = [ctx.device
            newBufferWithBytes:W
                        length:w_count * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufRho = [ctx.device
            newBufferWithLength:(size_t)n_rows * sizeof(float)
                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufMass = [ctx.device
            newBufferWithBytes:&mass length:sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNc = [ctx.device
            newBufferWithBytes:&n_cols length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNr = [ctx.device
            newBufferWithBytes:&n_rows length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];

        // 2D dispatch: X axis = threads within threadgroup (256),
        // Y axis = row index. Threadgroup size 256x1 means each row
        // gets its own threadgroup of 256 threads doing tree reduction.
        MTLSize threadgroup = MTLSizeMake(256, 1, 1);
        MTLSize grid = MTLSizeMake(256, n_rows, 1);

        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int it = 0; it < iters; it++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufW    offset:0 atIndex:0];
            [enc setBuffer:bufRho  offset:0 atIndex:1];
            [enc setBuffer:bufMass offset:0 atIndex:2];
            [enc setBuffer:bufNc   offset:0 atIndex:3];
            [enc setBuffer:bufNr   offset:0 atIndex:4];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_iter_ms = (t1 - t0) * 1000.0 / iters;
        fprintf(stderr, "kernel: %d iters, %.3f ms/iter (steady state)\n",
                iters, per_iter_ms);

        FILE *fo = fopen(path_density, "wb");
        if (!fo) { fprintf(stderr, "cannot open output\n"); return 1; }
        fwrite([bufRho contents], 1, (size_t)n_rows * sizeof(float), fo);
        fclose(fo);
    }
    free(W);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M6.3 — dist_active_active
// ──────────────────────────────────────────────────────────────────────
int run_dist_active_active(int argc, char **argv) {
    if (argc != 5 && argc != 6) {
        fprintf(stderr,
            "usage: sib_metal dist_active_active "
            "<n_active> <active.bin> <out.bin> [iters]\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    const char *path_active = argv[3];
    const char *path_out    = argv[4];
    int iters = (argc == 6) ? atoi(argv[5]) : 1;
    if (iters < 1) iters = 1;

    FILE *f = fopen(path_active, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path_active); return 1; }
    float *active = (float *)malloc((size_t)n_active * 3 * sizeof(float));
    if (fread(active, sizeof(float), n_active * 3, f) != n_active * 3) {
        fprintf(stderr, "short read on %s\n", path_active); return 1;
    }
    fclose(f);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso =
            make_pso(ctx, "dist_active_active");

        id<MTLBuffer> bufActive = [ctx.device
            newBufferWithBytes:active
                        length:(size_t)n_active * 3 * sizeof(float)
                       options:MTLResourceStorageModeShared];
        size_t out_bytes = (size_t)n_active * n_active * sizeof(float);
        id<MTLBuffer> bufDist = [ctx.device
            newBufferWithLength:out_bytes
                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufN = [ctx.device
            newBufferWithBytes:&n_active length:sizeof(uint32_t)
                       options:MTLResourceStorageModeShared];

        MTLSize threadgroup = MTLSizeMake(16, 16, 1);
        MTLSize grid = MTLSizeMake(n_active, n_active, 1);

        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int it = 0; it < iters; it++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufActive offset:0 atIndex:0];
            [enc setBuffer:bufDist   offset:0 atIndex:1];
            [enc setBuffer:bufN      offset:0 atIndex:2];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_iter_ms = (t1 - t0) * 1000.0 / iters;
        fprintf(stderr, "kernel: %d iters, %.3f ms/iter (steady state)\n",
                iters, per_iter_ms);

        FILE *fo = fopen(path_out, "wb");
        if (!fo) { fprintf(stderr, "cannot open output\n"); return 1; }
        fwrite([bufDist contents], 1, out_bytes, fo);
        fclose(fo);
    }
    free(active);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// M6.4 — density_constraint_grad (also outputs denom_helper)
// ──────────────────────────────────────────────────────────────────────
int run_density_constraint_grad(int argc, char **argv) {
    if (argc != 13 && argc != 14) {
        fprintf(stderr,
            "usage: sib_metal density_constraint_grad "
            "<n_active> <n_static> <h> <mass> <rho_rest> "
            "<active.bin> <static.bin> <r2_aa.bin> <r2_as.bin> "
            "<gradC.bin> <denom_helper.bin> [iters]\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h          = (float)atof(argv[4]);
    float mass       = (float)atof(argv[5]);
    float rho_rest   = (float)atof(argv[6]);
    const char *path_active = argv[7];
    const char *path_static = argv[8];
    const char *path_r2_aa  = argv[9];
    const char *path_r2_as  = argv[10];
    const char *path_gradC  = argv[11];
    const char *path_denomH = argv[12];
    int iters = (argc == 14) ? atoi(argv[13]) : 1;
    if (iters < 1) iters = 1;

    float spiky_const = -45.0f / ((float)M_PI * powf(h, 6.0f));

    // Load all four input buffers.
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

    float *active   = read_floats(path_active, (size_t)n_active * 3);
    float *static_p = read_floats(path_static, (size_t)n_static * 3);
    float *r2_aa    = read_floats(path_r2_aa,  (size_t)n_active * n_active);
    float *r2_as    = read_floats(path_r2_as,  (size_t)n_active * n_static);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLComputePipelineState> pso =
            make_pso(ctx, "density_constraint_grad");

        id<MTLBuffer> bufActive = [ctx.device
            newBufferWithBytes:active
                        length:(size_t)n_active * 3 * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufStatic = [ctx.device
            newBufferWithBytes:static_p
                        length:(size_t)n_static * 3 * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufR2aa = [ctx.device
            newBufferWithBytes:r2_aa
                        length:(size_t)n_active * n_active * sizeof(float)
                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufR2as = [ctx.device
            newBufferWithBytes:r2_as
                        length:(size_t)n_active * n_static * sizeof(float)
                       options:MTLResourceStorageModeShared];
        // packed_float3 = 12 bytes per element.
        size_t out_bytes = (size_t)n_active * 3 * sizeof(float);
        size_t denom_bytes = (size_t)n_active * sizeof(float);
        id<MTLBuffer> bufGradC = [ctx.device
            newBufferWithLength:out_bytes
                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufDenomH = [ctx.device
            newBufferWithLength:denom_bytes
                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufH        = [ctx.device newBufferWithBytes:&h
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufSpiky    = [ctx.device newBufferWithBytes:&spiky_const
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufMass     = [ctx.device newBufferWithBytes:&mass
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufRho      = [ctx.device newBufferWithBytes:&rho_rest
            length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNa       = [ctx.device newBufferWithBytes:&n_active
            length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bufNs       = [ctx.device newBufferWithBytes:&n_static
            length:sizeof(uint32_t) options:MTLResourceStorageModeShared];

        MTLSize threadgroup = MTLSizeMake(256, 1, 1);
        MTLSize grid = MTLSizeMake(256, n_active, 1);

        CFAbsoluteTime t0 = CFAbsoluteTimeGetCurrent();
        for (int it = 0; it < iters; it++) {
            id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bufActive  offset:0 atIndex:0];
            [enc setBuffer:bufStatic  offset:0 atIndex:1];
            [enc setBuffer:bufR2aa    offset:0 atIndex:2];
            [enc setBuffer:bufR2as    offset:0 atIndex:3];
            [enc setBuffer:bufGradC   offset:0 atIndex:4];
            [enc setBuffer:bufDenomH  offset:0 atIndex:5];
            [enc setBuffer:bufH       offset:0 atIndex:6];
            [enc setBuffer:bufSpiky   offset:0 atIndex:7];
            [enc setBuffer:bufMass    offset:0 atIndex:8];
            [enc setBuffer:bufRho     offset:0 atIndex:9];
            [enc setBuffer:bufNa      offset:0 atIndex:10];
            [enc setBuffer:bufNs      offset:0 atIndex:11];
            [enc dispatchThreads:grid threadsPerThreadgroup:threadgroup];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
        CFAbsoluteTime t1 = CFAbsoluteTimeGetCurrent();
        double per_iter_ms = (t1 - t0) * 1000.0 / iters;
        fprintf(stderr, "kernel: %d iters, %.3f ms/iter (steady state)\n",
                iters, per_iter_ms);

        FILE *fo = fopen(path_gradC, "wb");
        if (!fo) { fprintf(stderr, "cannot open output\n"); return 1; }
        fwrite([bufGradC contents], 1, out_bytes, fo);
        fclose(fo);
        FILE *fd = fopen(path_denomH, "wb");
        if (!fd) { fprintf(stderr, "cannot open denom output\n"); return 1; }
        fwrite([bufDenomH contents], 1, denom_bytes, fd);
        fclose(fd);
    }
    free(active); free(static_p); free(r2_aa); free(r2_as);
    return 0;
}

