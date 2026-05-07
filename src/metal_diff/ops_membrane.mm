// ops_membrane.mm — M10 membrane test ops (forward + backward).
//
// Each op is a thin host-side wrapper that loads its inputs from
// /tmp/*.bin, dispatches the kernel once, writes outputs. Used by
// the FD-validation tests in test_membrane_*.py.
//
// CLI conventions match ops_pair_spring.mm: positional args, fixed
// /tmp paths for outputs. Keeps the test harness simple.

#include "metal_common.h"

// ──────────────────────────────────────────────────────────────────────
// accumulate_membrane_fwd — runs the forward kernel and writes mem_corr.
//
// Args (12 total incl. op name):
//   <n_active> <n_elastic> <n_membranes> <r0>
//   <pos.bin>           [n_active*3 fp32]
//   <membranes.bin>     [n_membranes*3 int32]
//   <pmem_idx.bin>      [n_elastic*7 int32]
//   <mem_corr_init.bin> [n_active*3 fp32]   pre-existing accumulator value
// Output: /tmp/membrane_corr_out.bin  ([n_active*3] fp32)
// ──────────────────────────────────────────────────────────────────────
int run_accumulate_membrane_fwd(int argc, char **argv) {
    if (argc != 10) {
        fprintf(stderr, "usage: sib_metal accumulate_membrane_fwd "
                "<n_active> <n_elastic> <n_membranes> <r0> "
                "<pos.bin> <membranes.bin> <pmem_idx.bin> "
                "<mem_corr_init.bin>\n");
        return 1;
    }
    uint32_t n_active    = (uint32_t)atoi(argv[2]);
    uint32_t n_elastic   = (uint32_t)atoi(argv[3]);
    uint32_t n_membranes = (uint32_t)atoi(argv[4]);
    float    r0          = (float)atof(argv[5]);
    const char *p_pos    = argv[6];
    const char *p_mem    = argv[7];
    const char *p_pmi    = argv[8];
    const char *p_init   = argv[9];

    auto rd = ^(const char *p, size_t bytes) {
        FILE *f = fopen(p, "rb"); if (!f) { fprintf(stderr,"open %s\n",p); exit(1); }
        void *b = malloc(bytes);
        if (fread(b, 1, bytes, f) != bytes) { fprintf(stderr,"read %s\n",p); exit(1); }
        fclose(f); return b;
    };
    size_t pos_b  = (size_t)n_active * 3 * sizeof(float);
    size_t mem_b  = (size_t)n_membranes * 3 * sizeof(int32_t);
    size_t pmi_b  = (size_t)n_elastic * 7 * sizeof(int32_t);
    float   *pos     = (float   *)rd(p_pos,  pos_b);
    int32_t *mem     = (int32_t *)rd(p_mem,  mem_b);
    int32_t *pmi     = (int32_t *)rd(p_pmi,  pmi_b);
    float   *init    = (float   *)rd(p_init, pos_b);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLBuffer> bP   = [ctx.device newBufferWithBytes:pos length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM   = [ctx.device newBufferWithBytes:mem length:mem_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bPI  = [ctx.device newBufferWithBytes:pmi length:pmi_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bMC  = [ctx.device newBufferWithBytes:init length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa  = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNe  = [ctx.device newBufferWithBytes:&n_elastic length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR0  = [ctx.device newBufferWithBytes:&r0 length:sizeof(float) options:MTLResourceStorageModeShared];

        id<MTLComputePipelineState> pso = make_pso(ctx, "accumulate_membrane_correction");
        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pso];
        [enc setBuffer:bP  offset:0 atIndex:0];
        [enc setBuffer:bM  offset:0 atIndex:1];
        [enc setBuffer:bPI offset:0 atIndex:2];
        [enc setBuffer:bMC offset:0 atIndex:3];
        [enc setBuffer:bNa offset:0 atIndex:4];
        [enc setBuffer:bNe offset:0 atIndex:5];
        [enc setBuffer:bR0 offset:0 atIndex:6];
        [enc dispatchThreads:MTLSizeMake(n_active,1,1)
            threadsPerThreadgroup:MTLSizeMake(64,1,1)];
        [enc endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fo = fopen("/tmp/membrane_corr_out.bin","wb");
        fwrite([bMC contents], 1, pos_b, fo); fclose(fo);
    }
    free(pos); free(mem); free(pmi); free(init);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// accumulate_membrane_bwd — runs the backward kernel and writes grad_pos.
//
// Args (11 total incl. op name):
//   <n_active> <n_elastic> <n_membranes> <r0>
//   <pos.bin>          [n_active*3 fp32]
//   <membranes.bin>    [n_membranes*3 int32]
//   <pmem_idx.bin>     [n_elastic*7 int32]
//   <grad_corr.bin>    [n_active*3 fp32] upstream grad of mem_corr
// Output: /tmp/membrane_grad_pos.bin  ([n_active*3] fp32)
// ──────────────────────────────────────────────────────────────────────
int run_accumulate_membrane_bwd(int argc, char **argv) {
    if (argc != 10) {
        fprintf(stderr, "usage: sib_metal accumulate_membrane_bwd "
                "<n_active> <n_elastic> <n_membranes> <r0> "
                "<pos.bin> <membranes.bin> <pmem_idx.bin> "
                "<grad_corr.bin>\n");
        return 1;
    }
    uint32_t n_active    = (uint32_t)atoi(argv[2]);
    uint32_t n_elastic   = (uint32_t)atoi(argv[3]);
    uint32_t n_membranes = (uint32_t)atoi(argv[4]);
    float    r0          = (float)atof(argv[5]);
    const char *p_pos    = argv[6];
    const char *p_mem    = argv[7];
    const char *p_pmi    = argv[8];
    const char *p_grad   = argv[9];

    auto rd = ^(const char *p, size_t bytes) {
        FILE *f = fopen(p, "rb"); if (!f) { fprintf(stderr,"open %s\n",p); exit(1); }
        void *b = malloc(bytes);
        if (fread(b, 1, bytes, f) != bytes) { fprintf(stderr,"read %s\n",p); exit(1); }
        fclose(f); return b;
    };
    size_t pos_b  = (size_t)n_active * 3 * sizeof(float);
    size_t mem_b  = (size_t)n_membranes * 3 * sizeof(int32_t);
    size_t pmi_b  = (size_t)n_elastic * 7 * sizeof(int32_t);
    float   *pos     = (float   *)rd(p_pos,  pos_b);
    int32_t *mem     = (int32_t *)rd(p_mem,  mem_b);
    int32_t *pmi     = (int32_t *)rd(p_pmi,  pmi_b);
    float   *grad    = (float   *)rd(p_grad, pos_b);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLBuffer> bP   = [ctx.device newBufferWithBytes:pos length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM   = [ctx.device newBufferWithBytes:mem length:mem_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bPI  = [ctx.device newBufferWithBytes:pmi length:pmi_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGc  = [ctx.device newBufferWithBytes:grad length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGp  = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        memset([bGp contents], 0, pos_b);
        id<MTLBuffer> bNa  = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNe  = [ctx.device newBufferWithBytes:&n_elastic length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bR0  = [ctx.device newBufferWithBytes:&r0 length:sizeof(float) options:MTLResourceStorageModeShared];

        id<MTLComputePipelineState> pso = make_pso(ctx, "accumulate_membrane_correction_backward");
        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pso];
        [enc setBuffer:bP  offset:0 atIndex:0];
        [enc setBuffer:bM  offset:0 atIndex:1];
        [enc setBuffer:bPI offset:0 atIndex:2];
        [enc setBuffer:bGc offset:0 atIndex:3];
        [enc setBuffer:bGp offset:0 atIndex:4];
        [enc setBuffer:bNa offset:0 atIndex:5];
        [enc setBuffer:bNe offset:0 atIndex:6];
        [enc setBuffer:bR0 offset:0 atIndex:7];
        [enc dispatchThreads:MTLSizeMake(n_active,1,1)
            threadsPerThreadgroup:MTLSizeMake(64,1,1)];
        [enc endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fo = fopen("/tmp/membrane_grad_pos.bin","wb");
        fwrite([bGp contents], 1, pos_b, fo); fclose(fo);
    }
    free(pos); free(mem); free(pmi); free(grad);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// apply_membrane_fwd / _bwd — trivial forward/backward of pos += mem_corr.
// Mostly here so the dispatch table is symmetric with the other M10 ops.
// ──────────────────────────────────────────────────────────────────────
int run_apply_membrane_fwd(int argc, char **argv) {
    if (argc != 5) {
        fprintf(stderr, "usage: sib_metal apply_membrane_fwd "
                "<n_active> <pos.bin> <mem_corr.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    auto rd = ^(const char *p, size_t bytes) {
        FILE *f = fopen(p, "rb"); if (!f) exit(1);
        void *b = malloc(bytes);
        fread(b, 1, bytes, f); fclose(f); return b;
    };
    size_t pos_b = (size_t)n_active * 3 * sizeof(float);
    float *pos = (float*)rd(argv[3], pos_b);
    float *mc  = (float*)rd(argv[4], pos_b);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLBuffer> bP = [ctx.device newBufferWithBytes:pos length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM = [ctx.device newBufferWithBytes:mc  length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLComputePipelineState> pso = make_pso(ctx, "apply_membrane_correction");
        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pso];
        [enc setBuffer:bP offset:0 atIndex:0];
        [enc setBuffer:bM offset:0 atIndex:1];
        [enc setBuffer:bN offset:0 atIndex:2];
        [enc dispatchThreads:MTLSizeMake(n_active,1,1)
            threadsPerThreadgroup:MTLSizeMake(64,1,1)];
        [enc endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];
        FILE *fo = fopen("/tmp/membrane_apply_out.bin","wb");
        fwrite([bP contents], 1, pos_b, fo); fclose(fo);
    }
    free(pos); free(mc);
    return 0;
}

int run_apply_membrane_bwd(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: sib_metal apply_membrane_bwd "
                "<n_active> <grad_pos_post.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    auto rd = ^(const char *p, size_t bytes) {
        FILE *f = fopen(p, "rb"); if (!f) exit(1);
        void *b = malloc(bytes);
        fread(b, 1, bytes, f); fclose(f); return b;
    };
    size_t pos_b = (size_t)n_active * 3 * sizeof(float);
    float *gp = (float*)rd(argv[3], pos_b);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLBuffer> bGpost = [ctx.device newBufferWithBytes:gp length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGpre  = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGmc   = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        memset([bGpre contents], 0, pos_b);
        memset([bGmc contents], 0, pos_b);
        id<MTLBuffer> bN = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLComputePipelineState> pso = make_pso(ctx, "apply_membrane_correction_backward");
        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pso];
        [enc setBuffer:bGpost offset:0 atIndex:0];
        [enc setBuffer:bGpre  offset:0 atIndex:1];
        [enc setBuffer:bGmc   offset:0 atIndex:2];
        [enc setBuffer:bN     offset:0 atIndex:3];
        [enc dispatchThreads:MTLSizeMake(n_active,1,1)
            threadsPerThreadgroup:MTLSizeMake(64,1,1)];
        [enc endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];
        FILE *fp = fopen("/tmp/membrane_apply_grad_pos.bin","wb");
        fwrite([bGpre contents], 1, pos_b, fp); fclose(fp);
        FILE *fm = fopen("/tmp/membrane_apply_grad_mc.bin","wb");
        fwrite([bGmc contents], 1, pos_b, fm); fclose(fm);
    }
    free(gp);
    return 0;
}
