// ops_pair_spring.mm — pair forces (visc + surf tension) and spring
// bonds (Hooke), each with a forward and backward op for FD validation.

#include "metal_common.h"

// ──────────────────────────────────────────────────────────────────────
// pair_forces_fwd / _bwd — standalone test ops for pair_forces_grid
// kernel and its backward. Used by test_pair_forces_bwd.py for
// finite-diff validation.
//
// Args (forward, 16 total incl. op name):
//   <n_active> <n_static> <h> <mass> <sim_scale> <visc_pair_coef>
//   <pos_active.bin> <vel_active.bin> <sorted_static.bin>
//   <cell_start.bin> <density.bin>
//   <grid_dim_x> <grid_dim_y> <grid_dim_z>
//   <grid_origin.bin>
// Output: /tmp/pair_ext_accel.bin  ([n_active*3] fp32)
//
// Args (backward, 17 total): same as forward, plus
//   <grad_ext_accel.bin>
// Output: /tmp/pair_grad_pos.bin, /tmp/pair_grad_vel.bin
// ──────────────────────────────────────────────────────────────────────
int run_pair_forces_fwd(int argc, char **argv) {
    if (argc != 16) {
        fprintf(stderr, "usage: sib_metal pair_forces_fwd "
                "<n_active> <n_static> <h> <mass> <sim_scale> <visc_pair_coef> "
                "<pos_active.bin> <vel_active.bin> <sorted_static.bin> "
                "<cell_start.bin> <density.bin> "
                "<grid_dim_x> <grid_dim_y> <grid_dim_z> <grid_origin.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h           = (float)atof(argv[4]);
    float mass        = (float)atof(argv[5]);
    float sim_scale   = (float)atof(argv[6]);
    float visc_pair_coef = (float)atof(argv[7]);
    const char *path_pa = argv[8];
    const char *path_va = argv[9];
    const char *path_ss = argv[10];
    const char *path_cs = argv[11];
    const char *path_d  = argv[12];
    int grid_dim_x = atoi(argv[13]);
    int grid_dim_y = atoi(argv[14]);
    int grid_dim_z = atoi(argv[15]);

    auto read_bytes = ^(const char *path, size_t n) {
        FILE *f = fopen(path, "rb"); if (!f) { fprintf(stderr,"open %s\n",path); exit(1); }
        void *buf = malloc(n);
        if (fread(buf, 1, n, f) != n) { fprintf(stderr,"read %s\n",path); exit(1); }
        fclose(f); return buf;
    };
    float *pos_a = (float*)read_bytes(path_pa, (size_t)n_active*3*sizeof(float));
    float *vel_a = (float*)read_bytes(path_va, (size_t)n_active*3*sizeof(float));
    float *sort_s = (float*)read_bytes(path_ss, (size_t)n_static*3*sizeof(float));
    int n_cells = grid_dim_x * grid_dim_y * grid_dim_z;
    int *cs = (int*)read_bytes(path_cs, (size_t)(n_cells+1)*sizeof(int));
    float *dens = (float*)read_bytes(path_d, (size_t)n_active*sizeof(float));

    float h2 = h*h;
    double h_scaled = (double)h * (double)sim_scale;
    double h_s6 = pow(h_scaled, 6.0);
    double h_s9 = pow(h_scaled, 9.0);
    double divgradWvisco = 45.0 / (M_PI * h_s6);
    float visc_amp = (float)(1.5 * (double)mass * divgradWvisco * pow((double)sim_scale, 3.0));
    double wpoly6_si = 315.0 / (64.0 * M_PI * h_s9);
    float surf_amp = (float)(-1.7e-9 * (double)mass * wpoly6_si * (double)sim_scale / (double)mass);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        size_t pos_a_b = (size_t)n_active*3*sizeof(float);
        id<MTLBuffer> bPa = [ctx.device newBufferWithBytes:pos_a length:pos_a_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bVa = [ctx.device newBufferWithBytes:vel_a length:pos_a_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSs = [ctx.device newBufferWithBytes:sort_s length:(size_t)n_static*3*sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bCs = [ctx.device newBufferWithBytes:cs length:(size_t)(n_cells+1)*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD  = [ctx.device newBufferWithBytes:dens length:(size_t)n_active*sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bEa = [ctx.device newBufferWithLength:pos_a_b options:MTLResourceStorageModeShared];

        id<MTLBuffer> bH  = [ctx.device newBufferWithBytes:&h length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM  = [ctx.device newBufferWithBytes:&mass length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSS = [ctx.device newBufferWithBytes:&sim_scale length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bVp = [ctx.device newBufferWithBytes:&visc_pair_coef length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bVa2= [ctx.device newBufferWithBytes:&visc_amp length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSa = [ctx.device newBufferWithBytes:&surf_amp length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN  = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];

        int grid_dim_packed[3] = {grid_dim_x, grid_dim_y, grid_dim_z};
        id<MTLBuffer> bGd = [ctx.device newBufferWithBytes:grid_dim_packed length:3*sizeof(int) options:MTLResourceStorageModeShared];
        // Test op reads grid_origin from a fixed /tmp path so the CLI
        // arglist stays manageable.
        float gox[3] = {0,0,0};
        FILE *gf = fopen("/tmp/pair_grid_origin.bin","rb");
        if (gf) { fread(gox, sizeof(float), 3, gf); fclose(gf); }
        id<MTLBuffer> bGo = [ctx.device newBufferWithBytes:gox length:3*sizeof(float) options:MTLResourceStorageModeShared];

        id<MTLComputePipelineState> pso = make_pso(ctx, "pair_forces_grid");
        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pso];
        [enc setBuffer:bPa offset:0 atIndex:0];
        [enc setBuffer:bVa offset:0 atIndex:1];
        [enc setBuffer:bSs offset:0 atIndex:2];
        [enc setBuffer:bCs offset:0 atIndex:3];
        [enc setBuffer:bD  offset:0 atIndex:4];
        [enc setBuffer:bEa offset:0 atIndex:5];
        [enc setBuffer:bH  offset:0 atIndex:6];
        [enc setBuffer:bH2 offset:0 atIndex:7];
        [enc setBuffer:bM  offset:0 atIndex:8];
        [enc setBuffer:bSS offset:0 atIndex:9];
        [enc setBuffer:bVp offset:0 atIndex:10];
        [enc setBuffer:bVa2 offset:0 atIndex:11];
        [enc setBuffer:bSa offset:0 atIndex:12];
        [enc setBuffer:bN  offset:0 atIndex:13];
        [enc setBuffer:bGd offset:0 atIndex:14];
        [enc setBuffer:bGo offset:0 atIndex:15];
        [enc dispatchThreads:MTLSizeMake(32, n_active, 1)
            threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
        [enc endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fo = fopen("/tmp/pair_ext_accel.bin","wb");
        fwrite([bEa contents], 1, pos_a_b, fo); fclose(fo);
    }
    free(pos_a); free(vel_a); free(sort_s); free(cs); free(dens);
    return 0;
}

int run_pair_forces_bwd(int argc, char **argv) {
    if (argc != 17) {
        fprintf(stderr, "usage: sib_metal pair_forces_bwd "
                "<n_active> <n_static> <h> <mass> <sim_scale> <visc_pair_coef> "
                "<pos_active.bin> <vel_active.bin> <sorted_static.bin> "
                "<cell_start.bin> <density.bin> "
                "<grid_dim_x> <grid_dim_y> <grid_dim_z> <grid_origin.bin> "
                "<grad_ext_accel.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_static = (uint32_t)atoi(argv[3]);
    float h           = (float)atof(argv[4]);
    float mass        = (float)atof(argv[5]);
    float sim_scale   = (float)atof(argv[6]);
    float visc_pair_coef = (float)atof(argv[7]);
    const char *path_pa = argv[8];
    const char *path_va = argv[9];
    const char *path_ss = argv[10];
    const char *path_cs = argv[11];
    const char *path_d  = argv[12];
    int grid_dim_x = atoi(argv[13]);
    int grid_dim_y = atoi(argv[14]);
    int grid_dim_z = atoi(argv[15]);
    const char *path_ge = argv[16];

    auto read_bytes = ^(const char *path, size_t n) {
        FILE *f = fopen(path, "rb"); if (!f) { fprintf(stderr,"open %s\n",path); exit(1); }
        void *buf = malloc(n);
        if (fread(buf, 1, n, f) != n) { fprintf(stderr,"read %s\n",path); exit(1); }
        fclose(f); return buf;
    };
    size_t pos_a_b = (size_t)n_active*3*sizeof(float);
    float *pos_a = (float*)read_bytes(path_pa, pos_a_b);
    float *vel_a = (float*)read_bytes(path_va, pos_a_b);
    float *sort_s = (float*)read_bytes(path_ss, (size_t)n_static*3*sizeof(float));
    int n_cells = grid_dim_x * grid_dim_y * grid_dim_z;
    int *cs = (int*)read_bytes(path_cs, (size_t)(n_cells+1)*sizeof(int));
    float *dens = (float*)read_bytes(path_d, (size_t)n_active*sizeof(float));
    float *gea = (float*)read_bytes(path_ge, pos_a_b);

    float h2 = h*h;
    double h_scaled = (double)h * (double)sim_scale;
    double h_s6 = pow(h_scaled, 6.0);
    double h_s9 = pow(h_scaled, 9.0);
    double divgradWvisco = 45.0 / (M_PI * h_s6);
    float visc_amp = (float)(1.5 * (double)mass * divgradWvisco * pow((double)sim_scale, 3.0));
    double wpoly6_si = 315.0 / (64.0 * M_PI * h_s9);
    float surf_amp = (float)(-1.7e-9 * (double)mass * wpoly6_si * (double)sim_scale / (double)mass);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLBuffer> bPa = [ctx.device newBufferWithBytes:pos_a length:pos_a_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bVa = [ctx.device newBufferWithBytes:vel_a length:pos_a_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSs = [ctx.device newBufferWithBytes:sort_s length:(size_t)n_static*3*sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bCs = [ctx.device newBufferWithBytes:cs length:(size_t)(n_cells+1)*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bD  = [ctx.device newBufferWithBytes:dens length:(size_t)n_active*sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGe = [ctx.device newBufferWithBytes:gea length:pos_a_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGp = [ctx.device newBufferWithLength:pos_a_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGv = [ctx.device newBufferWithLength:pos_a_b options:MTLResourceStorageModeShared];
        memset([bGp contents], 0, pos_a_b);
        memset([bGv contents], 0, pos_a_b);

        id<MTLBuffer> bH  = [ctx.device newBufferWithBytes:&h length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bH2 = [ctx.device newBufferWithBytes:&h2 length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bM  = [ctx.device newBufferWithBytes:&mass length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSS = [ctx.device newBufferWithBytes:&sim_scale length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bVp = [ctx.device newBufferWithBytes:&visc_pair_coef length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bVa2= [ctx.device newBufferWithBytes:&visc_amp length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bSa = [ctx.device newBufferWithBytes:&surf_amp length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bN  = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];

        int grid_dim_packed[3] = {grid_dim_x, grid_dim_y, grid_dim_z};
        id<MTLBuffer> bGd = [ctx.device newBufferWithBytes:grid_dim_packed length:3*sizeof(int) options:MTLResourceStorageModeShared];
        float gox[3] = {0,0,0};
        FILE *gf = fopen("/tmp/pair_grid_origin.bin","rb");
        if (gf) { fread(gox, sizeof(float), 3, gf); fclose(gf); }
        id<MTLBuffer> bGo = [ctx.device newBufferWithBytes:gox length:3*sizeof(float) options:MTLResourceStorageModeShared];

        id<MTLComputePipelineState> pso = make_pso(ctx, "pair_forces_grid_backward");
        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pso];
        [enc setBuffer:bPa offset:0 atIndex:0];
        [enc setBuffer:bVa offset:0 atIndex:1];
        [enc setBuffer:bSs offset:0 atIndex:2];
        [enc setBuffer:bCs offset:0 atIndex:3];
        [enc setBuffer:bD  offset:0 atIndex:4];
        [enc setBuffer:bGe offset:0 atIndex:5];
        [enc setBuffer:bGp offset:0 atIndex:6];
        [enc setBuffer:bGv offset:0 atIndex:7];
        [enc setBuffer:bH  offset:0 atIndex:8];
        [enc setBuffer:bH2 offset:0 atIndex:9];
        [enc setBuffer:bM  offset:0 atIndex:10];
        [enc setBuffer:bSS offset:0 atIndex:11];
        [enc setBuffer:bVp offset:0 atIndex:12];
        [enc setBuffer:bVa2 offset:0 atIndex:13];
        [enc setBuffer:bSa offset:0 atIndex:14];
        [enc setBuffer:bN  offset:0 atIndex:15];
        [enc setBuffer:bGd offset:0 atIndex:16];
        [enc setBuffer:bGo offset:0 atIndex:17];
        [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
        [enc endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fp = fopen("/tmp/pair_grad_pos.bin","wb");
        fwrite([bGp contents], 1, pos_a_b, fp); fclose(fp);
        FILE *fv = fopen("/tmp/pair_grad_vel.bin","wb");
        fwrite([bGv contents], 1, pos_a_b, fv); fclose(fv);
    }
    free(pos_a); free(vel_a); free(sort_s); free(cs); free(dens); free(gea);
    return 0;
}

// ──────────────────────────────────────────────────────────────────────
// spring_bonds_fwd / _bwd — standalone test ops for spring kernel.
// Used by test_spring_bonds_bwd.py for finite-diff validation.
// ──────────────────────────────────────────────────────────────────────
int run_spring_bonds_fwd(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr, "usage: sib_metal spring_bonds_fwd "
                "<n_active> <n_bonds> <spring_K> "
                "<pos_active.bin> <bond_ij.bin> <bond_rest.bin> "
                "<ext_accel_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_bonds  = (uint32_t)atoi(argv[3]);
    float spring_K    = (float)atof(argv[4]);
    const char *p_pa  = argv[5];
    const char *p_bij = argv[6];
    const char *p_brl = argv[7];

    auto rd = ^(const char *p, size_t bytes) {
        FILE *f = fopen(p, "rb"); if (!f) { fprintf(stderr,"open %s\n",p); exit(1); }
        void *b = malloc(bytes);
        fread(b, 1, bytes, f); fclose(f); return b;
    };
    size_t pos_b = (size_t)n_active * 3 * sizeof(float);
    float *pos = (float *)rd(p_pa, pos_b);
    int *bij_raw = (int *)rd(p_bij, (size_t)n_bonds * 2 * sizeof(int));
    float *bri = (float *)rd(p_brl, (size_t)n_bonds * sizeof(float));

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLBuffer> bP  = [ctx.device newBufferWithBytes:pos length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bB  = [ctx.device newBufferWithBytes:bij_raw length:(size_t)n_bonds*2*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bL  = [ctx.device newBufferWithBytes:bri length:(size_t)n_bonds*sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bE  = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        memset([bE contents], 0, pos_b);  // start from zero so test sees just spring
        id<MTLBuffer> bK  = [ctx.device newBufferWithBytes:&spring_K length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNb = [ctx.device newBufferWithBytes:&n_bonds  length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];

        id<MTLComputePipelineState> pso = make_pso(ctx, "spring_bonds_force");
        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pso];
        [enc setBuffer:bP offset:0 atIndex:0];
        [enc setBuffer:bB offset:0 atIndex:1];
        [enc setBuffer:bL offset:0 atIndex:2];
        [enc setBuffer:bE offset:0 atIndex:3];
        [enc setBuffer:bK offset:0 atIndex:4];
        [enc setBuffer:bNb offset:0 atIndex:5];
        [enc setBuffer:bNa offset:0 atIndex:6];
        [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
        [enc endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fo = fopen("/tmp/spring_ext_accel.bin", "wb");
        fwrite([bE contents], 1, pos_b, fo); fclose(fo);
    }
    free(pos); free(bij_raw); free(bri);
    return 0;
}

int run_spring_bonds_bwd(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr, "usage: sib_metal spring_bonds_bwd "
                "<n_active> <n_bonds> <spring_K> "
                "<pos_active.bin> <bond_ij.bin> <bond_rest.bin> "
                "<grad_ext_accel.bin> <grad_pos_out.bin>\n");
        return 1;
    }
    uint32_t n_active = (uint32_t)atoi(argv[2]);
    uint32_t n_bonds  = (uint32_t)atoi(argv[3]);
    float spring_K    = (float)atof(argv[4]);
    const char *p_pa  = argv[5];
    const char *p_bij = argv[6];
    const char *p_brl = argv[7];
    const char *p_ge  = argv[8];

    auto rd = ^(const char *p, size_t bytes) {
        FILE *f = fopen(p, "rb"); if (!f) { fprintf(stderr,"open %s\n",p); exit(1); }
        void *b = malloc(bytes);
        fread(b, 1, bytes, f); fclose(f); return b;
    };
    size_t pos_b = (size_t)n_active * 3 * sizeof(float);
    float *pos = (float *)rd(p_pa, pos_b);
    int *bij_raw = (int *)rd(p_bij, (size_t)n_bonds * 2 * sizeof(int));
    float *bri = (float *)rd(p_brl, (size_t)n_bonds * sizeof(float));
    float *gea = (float *)rd(p_ge, pos_b);

    @autoreleasepool {
        MetalCtx ctx = make_ctx();
        id<MTLBuffer> bP  = [ctx.device newBufferWithBytes:pos length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bB  = [ctx.device newBufferWithBytes:bij_raw length:(size_t)n_bonds*2*sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bL  = [ctx.device newBufferWithBytes:bri length:(size_t)n_bonds*sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGe = [ctx.device newBufferWithBytes:gea length:pos_b options:MTLResourceStorageModeShared];
        id<MTLBuffer> bGp = [ctx.device newBufferWithLength:pos_b options:MTLResourceStorageModeShared];
        memset([bGp contents], 0, pos_b);
        id<MTLBuffer> bK  = [ctx.device newBufferWithBytes:&spring_K length:sizeof(float) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNb = [ctx.device newBufferWithBytes:&n_bonds  length:sizeof(uint32_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> bNa = [ctx.device newBufferWithBytes:&n_active length:sizeof(uint32_t) options:MTLResourceStorageModeShared];

        id<MTLComputePipelineState> pso = make_pso(ctx, "spring_bonds_force_backward");
        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pso];
        [enc setBuffer:bP  offset:0 atIndex:0];
        [enc setBuffer:bB  offset:0 atIndex:1];
        [enc setBuffer:bL  offset:0 atIndex:2];
        [enc setBuffer:bGe offset:0 atIndex:3];
        [enc setBuffer:bGp offset:0 atIndex:4];
        [enc setBuffer:bK  offset:0 atIndex:5];
        [enc setBuffer:bNb offset:0 atIndex:6];
        [enc setBuffer:bNa offset:0 atIndex:7];
        [enc dispatchThreads:MTLSizeMake(n_active, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
        [enc endEncoding];
        [cmd commit]; [cmd waitUntilCompleted];

        FILE *fp = fopen("/tmp/spring_grad_pos.bin", "wb");
        fwrite([bGp contents], 1, pos_b, fp); fclose(fp);
    }
    free(pos); free(bij_raw); free(bri); free(gea);
    return 0;
}
