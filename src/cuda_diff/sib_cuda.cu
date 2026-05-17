// AUTO-GENERATED from sib_cuda.cu via _split_sib_cuda.py
// Part of the #15 refactor: separates kernel definitions, op drivers,
// and the dispatcher into per-concern translation units linked via
// nvcc -rdc=true.

#include "cuda_common.h"
#include "shaders.h"
#include "ops.h"

// ──────────────────────────────────────────────────────────────────────
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: sib_cuda <op> [args...]\n"
                        "ops: dist_active_static, dist_active_active, "
                        "wpoly6_inplace, rowsum_density, "
                        "density_constraint_grad, xpbd_step\n");
        return 1;
    }
    if (std::strcmp(argv[1], "dist_active_static") == 0)
        return run_dist_active_static(argc, argv);
    if (std::strcmp(argv[1], "dist_active_active") == 0)
        return run_dist_active_active(argc, argv);
    if (std::strcmp(argv[1], "wpoly6_inplace") == 0)
        return run_wpoly6_inplace(argc, argv);
    if (std::strcmp(argv[1], "rowsum_density") == 0)
        return run_rowsum_density(argc, argv);
    if (std::strcmp(argv[1], "density_constraint_grad") == 0)
        return run_density_constraint_grad(argc, argv);
    if (std::strcmp(argv[1], "xpbd_step") == 0)
        return run_xpbd_step(argc, argv);
    if (std::strcmp(argv[1], "dist_active_static_bwd") == 0)
        return run_dist_active_static_bwd(argc, argv);
    if (std::strcmp(argv[1], "dist_active_active_bwd") == 0)
        return run_dist_active_active_bwd(argc, argv);
    if (std::strcmp(argv[1], "wpoly6_inplace_bwd") == 0)
        return run_wpoly6_inplace_bwd(argc, argv);
    if (std::strcmp(argv[1], "rowsum_density_bwd") == 0)
        return run_rowsum_density_bwd(argc, argv);
    if (std::strcmp(argv[1], "density_constraint_grad_bwd") == 0)
        return run_density_constraint_grad_bwd(argc, argv);
    if (std::strcmp(argv[1], "solve_density_constraint_bwd") == 0)
        return run_solve_density_constraint_bwd(argc, argv);
    if (std::strcmp(argv[1], "spring_bonds_force") == 0)
        return run_spring_bonds_force(argc, argv);
    if (std::strcmp(argv[1], "apply_ext_accel") == 0)
        return run_apply_ext_accel(argc, argv);
    if (std::strcmp(argv[1], "spring_bonds_force_bwd") == 0)
        return run_spring_bonds_force_bwd(argc, argv);
    if (std::strcmp(argv[1], "spring_K_partial") == 0)
        return run_spring_K_partial(argc, argv);
    if (std::strcmp(argv[1], "apply_ext_accel_bwd") == 0)
        return run_apply_ext_accel_bwd(argc, argv);
    if (std::strcmp(argv[1], "pair_forces_grid_fwd") == 0)
        return run_pair_forces_grid_fwd(argc, argv);
    if (std::strcmp(argv[1], "pair_forces_grid_bwd") == 0)
        return run_pair_forces_grid_bwd(argc, argv);
    if (std::strcmp(argv[1], "visc_K_partial") == 0)
        return run_visc_K_partial(argc, argv);
    if (std::strcmp(argv[1], "xpbd_full_fwd") == 0)
        return run_xpbd_full_fwd(argc, argv);
    fprintf(stderr, "unknown op: %s\n", argv[1]);
    return 1;
}
