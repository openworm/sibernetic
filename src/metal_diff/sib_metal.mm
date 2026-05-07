// sib_metal.mm — main + op dispatcher.
//
// All op implementations live in the various ops_*.mm files; this file
// just routes argv[1] to the right run_* function.

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <stdio.h>
#include <string.h>

// Forward declarations for the run_* functions exposed by ops_*.mm.
int run_dist_active_static(int argc, char **argv);
int run_dist_active_active(int argc, char **argv);
int run_wpoly6_inplace(int argc, char **argv);
int run_rowsum_density(int argc, char **argv);
int run_density_constraint_grad(int argc, char **argv);
int run_xpbd_step(int argc, char **argv);
int run_step_simple_fwd(int argc, char **argv);
int run_step_simple_bwd(int argc, char **argv);
int run_step_floor_fwd(int argc, char **argv);
int run_step_floor_bwd(int argc, char **argv);
int run_step_bond_fwd(int argc, char **argv);
int run_step_bond_bwd(int argc, char **argv);
int run_density_as_fwd(int argc, char **argv);
int run_density_as_bwd(int argc, char **argv);
int run_density_aa_fwd(int argc, char **argv);
int run_density_aa_bwd(int argc, char **argv);
int run_density_constraint_grad_bwd(int argc, char **argv);
int run_solve_density_constraint_fwd(int argc, char **argv);
int run_solve_density_constraint_bwd(int argc, char **argv);
int run_xpbd_full_fwd(int argc, char **argv);
int run_xpbd_full_bwd(int argc, char **argv);
int run_pair_forces_fwd(int argc, char **argv);
int run_pair_forces_bwd(int argc, char **argv);
int run_spring_bonds_fwd(int argc, char **argv);
int run_spring_bonds_bwd(int argc, char **argv);
int run_accumulate_membrane_fwd(int argc, char **argv);
int run_accumulate_membrane_bwd(int argc, char **argv);
int run_apply_membrane_fwd(int argc, char **argv);
int run_apply_membrane_bwd(int argc, char **argv);


int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: sib_metal <op> [args...]\n"
                        "ops: dist_active_static, dist_active_active, "
                        "wpoly6_inplace, rowsum_density, "
                        "density_constraint_grad, density_constraint_grad_bwd, "
                        "solve_density_constraint_fwd, solve_density_constraint_bwd, "
                        "xpbd_step, xpbd_full_fwd, xpbd_full_bwd, "
                        "step_simple_fwd, step_simple_bwd, "
                        "step_floor_fwd, step_floor_bwd, "
                        "step_bond_fwd, step_bond_bwd, "
                        "density_as_fwd, density_as_bwd, "
                        "density_aa_fwd, density_aa_bwd\n");
        return 1;
    }
    if (strcmp(argv[1], "dist_active_static") == 0)
        return run_dist_active_static(argc, argv);
    if (strcmp(argv[1], "dist_active_active") == 0)
        return run_dist_active_active(argc, argv);
    if (strcmp(argv[1], "wpoly6_inplace") == 0)
        return run_wpoly6_inplace(argc, argv);
    if (strcmp(argv[1], "rowsum_density") == 0)
        return run_rowsum_density(argc, argv);
    if (strcmp(argv[1], "density_constraint_grad") == 0)
        return run_density_constraint_grad(argc, argv);
    if (strcmp(argv[1], "xpbd_step") == 0)
        return run_xpbd_step(argc, argv);
    if (strcmp(argv[1], "step_simple_fwd") == 0)
        return run_step_simple_fwd(argc, argv);
    if (strcmp(argv[1], "step_simple_bwd") == 0)
        return run_step_simple_bwd(argc, argv);
    if (strcmp(argv[1], "step_floor_fwd") == 0)
        return run_step_floor_fwd(argc, argv);
    if (strcmp(argv[1], "step_floor_bwd") == 0)
        return run_step_floor_bwd(argc, argv);
    if (strcmp(argv[1], "step_bond_fwd") == 0)
        return run_step_bond_fwd(argc, argv);
    if (strcmp(argv[1], "step_bond_bwd") == 0)
        return run_step_bond_bwd(argc, argv);
    if (strcmp(argv[1], "density_as_fwd") == 0)
        return run_density_as_fwd(argc, argv);
    if (strcmp(argv[1], "density_as_bwd") == 0)
        return run_density_as_bwd(argc, argv);
    if (strcmp(argv[1], "density_aa_fwd") == 0)
        return run_density_aa_fwd(argc, argv);
    if (strcmp(argv[1], "density_aa_bwd") == 0)
        return run_density_aa_bwd(argc, argv);
    if (strcmp(argv[1], "density_constraint_grad_bwd") == 0)
        return run_density_constraint_grad_bwd(argc, argv);
    if (strcmp(argv[1], "solve_density_constraint_fwd") == 0)
        return run_solve_density_constraint_fwd(argc, argv);
    if (strcmp(argv[1], "solve_density_constraint_bwd") == 0)
        return run_solve_density_constraint_bwd(argc, argv);
    if (strcmp(argv[1], "xpbd_full_fwd") == 0)
        return run_xpbd_full_fwd(argc, argv);
    if (strcmp(argv[1], "xpbd_full_bwd") == 0)
        return run_xpbd_full_bwd(argc, argv);
    if (strcmp(argv[1], "pair_forces_fwd") == 0)
        return run_pair_forces_fwd(argc, argv);
    if (strcmp(argv[1], "pair_forces_bwd") == 0)
        return run_pair_forces_bwd(argc, argv);
    if (strcmp(argv[1], "spring_bonds_fwd") == 0)
        return run_spring_bonds_fwd(argc, argv);
    if (strcmp(argv[1], "spring_bonds_bwd") == 0)
        return run_spring_bonds_bwd(argc, argv);
    if (strcmp(argv[1], "accumulate_membrane_fwd") == 0)
        return run_accumulate_membrane_fwd(argc, argv);
    if (strcmp(argv[1], "accumulate_membrane_bwd") == 0)
        return run_accumulate_membrane_bwd(argc, argv);
    if (strcmp(argv[1], "apply_membrane_fwd") == 0)
        return run_apply_membrane_fwd(argc, argv);
    if (strcmp(argv[1], "apply_membrane_bwd") == 0)
        return run_apply_membrane_bwd(argc, argv);
    fprintf(stderr, "unknown op: %s\n", argv[1]);
    return 1;
}
