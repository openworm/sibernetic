// ops.h — forward declarations of all CLI op drivers.
//
// Each `run_<op>` consumes the CLI argv slice for its op and returns an
// exit code. Defined across ops_kernels_m6.cu, ops_xpbd_step.cu, and
// ops_xpbd_full.cu; called from main() in sib_cuda.cu.
#pragma once

// ── M6 forward + backward standalone ops ─────────────────────────────
int run_wpoly6_inplace(int argc, char **argv);
int run_dist_active_static(int argc, char **argv);
int run_dist_active_active(int argc, char **argv);
int run_rowsum_density(int argc, char **argv);
int run_density_constraint_grad(int argc, char **argv);

int run_dist_active_static_bwd(int argc, char **argv);
int run_dist_active_active_bwd(int argc, char **argv);
int run_wpoly6_inplace_bwd(int argc, char **argv);
int run_rowsum_density_bwd(int argc, char **argv);
int run_density_constraint_grad_bwd(int argc, char **argv);
int run_solve_density_constraint_bwd(int argc, char **argv);

// ── xpbd_step (all-in-one orchestrator) ──────────────────────────────
int run_xpbd_step(int argc, char **argv);

// ── xpbd_full_fwd / xpbd_full_bwd (ops_xpbd_full.cu) ─────────────────
// Per spec Phase 4. Not yet implemented; pair_forces_grid + Metal-parity
// CLI signature are prerequisites. Add prototypes here when those land.

// ── Pair forces + spring bonds (ops_pair_spring.cu) ──────────────────
int run_spring_bonds_force(int argc, char **argv);
int run_apply_ext_accel(int argc, char **argv);
int run_spring_bonds_force_bwd(int argc, char **argv);
int run_spring_K_partial(int argc, char **argv);
int run_apply_ext_accel_bwd(int argc, char **argv);
int run_pair_forces_grid_fwd(int argc, char **argv);
