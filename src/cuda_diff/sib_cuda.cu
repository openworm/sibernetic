/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013, 2026 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/

// sib_cuda.cu — entry point + CLI dispatcher.
// All op handlers live in ops_*.cu; kernels in shaders.cu.

#include "cuda_common.h"
#include "shaders.h"
#include "ops.h"

// ──────────────────────────────────────────────────────────────────────
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr,
            "usage: sib_cuda <op> [args...]\n"
            "\n"
            "Integrated XPBD ops (shared CLI with sib_metal):\n"
            "  xpbd_step          imperative multi-step driver\n"
            "  xpbd_full_fwd      differentiable forward (writes tape)\n"
            "  xpbd_full_bwd      differentiable backward (reads tape)\n"
            "\n"
            "M6 / density-chain primitives (fwd):\n"
            "  dist_active_static, dist_active_active,\n"
            "  wpoly6_inplace, rowsum_density,\n"
            "  density_constraint_grad\n"
            "\n"
            "M6 / density-chain primitives (bwd, FD-validated):\n"
            "  dist_active_static_bwd, dist_active_active_bwd,\n"
            "  wpoly6_inplace_bwd, rowsum_density_bwd,\n"
            "  density_constraint_grad_bwd,\n"
            "  solve_density_constraint_bwd\n"
            "\n"
            "Worm physics (pair forces, springs, external acceleration):\n"
            "  pair_forces_grid_fwd, pair_forces_grid_bwd,\n"
            "  spring_bonds_force, spring_bonds_force_bwd,\n"
            "  apply_ext_accel, apply_ext_accel_bwd,\n"
            "  visc_K_partial, spring_K_partial\n"
            "\n"
            "Run any op with no further args to see its argv layout.\n");
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
    if (std::strcmp(argv[1], "xpbd_full_bwd") == 0)
        return run_xpbd_full_bwd(argc, argv);
    fprintf(stderr, "unknown op: %s\n", argv[1]);
    return 1;
}
