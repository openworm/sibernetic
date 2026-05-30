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
// Per spec Phase 4. Differentiable XPBD forward + backward.
int run_xpbd_full_fwd(int argc, char **argv);
int run_xpbd_full_bwd(int argc, char **argv);

// ── Pair forces + spring bonds (ops_pair_spring.cu) ──────────────────
int run_spring_bonds_force(int argc, char **argv);
int run_apply_ext_accel(int argc, char **argv);
int run_spring_bonds_force_bwd(int argc, char **argv);
int run_spring_K_partial(int argc, char **argv);
int run_apply_ext_accel_bwd(int argc, char **argv);
int run_pair_forces_grid_fwd(int argc, char **argv);
int run_pair_forces_grid_bwd(int argc, char **argv);
int run_visc_K_partial(int argc, char **argv);
