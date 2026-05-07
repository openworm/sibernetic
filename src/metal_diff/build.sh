#!/bin/bash
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

# Modules:
#   metal_common.{h,mm}   shared MetalCtx/make_ctx/make_pso + spatial grid
#   ops_kernels_m6.mm     atomic M6 kernel ops (FD-validated building blocks)
#   ops_xpbd_step.mm      M7 imperative XPBD pipeline
#   ops_xpbd_full.mm      differentiable xpbd_full_fwd + xpbd_full_bwd
#   ops_pair_spring.mm    pair_forces + spring_bonds (fwd/bwd, FD-validated)
#   ops_test_steps.mm     step_simple/floor/bond test ops
#   ops_test_density.mm   density-pipeline test ops
#   sib_metal.mm          slim main + op dispatcher

clang++ -std=c++17 -O2 \
    -framework Metal -framework Foundation \
    metal_common.mm \
    ops_kernels_m6.mm \
    ops_xpbd_step.mm \
    ops_xpbd_full.mm \
    ops_pair_spring.mm \
    ops_membrane.mm \
    ops_test_steps.mm \
    ops_test_density.mm \
    sib_metal.mm \
    -o sib_metal

echo "Built $DIR/sib_metal"
