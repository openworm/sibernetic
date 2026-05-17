#!/bin/bash
# Linux CUDA build for the differentiable Sibernetic substrate.
# Mirrors src/metal_diff/build.sh, swapping clang++ + Metal for nvcc.
#
# Requirements:
#   - CUDA Toolkit 12.0+
#   - A host C++17 compiler nvcc can drive (gcc/clang)
#
# Architecture: sm_75 (Turing) floor. Bump via:
#   CUDA_ARCH=sm_86 ./build.sh    # Ampere / RTX 30-series
#   CUDA_ARCH=sm_89 ./build.sh    # Ada / RTX 40-series / L4
#
# Layout (post #15 refactor): six translation units linked via
# nvcc -rdc=true. Device-link is performed by nvcc automatically in
# the final invocation.
#   cuda_common.cu         host I/O helpers
#   shaders.cu             all __global__ kernel definitions
#   ops_kernels_m6.cu      M6 fwd/bwd + M7.1 bwd standalone drivers
#   ops_xpbd_step.cu       run_xpbd_step (all-in-one orchestrator)
#   ops_xpbd_full.cu       xpbd_full / xpbd_step_diff / xpbd_step_full_diff
#   sib_cuda.cu            main() + CLI dispatcher

set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

CUDA_ARCH="${CUDA_ARCH:-sm_75}"

nvcc -std=c++17 -O2 -arch="$CUDA_ARCH" -rdc=true \
    cuda_common.cu \
    shaders.cu \
    ops_kernels_m6.cu \
    ops_xpbd_step.cu \
    ops_xpbd_full.cu \
    ops_pair_spring.cu \
    sib_cuda.cu \
    -o sib_cuda

echo "Built $DIR/sib_cuda (arch=$CUDA_ARCH)"
