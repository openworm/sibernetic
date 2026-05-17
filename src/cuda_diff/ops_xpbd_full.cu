// ops_xpbd_full.cu — differentiable forward + backward of the full
// XPBD pipeline (predict + density + bonds + pair-forces + floor + update).
//
// Per src/cuda/README.md Phase 4: this file exposes `xpbd_full_fwd` and
// `xpbd_full_bwd` CLI ops whose signatures match metal_diff/sib_metal so
// the shared sgd_true.py runs against either substrate by only changing
// the BIN path. Force-based spring bonds via spring_K (NOT XPBD distance
// constraints), pair forces via visc_pair_coef, and the optional floor
// + restitution flags follow the Metal driver in ops_xpbd_full.mm.
//
// TODO: implementation pending. The previous monolithic sib_cuda.cu
// shipped non-spec-compliant variants (xpbd_full_fwd_step, xpbd_step_diff_*,
// xpbd_step_full_diff_*) that used XPBD distance constraints with
// alpha_dist + had different argv layouts than Metal. Those were removed
// during the spec-alignment pass. The proper xpbd_full_fwd/_bwd will
// land alongside pair_forces_grid + grid construction kernels.

#include "cuda_common.h"
#include "shaders.h"
#include "ops.h"
