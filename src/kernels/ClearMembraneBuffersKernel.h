#pragma once

// Kernel argument abstraction for the `clearMembraneBuffers` kernel.
//
// Zeros the delta accumulator region [N..2N) of the position and velocity
// buffers. These regions store per-particle membrane interaction deltas
// that must be cleared at the start of each timestep.
// Metal signature (sphFluid.metal):
//   kernel void clearMembraneBuffers(
//       device float4 *position          [[buffer(0)]],  (2×N)
//       device float4 *velocity          [[buffer(1)]],  (2×N)
//       constant uint &particleCount     [[buffer(2)]],
//       uint id [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct ClearMembraneBuffersMetalArgs {
  MTL::Buffer *position;  // [[buffer(0)]]  in/out (2×N)
  MTL::Buffer *velocity;  // [[buffer(1)]]  in/out (2×N)
  uint32_t particleCount; // [[buffer(2)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, position, 0);
    bindBuffer(enc, velocity, 1);
    bindScalar(enc, particleCount, 2);
  }
};

// ============ Constants ============
inline constexpr const char *kClearMembraneBuffersKernelName =
    "clearMembraneBuffers";

} // namespace Sibernetic
