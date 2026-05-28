#pragma once

// Kernel argument abstraction for the `clearBuffers` kernel.
//
// Clears the neighborMap buffer to "no neighbor" sentinel values (-1, -1)
// for all particles. Each thread handles one particle's kMaxNeighborCount
// (32) float2 entries.
// Metal signature (sphFluid.metal):
//   kernel void clearBuffers(
//       device float4 *neighborMap      [[buffer(0)]],  // cast for efficiency
//       constant uint &particleCount    [[buffer(1)]],
//       uint id [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct ClearBuffersMetalArgs {
  MTL::Buffer *neighborMap; // [[buffer(0)]]  in/out
  uint32_t particleCount;   // [[buffer(1)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, neighborMap, 0);
    bindScalar(enc, particleCount, 1);
  }
};

// ============ Constants ============
inline constexpr const char *kClearBuffersKernelName = "clearBuffers";

} // namespace Sibernetic
