#pragma once

// Kernel argument abstraction for the `indexx` kernel.
// Metal signature (sphFluid.metal):
//   kernel void indexx(
//       const device uint2 *particleIndex [[buffer(0)]],
//       constant uint &gridCellCount     [[buffer(1)]],
//       device uint *gridCellIndex       [[buffer(2)]],  (output)
//       constant uint &particleCount     [[buffer(3)]],
//       uint targetCellId [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct IndexxMetalArgs {
  MTL::Buffer *particleIndex; // [[buffer(0)]]
  uint32_t gridCellCount;     // [[buffer(1)]]
  MTL::Buffer *gridCellIndex; // [[buffer(2)]] output
  uint32_t particleCount;     // [[buffer(3)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, particleIndex, 0);
    bindScalar(enc, gridCellCount, 1);
    bindBuffer(enc, gridCellIndex, 2);
    bindScalar(enc, particleCount, 3);
  }
};

// ============ Constants ============
inline constexpr const char *kIndexxKernelName = "indexx";

} // namespace Sibernetic
