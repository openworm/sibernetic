#pragma once

// Kernel argument abstraction for the `sortPostPass` kernel.
// Metal signature (sphFluid.metal):
//   kernel void sortPostPass(
//       const device uint2 *particleIndex     [[buffer(0)]],
//       device uint *sortedParticleIdBySerialId        [[buffer(1)]],  (output)
//       const device float4 *position         [[buffer(2)]],
//       const device float4 *velocity         [[buffer(3)]],
//       device float4 *sortedPosition         [[buffer(4)]],  (output)
//       device float4 *sortedVelocity         [[buffer(5)]],  (output)
//       constant uint &particleCount          [[buffer(6)]],
//       uint particleId [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct SortPostPassMetalArgs {
  MTL::Buffer *particleIndex;              // [[buffer(0)]]
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(1)]] output
  MTL::Buffer *position;                   // [[buffer(2)]]
  MTL::Buffer *velocity;                   // [[buffer(3)]]
  MTL::Buffer *sortedPosition;             // [[buffer(4)]] output
  MTL::Buffer *sortedVelocity;             // [[buffer(5)]] output
  uint32_t particleCount;                  // [[buffer(6)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, particleIndex, 0);
    bindBuffer(enc, sortedParticleIdBySerialId, 1);
    bindBuffer(enc, position, 2);
    bindBuffer(enc, velocity, 3);
    bindBuffer(enc, sortedPosition, 4);
    bindBuffer(enc, sortedVelocity, 5);
    bindScalar(enc, particleCount, 6);
  }
};

// ============ Constants ============
inline constexpr const char *kSortPostPassKernelName = "sortPostPass";

} // namespace Sibernetic
