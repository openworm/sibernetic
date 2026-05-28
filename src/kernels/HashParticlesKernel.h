#pragma once

// Kernel argument abstraction for the `hashParticles` kernel.
// Metal signature (sphFluid.metal):
//   kernel void hashParticles(
//       const device float4 *position        [[buffer(0)]],
//       constant uint &gridCellsX            [[buffer(1)]],
//       constant uint &gridCellsY            [[buffer(2)]],
//       constant uint &gridCellsZ            [[buffer(3)]],
//       constant float &hashGridCellSizeInv  [[buffer(4)]],
//       constant float &xmin                 [[buffer(5)]],
//       constant float &ymin                 [[buffer(6)]],
//       constant float &zmin                 [[buffer(7)]],
//       device uint2 *particleIndex          [[buffer(8)]],  (output)
//       constant uint &particleCount         [[buffer(9)]],
//       uint particleId [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct HashParticlesMetalArgs {
  MTL::Buffer *position;      // [[buffer(0)]]
  uint32_t gridCellsX;        // [[buffer(1)]]
  uint32_t gridCellsY;        // [[buffer(2)]]
  uint32_t gridCellsZ;        // [[buffer(3)]]
  float hashGridCellSizeInv;  // [[buffer(4)]]
  float xmin;                 // [[buffer(5)]]
  float ymin;                 // [[buffer(6)]]
  float zmin;                 // [[buffer(7)]]
  MTL::Buffer *particleIndex; // [[buffer(8)]] output
  uint32_t particleCount;     // [[buffer(9)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, position, 0);
    bindScalar(enc, gridCellsX, 1);
    bindScalar(enc, gridCellsY, 2);
    bindScalar(enc, gridCellsZ, 3);
    bindScalar(enc, hashGridCellSizeInv, 4);
    bindScalar(enc, xmin, 5);
    bindScalar(enc, ymin, 6);
    bindScalar(enc, zmin, 7);
    bindBuffer(enc, particleIndex, 8);
    bindScalar(enc, particleCount, 9);
  }
};

// ============ Constants ============
inline constexpr const char *kHashParticlesKernelName = "hashParticles";

} // namespace Sibernetic
