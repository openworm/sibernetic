#pragma once

// Kernel argument abstraction for the `pcisph_predictPositions` kernel.
// Metal signature (sphFluid.metal):
//   kernel void pcisph_predictPositions(
//       device float4 *acceleration             [[buffer(0)]],
//       device float4 *sortedPosition           [[buffer(1)]],
//       const device float4 *sortedVelocity     [[buffer(2)]],
//       const device uint2 *sortedCellAndSerialId [[buffer(3)]],
//       const device uint *sortedParticleIdBySerialId [[buffer(4)]],
//       constant float &gravitationalAccelerationX               [[buffer(5)]],
//       constant float &gravitationalAccelerationY               [[buffer(6)]],
//       constant float &gravitationalAccelerationZ               [[buffer(7)]],
//       constant float &simulationScaleInv      [[buffer(8)]],
//       constant float &deltaTime               [[buffer(9)]],
//       const device float4 *originalPosition           [[buffer(10)]],
//       const device float4 *velocity           [[buffer(11)]],
//       constant float &r0                      [[buffer(12)]],
//       const device float2 *neighborMap        [[buffer(13)]],
//       constant uint &particleCount            [[buffer(14)]],
//       uint serialId [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct PcisphPredictPositionsMetalArgs {
  MTL::Buffer *acceleration;               // [[buffer(0)]]  in/out (3×N)
  MTL::Buffer *sortedPosition;             // [[buffer(1)]]  in/out (2×N)
  MTL::Buffer *sortedVelocity;             // [[buffer(2)]]
  MTL::Buffer *sortedCellAndSerialId;      // [[buffer(3)]]
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(4)]]
  float gravitationalAccelerationX;        // [[buffer(5)]]
  float gravitationalAccelerationY;        // [[buffer(6)]]
  float gravitationalAccelerationZ;        // [[buffer(7)]]
  float simulationScaleInv;                // [[buffer(8)]]
  float deltaTime;                         // [[buffer(9)]]
  MTL::Buffer *originalPosition;           // [[buffer(10)]]
  MTL::Buffer *velocity;                   // [[buffer(11)]]
  float r0;                                // [[buffer(12)]]
  MTL::Buffer *neighborMap;                // [[buffer(13)]]
  uint32_t particleCount;                  // [[buffer(14)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, acceleration, 0);
    bindBuffer(enc, sortedPosition, 1);
    bindBuffer(enc, sortedVelocity, 2);
    bindBuffer(enc, sortedCellAndSerialId, 3);
    bindBuffer(enc, sortedParticleIdBySerialId, 4);
    bindScalar(enc, gravitationalAccelerationX, 5);
    bindScalar(enc, gravitationalAccelerationY, 6);
    bindScalar(enc, gravitationalAccelerationZ, 7);
    bindScalar(enc, simulationScaleInv, 8);
    bindScalar(enc, deltaTime, 9);
    bindBuffer(enc, originalPosition, 10);
    bindBuffer(enc, velocity, 11);
    bindScalar(enc, r0, 12);
    bindBuffer(enc, neighborMap, 13);
    bindScalar(enc, particleCount, 14);
  }
};

// ============ Constants ============
inline constexpr const char *kPcisphPredictPositionsKernelName =
    "pcisph_predictPositions";

} // namespace Sibernetic
