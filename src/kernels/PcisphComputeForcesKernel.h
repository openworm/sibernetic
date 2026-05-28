#pragma once

// Kernel argument abstraction for the `pcisph_computeForcesAndInitPressure`
// kernel.
// Metal signature (sphFluid.metal):
//   kernel void pcisph_computeForcesAndInitPressure(
//       const device float2 *neighborMap          [[buffer(0)]],
//       const device float *rho                   [[buffer(1)]],
//       device float *pressure                    [[buffer(2)]],  (output)
//       const device float4 *sortedPosition       [[buffer(3)]],
//       const device float4 *sortedVelocity       [[buffer(4)]],
//       device float4 *acceleration               [[buffer(5)]],  (output)
//       const device uint *sortedParticleIdBySerialId      [[buffer(6)]],
//       constant float &surfTensCoeff             [[buffer(7)]],
//       constant float &massMultLaplacianWviscosityCoeff [[buffer(8)]],
//       constant float &hScaled                   [[buffer(9)]],
//       constant float &mu                        [[buffer(10)]],
//       constant float &gravity_x                 [[buffer(11)]],
//       constant float &gravity_y                 [[buffer(12)]],
//       constant float &gravity_z                 [[buffer(13)]],
//       const device float4 *position             [[buffer(14)]],
//       const device uint2 *particleIndex         [[buffer(15)]],
//       constant uint &particleCount              [[buffer(16)]],
//       constant float &mass                      [[buffer(17)]],
//       uint gid [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct PcisphComputeForcesMetalArgs {
  MTL::Buffer *neighborMap;                // [[buffer(0)]]
  MTL::Buffer *rho;                        // [[buffer(1)]]
  MTL::Buffer *pressure;                   // [[buffer(2)]] output
  MTL::Buffer *sortedPosition;             // [[buffer(3)]]
  MTL::Buffer *sortedVelocity;             // [[buffer(4)]]
  MTL::Buffer *acceleration;               // [[buffer(5)]] output
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(6)]]
  float surfTensCoeff;                     // [[buffer(7)]]
  float massMultLaplacianWviscosityCoeff;  // [[buffer(8)]]
  float hScaled;                           // [[buffer(9)]]
  float mu;                                // [[buffer(10)]]
  float gravity_x;                         // [[buffer(11)]]
  float gravity_y;                         // [[buffer(12)]]
  float gravity_z;                         // [[buffer(13)]]
  MTL::Buffer *position;                   // [[buffer(14)]]
  MTL::Buffer *particleIndex;              // [[buffer(15)]]
  uint32_t particleCount;                  // [[buffer(16)]]
  float mass;                              // [[buffer(17)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, neighborMap, 0);
    bindBuffer(enc, rho, 1);
    bindBuffer(enc, pressure, 2);
    bindBuffer(enc, sortedPosition, 3);
    bindBuffer(enc, sortedVelocity, 4);
    bindBuffer(enc, acceleration, 5);
    bindBuffer(enc, sortedParticleIdBySerialId, 6);
    bindScalar(enc, surfTensCoeff, 7);
    bindScalar(enc, massMultLaplacianWviscosityCoeff, 8);
    bindScalar(enc, hScaled, 9);
    bindScalar(enc, mu, 10);
    bindScalar(enc, gravity_x, 11);
    bindScalar(enc, gravity_y, 12);
    bindScalar(enc, gravity_z, 13);
    bindBuffer(enc, position, 14);
    bindBuffer(enc, particleIndex, 15);
    bindScalar(enc, particleCount, 16);
    bindScalar(enc, mass, 17);
  }
};

// ============ Constants ============
inline constexpr const char *kPcisphComputeForcesKernelName =
    "pcisph_computeForcesAndInitPressure";

} // namespace Sibernetic
