#pragma once

// Kernel argument abstraction for the `pcisph_integrate` kernel.
//
// Performs leapfrog time integration for PCISPH. Has three modes:
//
//   timestepIndex == 0:
//     Store combined acceleration: acceleration[2N + serialId] =
//       acceleration[sortedId] + acceleration[N + sortedId].
//     No position or velocity update.
//
//   mode == 0 (positions):
//     Leapfrog position step:
//       x(t+dt) = x(t) + v(t)*dt + a(t)*dt^2/2
//     Writes into sortedPosition[sortedId] in-place.
//
//   mode == 1 (velocities):
//     Leapfrog velocity step with boundary interaction:
//       v(t+dt) = v(t) + (a(t) + a(t+dt))*dt/2
//     Applies boundary correction, then writes velocity, position, and
//     stores combined acceleration for next step.
//
//   mode == 2 (semi-implicit Euler):
//       v(t+dt) = v(t) + a(t+dt)*dt
//       x(t+dt) = x(t) + v(t+dt)*dt
//     Applies boundary correction, then writes everything.
// Metal signature (sphFluid.metal):
//   kernel void pcisph_integrate(
//       device float4 *acceleration                      [[buffer(0)]],
//       device float4 *sortedPosition                    [[buffer(1)]],
//       const device float4 *sortedVelocity              [[buffer(2)]],
//       const device uint2 *sortedCellAndSerialId        [[buffer(3)]],
//       const device uint *sortedParticleIdBySerialId    [[buffer(4)]],
//       constant float &simulationScaleInv               [[buffer(5)]],
//       constant float &deltaTime                        [[buffer(6)]],
//       device float4 *originalPosition                  [[buffer(7)]],
//       device float4 *velocity                          [[buffer(8)]],
//       constant float &r0                               [[buffer(9)]],
//       const device float2 *neighborMap                 [[buffer(10)]],
//       constant uint &particleCount                     [[buffer(11)]],
//       constant int &timestepIndex                      [[buffer(12)]],
//       constant int &mode                               [[buffer(13)]],
//       uint serialId [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct PcisphIntegrateMetalArgs {
  MTL::Buffer *acceleration;               // [[buffer(0)]]  in/out (3×N)
  MTL::Buffer *sortedPosition;             // [[buffer(1)]]  in/out
  MTL::Buffer *sortedVelocity;             // [[buffer(2)]]
  MTL::Buffer *sortedCellAndSerialId;      // [[buffer(3)]]
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(4)]]
  float simulationScaleInv;                // [[buffer(5)]]
  float deltaTime;                         // [[buffer(6)]]
  MTL::Buffer *originalPosition;           // [[buffer(7)]]  in/out
  MTL::Buffer *velocity;                   // [[buffer(8)]]  in/out
  float r0;                                // [[buffer(9)]]
  MTL::Buffer *neighborMap;                // [[buffer(10)]]
  uint32_t particleCount;                  // [[buffer(11)]]
  int32_t timestepIndex;                   // [[buffer(12)]]
  int32_t mode;                            // [[buffer(13)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, acceleration, 0);
    bindBuffer(enc, sortedPosition, 1);
    bindBuffer(enc, sortedVelocity, 2);
    bindBuffer(enc, sortedCellAndSerialId, 3);
    bindBuffer(enc, sortedParticleIdBySerialId, 4);
    bindScalar(enc, simulationScaleInv, 5);
    bindScalar(enc, deltaTime, 6);
    bindBuffer(enc, originalPosition, 7);
    bindBuffer(enc, velocity, 8);
    bindScalar(enc, r0, 9);
    bindBuffer(enc, neighborMap, 10);
    bindScalar(enc, particleCount, 11);
    bindScalar(enc, timestepIndex, 12);
    bindScalar(enc, mode, 13);
  }
};

// ============ Constants ============
inline constexpr const char *kPcisphIntegrateKernelName = "pcisph_integrate";

} // namespace Sibernetic
