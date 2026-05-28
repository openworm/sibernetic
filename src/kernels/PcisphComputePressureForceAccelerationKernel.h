#pragma once

// Kernel argument abstraction for the
// `pcisph_computePressureForceAcceleration` kernel.
//
// Computes the pressure-gradient force acceleration for each fluid particle
// and writes it into acceleration[N..2N). Boundary particles get zero.
// This kernel runs inside the PCISPH predict-correct loop, after
// pcisph_correctPressure has updated the pressure field.
//
// The pressure force follows Solenthaler & Pajarola (2009) formula (5):
//   F_pressure_i = -m * sum_j (p_i + p_j) / (2 * rho_j) * grad_W_spiky(r_ij)
// with a close-range correction when r_ij < hScaled/4 for stability.
// Metal signature (sphFluid.metal):
//   kernel void pcisph_computePressureForceAcceleration(
//       const device float2 *neighborMap                       [[buffer(0)]],
//       const device float  *pressure                          [[buffer(1)]],
//       const device float  *rho                               [[buffer(2)]],
//       const device float4 *sortedPosition                    [[buffer(3)]],
//       const device uint   *sortedParticleIdBySerialId        [[buffer(4)]],
//       constant float &delta                                  [[buffer(5)]],
//       constant float &massMultGradWspikyCoefficient          [[buffer(6)]],
//       constant float &h                                      [[buffer(7)]],
//       constant float &simulationScale                        [[buffer(8)]],
//       constant float &restDensity                            [[buffer(9)]],
//       device float4 *acceleration                            [[buffer(10)]],
//       const device float4 *originalPosition                  [[buffer(11)]],
//       const device uint2  *sortedCellAndSerialId             [[buffer(12)]],
//       constant uint  &particleCount                          [[buffer(13)]],
//       uint serialId [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct PcisphComputePressureForceAccelerationMetalArgs {
  MTL::Buffer *neighborMap;                // [[buffer(0)]]
  MTL::Buffer *pressure;                   // [[buffer(1)]]
  MTL::Buffer *rho;                        // [[buffer(2)]]
  MTL::Buffer *sortedPosition;             // [[buffer(3)]]
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(4)]]
  float delta;                             // [[buffer(5)]]
  float massMultGradWspikyCoefficient;     // [[buffer(6)]]
  float h;                                 // [[buffer(7)]]
  float simulationScale;                   // [[buffer(8)]]
  float restDensity;                       // [[buffer(9)]]
  MTL::Buffer *acceleration;               // [[buffer(10)]] output
  MTL::Buffer *originalPosition;           // [[buffer(11)]]
  MTL::Buffer *sortedCellAndSerialId;      // [[buffer(12)]]
  uint32_t particleCount;                  // [[buffer(13)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, neighborMap, 0);
    bindBuffer(enc, pressure, 1);
    bindBuffer(enc, rho, 2);
    bindBuffer(enc, sortedPosition, 3);
    bindBuffer(enc, sortedParticleIdBySerialId, 4);
    bindScalar(enc, delta, 5);
    bindScalar(enc, massMultGradWspikyCoefficient, 6);
    bindScalar(enc, h, 7);
    bindScalar(enc, simulationScale, 8);
    bindScalar(enc, restDensity, 9);
    bindBuffer(enc, acceleration, 10);
    bindBuffer(enc, originalPosition, 11);
    bindBuffer(enc, sortedCellAndSerialId, 12);
    bindScalar(enc, particleCount, 13);
  }
};

// ============ Constants ============
inline constexpr const char *kPcisphComputePressureForceAccelerationKernelName =
    "pcisph_computePressureForceAcceleration";

} // namespace Sibernetic
