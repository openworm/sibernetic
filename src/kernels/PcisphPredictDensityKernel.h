#pragma once

// Kernel argument abstraction for the `pcisph_predictDensity` kernel.
//
// Computes predicted density from predicted positions (sortedPosition[N..2N))
// and writes the result to rho[N..2N). Used inside the PCISPH pressure
// correction loop after pcisph_predictPositions updates the predicted
// positions.
// Metal signature (sphFluid.metal):
//   kernel void pcisph_predictDensity(
//       const device float2 *neighborMap                  [[buffer(0)]],
//       const device uint   *sortedParticleIdBySerialId   [[buffer(1)]],
//       constant float &massMultWpoly6Coefficient         [[buffer(2)]],
//       constant float &h                                 [[buffer(3)]],
//       constant float &restDensity                       [[buffer(4)]],
//       constant float &simulationScale                   [[buffer(5)]],
//       const device float4 *sortedPosition               [[buffer(6)]],
//       device float *rho                                 [[buffer(7)]],
//       constant uint &particleCount                      [[buffer(8)]],
//       uint serialId [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct PcisphPredictDensityMetalArgs {
  MTL::Buffer *neighborMap;                // [[buffer(0)]]
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(1)]]
  float massMultWpoly6Coefficient;         // [[buffer(2)]]
  float h;                                 // [[buffer(3)]]
  float restDensity;                       // [[buffer(4)]]
  float simulationScale;                   // [[buffer(5)]]
  MTL::Buffer *sortedPosition;             // [[buffer(6)]]
  MTL::Buffer *rho;                        // [[buffer(7)]] output
  uint32_t particleCount;                  // [[buffer(8)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, neighborMap, 0);
    bindBuffer(enc, sortedParticleIdBySerialId, 1);
    bindScalar(enc, massMultWpoly6Coefficient, 2);
    bindScalar(enc, h, 3);
    bindScalar(enc, restDensity, 4);
    bindScalar(enc, simulationScale, 5);
    bindBuffer(enc, sortedPosition, 6);
    bindBuffer(enc, rho, 7);
    bindScalar(enc, particleCount, 8);
  }
};

// ============ Constants ============
inline constexpr const char *kPcisphPredictDensityKernelName =
    "pcisph_predictDensity";

} // namespace Sibernetic
