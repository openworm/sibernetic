#pragma once

// Kernel argument abstraction for the `pcisph_computeDensity` kernel.
// Metal signature (sphFluid.metal):
//   kernel void pcisph_computeDensity(
//       const device float2 *neighborMap           [[buffer(0)]],
//       constant float &massMultWpoly6Coefficient  [[buffer(1)]],
//       constant float &hScaled2                   [[buffer(2)]],
//       device float *rho                          [[buffer(3)]],  (output)
//       const device uint *sortedParticleIdBySerialId       [[buffer(4)]],
//       constant uint &particleCount               [[buffer(5)]],
//       uint serialId [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct ComputeDensityMetalArgs {
  MTL::Buffer *neighborMap;                // [[buffer(0)]]
  float massMultWpoly6Coefficient;         // [[buffer(1)]]
  float hScaled2;                          // [[buffer(2)]]
  MTL::Buffer *rho;                        // [[buffer(3)]] output
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(4)]]
  uint32_t particleCount;                  // [[buffer(5)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, neighborMap, 0);
    bindScalar(enc, massMultWpoly6Coefficient, 1);
    bindScalar(enc, hScaled2, 2);
    bindBuffer(enc, rho, 3);
    bindBuffer(enc, sortedParticleIdBySerialId, 4);
    bindScalar(enc, particleCount, 5);
  }
};

// ============ Constants ============
inline constexpr const char *kComputeDensityKernelName =
    "pcisph_computeDensity";

} // namespace Sibernetic
