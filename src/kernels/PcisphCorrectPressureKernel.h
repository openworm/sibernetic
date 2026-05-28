#pragma once

// Kernel argument abstraction for the `pcisph_correctPressure` kernel.
//
// Reads predicted density from rho[N..2N) (written by pcisph_predictDensity),
// computes the density error against the reference density rho0, scales by the
// precomputed correction factor delta, clamps to non-negative, and accumulates
// the correction into pressure[0..N).
// Metal signature (sphFluid.metal):
//   kernel void pcisph_correctPressure(
//       const device uint  *sortedParticleIdBySerialId  [[buffer(0)]],
//       constant float &restDensity               [[buffer(1)]],
//       device float *pressure                 [[buffer(2)]],
//       const device float *rho                [[buffer(3)]],
//       constant float &delta                  [[buffer(4)]],
//       constant uint  &particleCount          [[buffer(5)]],
//       uint gid [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct PcisphCorrectPressureMetalArgs {
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(0)]]
  float restDensity;                       // [[buffer(1)]]
  MTL::Buffer *pressure;                   // [[buffer(2)]] input/output
  MTL::Buffer *rho;                        // [[buffer(3)]]
  float delta;                             // [[buffer(4)]]
  uint32_t particleCount;                  // [[buffer(5)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, sortedParticleIdBySerialId, 0);
    bindScalar(enc, restDensity, 1);
    bindBuffer(enc, pressure, 2);
    bindBuffer(enc, rho, 3);
    bindScalar(enc, delta, 4);
    bindScalar(enc, particleCount, 5);
  }
};

// ============ Constants ============
inline constexpr const char *kPcisphCorrectPressureKernelName =
    "pcisph_correctPressure";

} // namespace Sibernetic
