#pragma once

// Kernel argument abstractions for the `computeInteractionWithMembranes` and
// `computeInteractionWithMembranes_finalize` kernels.
//
// computeInteractionWithMembranes: For each liquid particle, finds neighboring
// elastic (membrane) particles, projects the particle onto membrane planes,
// and accumulates a position correction into position[N + serialId].
//
// computeInteractionWithMembranes_finalize: Applies the accumulated delta
// from position[N + serialId] to position[serialId].
//   kernel void computeInteractionWithMembranes_finalize(
//       device float4 *position                       [[buffer(0)]],
//       const device uint *sortedParticleIdBySerialId [[buffer(1)]],
//       constant uint &particleCount                  [[buffer(2)]],
//       uint serialId [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

// Maximum membranes sharing a single elastic particle vertex.
inline constexpr int kMaxMembranesIncludingSameParticle = 7;

// ============ computeInteractionWithMembranes ============

struct ComputeInteractionWithMembranesMetalArgs {
  MTL::Buffer *position;                   // [[buffer(0)]]
  MTL::Buffer *velocity;                   // [[buffer(1)]]
  MTL::Buffer *sortedCellAndSerialId;      // [[buffer(2)]]
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(3)]]
  MTL::Buffer *neighborMap;                // [[buffer(4)]]
  MTL::Buffer *particleMembranesList;      // [[buffer(5)]]
  MTL::Buffer *membraneData;               // [[buffer(6)]]
  uint32_t particleCount;                  // [[buffer(7)]]
  float r0;                                // [[buffer(8)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, position, 0);
    bindBuffer(enc, velocity, 1);
    bindBuffer(enc, sortedCellAndSerialId, 2);
    bindBuffer(enc, sortedParticleIdBySerialId, 3);
    bindBuffer(enc, neighborMap, 4);
    bindBuffer(enc, particleMembranesList, 5);
    bindBuffer(enc, membraneData, 6);
    bindScalar(enc, particleCount, 7);
    bindScalar(enc, r0, 8);
  }
};

inline constexpr const char *kComputeInteractionWithMembranesKernelName =
    "computeInteractionWithMembranes";

// ============ computeInteractionWithMembranes_finalize ============

struct ComputeInteractionWithMembranesFinalizeMetalArgs {
  MTL::Buffer *position;                   // [[buffer(0)]]
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(1)]]
  uint32_t particleCount;                  // [[buffer(2)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, position, 0);
    bindBuffer(enc, sortedParticleIdBySerialId, 1);
    bindScalar(enc, particleCount, 2);
  }
};

inline constexpr const char
    *kComputeInteractionWithMembranesFinalizeKernelName =
        "computeInteractionWithMembranes_finalize";

} // namespace Sibernetic
