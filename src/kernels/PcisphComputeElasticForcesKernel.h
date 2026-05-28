#pragma once

// Kernel argument abstraction for the `pcisph_computeElasticForces` kernel.
//
// Computes elastic spring forces (Hooke's law) and muscle activation forces
// for the worm body simulation. Thread index is the elastic particle index
// (0..numOfElasticP-1).
// Metal signature (sphFluid.metal):
//   kernel void pcisph_computeElasticForces(
//       const device float4 *sortedPosition               [[buffer(0)]],
//       device float4 *acceleration                       [[buffer(1)]],
//       const device uint *sortedParticleIdBySerialId     [[buffer(2)]],
//       const device uint2 *sortedCellAndSerialId         [[buffer(3)]],
//       constant float &maxMuscleForce                    [[buffer(4)]],
//       constant float &simulationScale                   [[buffer(5)]],
//       constant uint &numOfElasticP                      [[buffer(6)]],
//       const device float4 *elasticConnectionsData       [[buffer(7)]],
//       constant uint &muscleCount                        [[buffer(8)]],
//       const device float *muscleActivationSignal        [[buffer(9)]],
//       const device float4 *originalPosition             [[buffer(10)]],
//       constant float &elasticityCoefficient             [[buffer(11)]],
//       uint index [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct PcisphComputeElasticForcesMetalArgs {
  MTL::Buffer *sortedPosition;             // [[buffer(0)]]
  MTL::Buffer *acceleration;               // [[buffer(1)]]  in/out
  MTL::Buffer *sortedParticleIdBySerialId; // [[buffer(2)]]
  MTL::Buffer *sortedCellAndSerialId;      // [[buffer(3)]]
  float maxMuscleForce;                    // [[buffer(4)]]
  float simulationScale;                   // [[buffer(5)]]
  uint32_t numOfElasticP;                  // [[buffer(6)]]
  MTL::Buffer *elasticConnectionsData;     // [[buffer(7)]]
  uint32_t muscleCount;                    // [[buffer(8)]]
  MTL::Buffer *muscleActivationSignal;     // [[buffer(9)]]
  MTL::Buffer *originalPosition;           // [[buffer(10)]]
  float elasticityCoefficient;             // [[buffer(11)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, sortedPosition, 0);
    bindBuffer(enc, acceleration, 1);
    bindBuffer(enc, sortedParticleIdBySerialId, 2);
    bindBuffer(enc, sortedCellAndSerialId, 3);
    bindScalar(enc, maxMuscleForce, 4);
    bindScalar(enc, simulationScale, 5);
    bindScalar(enc, numOfElasticP, 6);
    bindBuffer(enc, elasticConnectionsData, 7);
    bindScalar(enc, muscleCount, 8);
    bindBuffer(enc, muscleActivationSignal, 9);
    bindBuffer(enc, originalPosition, 10);
    bindScalar(enc, elasticityCoefficient, 11);
  }
};

// ============ Constants ============
inline constexpr const char *kPcisphComputeElasticForcesKernelName =
    "pcisph_computeElasticForces";

} // namespace Sibernetic
