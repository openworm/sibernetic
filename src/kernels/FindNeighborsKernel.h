#pragma once

// Kernel argument abstraction for the `findNeighbors` kernel.
// Metal signature (sphFluid.metal):
//   kernel void findNeighbors(
//       const device uint *gridCellIndicesFixedUp [[buffer(0)]],
//       const device float4 *sortedPosition      [[buffer(1)]],
//       constant uint &gridCellCount             [[buffer(2)]],
//       constant uint &gridCellsX                [[buffer(3)]],
//       constant uint &gridCellsY                [[buffer(4)]],
//       constant uint &gridCellsZ                [[buffer(5)]],
//       constant float &h                        [[buffer(6)]],
//       constant float &hashGridCellSize         [[buffer(7)]],
//       constant float &hashGridCellSizeInv      [[buffer(8)]],
//       constant float &simulationScale          [[buffer(9)]],
//       constant float &xmin                     [[buffer(10)]],
//       constant float &ymin                     [[buffer(11)]],
//       constant float &zmin                     [[buffer(12)]],
//       device float2 *neighborMap               [[buffer(13)]],  (output)
//       constant uint &particleCount             [[buffer(14)]],
//       uint particleId [[thread_position_in_grid]])

#include <cstdint>

#include "common/KernelArgs.h"

#include "Metal/MTLBuffer.hpp"
#include "Metal/MTLComputeCommandEncoder.hpp"

namespace Sibernetic {

struct FindNeighborsMetalArgs {
  MTL::Buffer *gridCellIndexFixedUp; // [[buffer(0)]]
  MTL::Buffer *sortedPosition;       // [[buffer(1)]]
  uint32_t gridCellCount;            // [[buffer(2)]]
  uint32_t gridCellsX;               // [[buffer(3)]]
  uint32_t gridCellsY;               // [[buffer(4)]]
  uint32_t gridCellsZ;               // [[buffer(5)]]
  float h;                           // [[buffer(6)]]
  float hashGridCellSize;            // [[buffer(7)]]
  float hashGridCellSizeInv;         // [[buffer(8)]]
  float simulationScale;             // [[buffer(9)]]
  float xmin;                        // [[buffer(10)]]
  float ymin;                        // [[buffer(11)]]
  float zmin;                        // [[buffer(12)]]
  MTL::Buffer *neighborMap;          // [[buffer(13)]] output
  uint32_t particleCount;            // [[buffer(14)]]

  void bind(MTL::ComputeCommandEncoder *enc) const {
    bindBuffer(enc, gridCellIndexFixedUp, 0);
    bindBuffer(enc, sortedPosition, 1);
    bindScalar(enc, gridCellCount, 2);
    bindScalar(enc, gridCellsX, 3);
    bindScalar(enc, gridCellsY, 4);
    bindScalar(enc, gridCellsZ, 5);
    bindScalar(enc, h, 6);
    bindScalar(enc, hashGridCellSize, 7);
    bindScalar(enc, hashGridCellSizeInv, 8);
    bindScalar(enc, simulationScale, 9);
    bindScalar(enc, xmin, 10);
    bindScalar(enc, ymin, 11);
    bindScalar(enc, zmin, 12);
    bindBuffer(enc, neighborMap, 13);
    bindScalar(enc, particleCount, 14);
  }
};

// ============ Constants ============
inline constexpr const char *kFindNeighborsKernelName = "findNeighbors";

} // namespace Sibernetic
