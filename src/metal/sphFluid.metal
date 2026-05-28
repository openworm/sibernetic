#include <metal_stdlib>
using namespace metal;

constant int kMaxNeighborCount = 32;
constant int kNoParticleId = -1;
// Each particle probes the 2×2×2 octant of cells overlapping its h-radius ball.
constant int kNeighborCellCount = 8;
// Maximum number of membrane triangles that share a single elastic particle.
constant int kMaxMembranesIncludingSameParticle = 7;
constant int kLiquidParticle = 1;
constant int kElasticParticle = 2;
constant int kBoundaryParticle = 3;

// Empirical surface tension scaling constant.
// See the surface tension force comment below for full attribution.
constant float kSurfaceTensionScale = -1.7e-09f;
// Empirical viscosity scale factor baked into this codebase.
// Not derived from the Müller (2003) viscosity formulation.
constant float kViscosityScale = 1.5f / 1000.0f;
// Default viscosity coefficient for fluid-fluid and fluid-boundary
// interactions.
constant float kViscosityCoeffDefault = 1.0e-4f;
// Reduced viscosity coefficient for worm-fluid interface interactions.
// 10x lower than default to allow worm to slide through the medium.
constant float kViscosityCoeffWormFluid = 1.0e-5f;
// Boundary friction damping factor (Ihmsen et al., 2010, eq. 12).
// Applied to velocity after removing the wall-normal component.
// 0.99 = 1% energy dissipation per correction.
constant float kBoundaryFriction = 0.99f;

// Particle type ranges are encoded per particle in configuration files
// ([position] section, 4th column) and loaded into originalPosition[].w.
constant float kWormTypeMin = 2.05f;
constant float kWormTypeMax = 2.25f;
constant float kFluidTypeMin = 2.25f;
constant float kFluidTypeMax = 2.35f;

inline bool isWormType(float particleType) {
  return particleType > kWormTypeMin && particleType < kWormTypeMax;
}

inline bool isFluidType(float particleType) {
  return particleType > kFluidTypeMin && particleType < kFluidTypeMax;
}

// Shared Poly6 radial core term: (hScaled^2 - r^2)^3.
// Reusing this in density and surface tension is correct because both
// operators are built from the same compact-support Poly6 distance weight;
// they differ in how this scalar is used (density sums scalar mass weights,
// while surface tension multiplies by a direction vector and a different
// physical coefficient).
inline float poly6Core(float distanceSquared, float hScaledSquared) {
  const float delta = hScaledSquared - distanceSquared;
  return delta * delta * delta;
}

// neighborMap: A packed table of neighbor information for SPH force kernels.
// Structure: float2 array with kMaxNeighborCount entries per particle.
//   Indexing: neighborMap[sortedParticleId * kMaxNeighborCount + slotIndex]
//   Each float2 tuple is:
//     .x = neighbor sorted particle ID (or kNoParticleId=-1 if empty)
//     .y = distance to neighbor in simulation units (or -1.0f if empty)
//
// Built by findNeighbors kernel:
//   - Probes neighboring spatial hash cells to find particles within h-radius
//   - Maintains a max-heap of kMaxNeighborCount closest neighbors
//   - Empty slots at the end have .x = kNoParticleId, .y = -1.0f
//
// Consumed by SPH kernels (pcisph_computeDensity,
// pcisph_computeForcesAndInitPressure):
//   - Iterate over slots; skip when .x == kNoParticleId
//   - Use .x to index into sorted position/velocity/density buffers
//   - Use .y for distance-dependent kernel evaluations (cutoff check at
//   hScaled)

// Flattens 3D grid coordinates (x,y,z) into a linear cell id.
inline int cellId(int3 cellCoordinates, uint gridCellsX, uint gridCellsY) {
  return cellCoordinates.x + cellCoordinates.y * static_cast<int>(gridCellsX) +
         cellCoordinates.z * static_cast<int>(gridCellsX) *
             static_cast<int>(gridCellsY);
}

// Converts a particle's world position into 3D grid cell coordinates.
// Divides each position component by the inverse cell size (effectively scaling
// by 1/cellSize to get which cell the particle belongs to).
// Only .x/.y/.z are used; .w is ignored, so this works for both
// originalPosition (where .w may be particle type) and sortedPosition
// (where .w may be cell ID).
inline int3 cellCoordinates(float4 position, float hashGridCellSizeInv) {
  return int3(static_cast<int>(position.x * hashGridCellSizeInv),
              static_cast<int>(position.y * hashGridCellSizeInv),
              static_cast<int>(position.z * hashGridCellSizeInv));
}

// Clears the neighborMap buffer to "no neighbor" sentinel values.
// Each thread handles one particle's kMaxNeighborCount float2 entries,
// writing (-1, -1) to every slot. Uses float4 writes for efficiency:
// each float4(-1) covers two float2 entries, so 16 writes clear 32 slots.
kernel void clearBuffers(device float4 *neighborMap [[buffer(0)]],
                         constant uint &particleCount [[buffer(1)]],
                         uint id [[thread_position_in_grid]]) {
  if (id >= particleCount)
    return;

  const int outIdx = (static_cast<int>(id) * kMaxNeighborCount) >> 1;
  const float4 noNeighbor = float4(-1.0f);

  for (int i = 0; i < kMaxNeighborCount / 2; ++i) {
    neighborMap[outIdx + i] = noNeighbor;
  }
}

// Zeros the delta accumulator region [N..2N) of position and velocity buffers.
// These regions store per-particle membrane interaction deltas that must be
// cleared at the start of each timestep.
kernel void clearMembraneBuffers(device float4 *position [[buffer(0)]],
                                 device float4 *velocity [[buffer(1)]],
                                 constant uint &particleCount [[buffer(2)]],
                                 uint id [[thread_position_in_grid]]) {
  if (id >= particleCount)
    return;

  position[particleCount + id] = float4(0.0f);
  velocity[particleCount + id] = float4(0.0f);
}

// Computes each particle's spatial hash cell and writes (cellId, serialId)
// into sortedCellAndSerialId. serialId is the original particle index in
// originalPosition[].
kernel void hashParticles(const device float4 *originalPosition [[buffer(0)]],
                          constant uint &gridCellsX [[buffer(1)]],
                          constant uint &gridCellsY [[buffer(2)]],
                          constant uint &gridCellsZ [[buffer(3)]],
                          constant float &hashGridCellSizeInv [[buffer(4)]],
                          constant float &xmin [[buffer(5)]],
                          constant float &ymin [[buffer(6)]],
                          constant float &zmin [[buffer(7)]],
                          device uint2 *sortedCellAndSerialId [[buffer(8)]],
                          constant uint &particleCount [[buffer(9)]],
                          uint serialId [[thread_position_in_grid]]) {
  if (serialId >= particleCount) {
    return;
  }

  (void)gridCellsZ;
  (void)xmin;
  (void)ymin;
  (void)zmin;

  const float4 p = originalPosition[serialId];
  int3 cellCoordinates;
  cellCoordinates.x = static_cast<int>(p.x * hashGridCellSizeInv);
  cellCoordinates.y = static_cast<int>(p.y * hashGridCellSizeInv);
  cellCoordinates.z = static_cast<int>(p.z * hashGridCellSizeInv);

  // Keep low 24 bits to match existing OpenCL behavior.
  // This limits the simulation to at most 2^24 = 16,777,216 grid cells
  // (e.g. a 256x256x256 grid). Typical Sibernetic runs use ~32^3-64^3
  // cells, so this ceiling is never approached in practice.
  const int cell = cellId(cellCoordinates, gridCellsX, gridCellsY) & 0x00ffffff;
  sortedCellAndSerialId[serialId] = uint2(static_cast<uint>(cell), serialId);
}

// After sortedCellAndSerialId is sorted by cell id, position/velocity are still
// in original serial order. This pass gathers them into sortedPosition/
// sortedVelocity and builds a reverse map (serialId -> sorted index).
//
// Example:
//   sorted sortedCellAndSerialId: [(1,1), (1,2), (2,3), (3,0), ...]
//   position:             [pos_0, pos_1, pos_2, pos_3, ...]
//   sortedPosition:       [pos_1, pos_2, pos_3, pos_0, ...]
kernel void sortPostPass(const device uint2 *sortedCellAndSerialId
                         [[buffer(0)]],
                         device uint *sortedParticleIdBySerialId [[buffer(1)]],
                         const device float4 *originalPosition [[buffer(2)]],
                         const device float4 *velocity [[buffer(3)]],
                         device float4 *sortedPosition [[buffer(4)]],
                         device float4 *sortedVelocity [[buffer(5)]],
                         constant uint &particleCount [[buffer(6)]],
                         uint sortedParticleId [[thread_position_in_grid]]) {
  if (sortedParticleId >= particleCount) {
    return;
  }

  const uint2 spi = sortedCellAndSerialId[sortedParticleId];
  const uint serialId = spi.y; // PI_SERIAL_ID
  const uint cellId = spi.x;   // PI_CELL_ID

  float4 pos = originalPosition[serialId];
  // Store cell ID in sorted position's .w component for neighbor lookups.
  // (Original position[].w still holds particle type; only sortedPosition[]
  // gets cell ID overwritten here.)
  pos.w = float(cellId); // POSITION_CELL_ID

  sortedPosition[sortedParticleId] = pos;
  sortedVelocity[sortedParticleId] = velocity[serialId];
  sortedParticleIdBySerialId[serialId] = sortedParticleId;
}

// For each cell id, finds the first index in sorted sortedCellAndSerialId whose
// particle belongs to that cell. Empty cells get UINT_MAX; slot gridCellCount
// stores PARTICLE_COUNT as the end sentinel.
// TODO(weiweng): Revisit this after full Metal kernel port. Consider replacing
// per-cell binary search with a linear boundary-scan build:
// 1) initialize gridCellIndex to UINT_MAX,
// 2) launch one thread per sorted particle index,
// 3) write start index when cell id changes,
// 4) set gridCellIndex[gridCellCount] = particleCount.
kernel void indexx(const device uint2 *sortedCellAndSerialId [[buffer(0)]],
                   constant uint &gridCellCount [[buffer(1)]],
                   device uint *gridCellIndex [[buffer(2)]],
                   constant uint &particleCount [[buffer(3)]],
                   uint targetCellId [[thread_position_in_grid]]) {
  if (targetCellId > gridCellCount) {
    return;
  }
  if (targetCellId == gridCellCount) {
    gridCellIndex[targetCellId] = particleCount;
    return;
  }
  if (targetCellId == 0) {
    gridCellIndex[targetCellId] = 0;
    return;
  }

  int low = 0;
  int high = static_cast<int>(particleCount) - 1;
  int cellIndex = -1;
  while (low <= high) {
    const int idx = low + (high - low) / 2;
    const int sampleCellId = static_cast<int>(sortedCellAndSerialId[idx].x);

    if (sampleCellId < static_cast<int>(targetCellId)) {
      low = idx + 1;
    } else if (sampleCellId > static_cast<int>(targetCellId)) {
      high = idx - 1;
    } else {
      // Found a match; check if it's the leftmost occurrence.
      // Guard idx-1 access: when idx==0 there is no previous element,
      // so this is trivially the left boundary.
      const bool isLeftBoundary =
          (idx == 0) ||
          (static_cast<int>(sortedCellAndSerialId[idx - 1].x) < sampleCellId);
      if (isLeftBoundary) {
        cellIndex = idx;
        break;
      }
      high = idx - 1; // keep searching left for the first occurrence
    }
  }

  // cellIndex == -1 when not found; cast to uint gives 0xFFFFFFFF == UINT_MAX.
  gridCellIndex[targetCellId] = static_cast<uint>(cellIndex);
}

// Linear scan for the index of the farthest (largest distance) neighbor.
// O(kMaxNeighborCount) but branch-free and GPU-friendly — uniform control flow
// across all threads in a SIMD group, unlike a heap which causes divergence.
inline int getMaxIndex(thread float *closestDistances) {
  int result = 0;
  float maxDist = closestDistances[0];
  for (int i = 1; i < kMaxNeighborCount; ++i) {
    if (closestDistances[i] > maxDist) {
      maxDist = closestDistances[i];
      result = i;
    }
  }
  return result;
}

// Scans all particles in probeCellId and maintains closestDistances/
// closestIndexes using a linear fill-then-scan strategy (matching the OpenCL
// kernel). The first kMaxNeighborCount candidates fill slots sequentially at
// O(1) each; once full, each new closer candidate replaces the current farthest
// entry found via getMaxIndex (O(kMaxNeighborCount) linear scan). This avoids
// the thread divergence caused by a heap-based approach.
inline int searchNeighborsInCell(
    int probeCellId, const device uint *gridCellIndices,
    float4 myParticlePosition, int mySortedParticleId,
    const device float4 *sortedPosition, thread int *closestIndexes,
    thread float *closestDistances, int lastFarthest, thread int *foundCount) {
  const int baseSortedParticleId =
      static_cast<int>(gridCellIndices[probeCellId]);
  const int nextSortedParticleId =
      static_cast<int>(gridCellIndices[probeCellId + 1]);
  const int particleCountThisCell = nextSortedParticleId - baseSortedParticleId;
  int farthestNeighbor = lastFarthest;

  for (int i = 0; i < particleCountThisCell; ++i) {
    const int neighborSortedParticleId = baseSortedParticleId + i;
    if (mySortedParticleId == neighborSortedParticleId)
      continue;

    const float4 d =
        myParticlePosition - sortedPosition[neighborSortedParticleId];
    const float distanceSquared = d.x * d.x + d.y * d.y + d.z * d.z;

    if (distanceSquared <= closestDistances[farthestNeighbor]) {
      closestDistances[farthestNeighbor] = distanceSquared;
      closestIndexes[farthestNeighbor] = neighborSortedParticleId;
      if (*foundCount < kMaxNeighborCount - 1) {
        (*foundCount)++;
        farthestNeighbor = *foundCount;
      } else {
        farthestNeighbor = getMaxIndex(closestDistances);
      }
    }
  }
  return farthestNeighbor;
}

// Computes the linear cell ID of a neighbor cell by applying (deltaX, deltaY,
// deltaZ) offsets to the given cell ID, with periodic boundary wrapping.
inline int neighborCellId(int cell, int deltaX, int deltaY, int deltaZ,
                          uint gridCellsX, uint gridCellsY,
                          uint gridCellCount) {
  const int dx = deltaX;
  const int dy = deltaY * static_cast<int>(gridCellsX);
  const int dz =
      deltaZ * static_cast<int>(gridCellsX) * static_cast<int>(gridCellsY);
  int newCell = cell + dx + dy + dz;
  const int gridCellCountInt = static_cast<int>(gridCellCount);
  newCell = (newCell < 0) ? (newCell + gridCellCountInt) : newCell;
  newCell =
      (newCell >= gridCellCountInt) ? (newCell - gridCellCountInt) : newCell;
  return newCell;
}

kernel void findNeighbors(const device uint *gridCellIndicesFixedUp
                          [[buffer(0)]],
                          const device float4 *sortedPosition [[buffer(1)]],
                          constant uint &gridCellCount [[buffer(2)]],
                          constant uint &gridCellsX [[buffer(3)]],
                          constant uint &gridCellsY [[buffer(4)]],
                          constant uint &gridCellsZ [[buffer(5)]],
                          constant float &h [[buffer(6)]],
                          constant float &hashGridCellSize [[buffer(7)]],
                          constant float &hashGridCellSizeInv [[buffer(8)]],
                          constant float &simulationScale [[buffer(9)]],
                          constant float &xmin [[buffer(10)]],
                          constant float &ymin [[buffer(11)]],
                          constant float &zmin [[buffer(12)]],
                          device float2 *neighborMap [[buffer(13)]],
                          constant uint &particleCount [[buffer(14)]],
                          uint sortedParticleId [[thread_position_in_grid]]) {
  if (sortedParticleId >= particleCount) {
    return;
  }

  (void)gridCellsZ;

  // "FixedUp" means empty-cell holes were filled in a post-pass so each cell
  // has a valid start pointer. Empty cells satisfy start[c] == start[c + 1],
  // producing a zero-length [start, end) range during neighbor scans.
  const device uint *gridCellIndices = gridCellIndicesFixedUp;
  const float4 position = sortedPosition[sortedParticleId];
  const int particleCellId = static_cast<int>(position.w) & 0x00ffffff;

  const float rThresholdSquared = h * h;
  float closestDistances[kMaxNeighborCount];
  int closestIndexes[kMaxNeighborCount];
  int foundCount = 0;
  for (int k = 0; k < kMaxNeighborCount; ++k) {
    closestDistances[k] = rThresholdSquared;
    closestIndexes[k] = kNoParticleId;
  }

  // Determine which neighboring cells to probe based on position inside cell.
  const float4 p = position - float4(xmin, ymin, zmin, 0.0f);
  const int3 cellCoords = cellCoordinates(position, hashGridCellSizeInv);
  const float3 cellMinCorner =
      float3(cellCoords.x * hashGridCellSize, cellCoords.y * hashGridCellSize,
             cellCoords.z * hashGridCellSize);
  const bool3 isInLowerHalf =
      bool3((p.x - cellMinCorner.x) < h, (p.y - cellMinCorner.y) < h,
            (p.z - cellMinCorner.z) < h);
  // Match OpenCL semantics: lower-half probes negative direction (-1),
  // upper-half probes positive direction (+1).
  const int3 delta = int3(isInLowerHalf.x ? -1 : 1, isInLowerHalf.y ? -1 : 1,
                          isInLowerHalf.z ? -1 : 1);

  // Build the kNeighborCellCount neighbor cell IDs. Each entry covers one
  // combination of (delta.x or 0, delta.y or 0, delta.z or 0), selected by
  // the 3 low bits of i: bit 0 → x, bit 1 → y, bit 2 → z.
  int neighborCellIds[kNeighborCellCount];
  for (int i = 0; i < kNeighborCellCount; ++i) {
    const int dx = (i & 1) ? delta.x : 0;
    const int dy = (i & 2) ? delta.y : 0;
    const int dz = (i & 4) ? delta.z : 0;
    neighborCellIds[i] = neighborCellId(particleCellId, dx, dy, dz, gridCellsX,
                                        gridCellsY, gridCellCount);
  }

  int lastFarthest = 0;
  for (int i = 0; i < kNeighborCellCount; ++i) {
    lastFarthest = searchNeighborsInCell(
        neighborCellIds[i], gridCellIndices, position,
        static_cast<int>(sortedParticleId), sortedPosition, closestIndexes,
        closestDistances, lastFarthest, &foundCount);
  }

  for (int i = 0; i < kMaxNeighborCount; ++i) {
    float2 neighborData;
    neighborData.x = static_cast<float>(closestIndexes[i]);
    if (closestIndexes[i] >= 0) {
      // Positions are in hash-grid space; multiply by simulationScale to
      // convert the distance into simulation units for the SPH kernels.
      neighborData.y = sqrt(closestDistances[i]) * simulationScale;
    } else {
      neighborData.y = -1.0f;
    }
    neighborMap[sortedParticleId * kMaxNeighborCount + static_cast<uint>(i)] =
        neighborData;
  }
}

// Computes elastic spring forces and muscle activation forces for the worm
// body simulation.
//
// Thread index: elastic particle index (0..numOfElasticP-1), which also
// equals the particle's serial ID because elastic particles occupy the
// first numOfElasticP slots in the originalPosition array (loaded before
// liquid and boundary particles).
//
// "Elastic" = any particle with int(type) == 2 (ELASTIC_PARTICLE).
// This includes both worm body particles (subtype 2.05–2.25) and agar
// particles (other 2.x subtypes). The kernel distinguishes them via the
// isWormType() check to apply different stiffness coefficients.
// Excluded: liquid particles (type 1) and boundary particles (type 3).
//
// For each elastic particle, iterates its elasticConnectionsData entries
// (up to kMaxNeighborCount):
//   - Computes Hooke's law: F = -k * (r - r0) * direction
//   - Worm body pairs (both particles in [2.05, 2.25]) use full elasticity
//   - Agar/other pairs use 25% elasticity
//   - If the connection is a muscle (muscleId > 0) and activated,
//     adds contraction force: F += -activation * maxMuscleForce * direction
//
// OpenCL signature (sphFluid.cl, lines 674-748):
//   __kernel void pcisph_computeElasticForces(
//       __global float2 *neighborMap,            // arg 0  (unused)
//       __global float4 *sortedPosition,         // arg 1
//       __global float4 *sortedVelocity,         // arg 2  (unused)
//       __global float4 *acceleration,           // arg 3
//       __global uint   *particleIndexBack,      // arg 4
//       __global uint2  *particleIndex,          // arg 5
//       float max_muscle_force,                  // arg 6
//       float mass,                              // arg 7  (unused)
//       float simulationScale,                   // arg 8
//       uint numOfElasticP,                      // arg 9
//       __global float4 *elasticConnectionsData, // arg 10
//       uint PARTICLE_COUNT,                     // arg 11 (unused)
//       uint MUSCLE_COUNT,                       // arg 12
//       __global float *muscle_activation_signal,// arg 13
//       __global float4 *position,               // arg 14
//       float elasticityCoefficient              // arg 15
//   )
kernel void pcisph_computeElasticForces(
    const device float4 *sortedPosition [[buffer(0)]],
    device float4 *acceleration [[buffer(1)]],
    const device uint *sortedParticleIdBySerialId [[buffer(2)]],
    const device uint2 *sortedCellAndSerialId [[buffer(3)]],
    constant float &maxMuscleForce [[buffer(4)]],
    constant float &simulationScale [[buffer(5)]],
    constant uint &numOfElasticP [[buffer(6)]],
    const device float4 *elasticConnectionsData [[buffer(7)]],
    constant uint &muscleCount [[buffer(8)]],
    const device float *muscleActivationSignal [[buffer(9)]],
    const device float4 *originalPosition [[buffer(10)]],
    constant float &elasticityCoefficient [[buffer(11)]],
    uint elasticIndex [[thread_position_in_grid]]) {
  if (elasticIndex >= numOfElasticP)
    return;

  // elasticIndex is the serial elastic particle index. Look up its sorted ID.
  const int sortedParticleId =
      static_cast<int>(sortedParticleIdBySerialId[elasticIndex]);
  // sortedCellAndSerialId[sortedParticleId].y is the inverse of
  // sortedParticleIdBySerialId[elasticIndex], so it always equals
  // elasticIndex. No separate serialId variable needed.
  const int connectionBaseIdx =
      static_cast<int>(elasticIndex) * kMaxNeighborCount;

  for (int connectionSlot = 0; connectionSlot < kMaxNeighborCount;
       ++connectionSlot) {
    const int connectedSerialId = static_cast<int>(
        elasticConnectionsData[connectionBaseIdx + connectionSlot].x);
    if (connectedSerialId == kNoParticleId)
      break;

    const int neighborSortedParticleId =
        static_cast<int>(sortedParticleIdBySerialId[connectedSerialId]);
    const int neighborSerialId =
        static_cast<int>(sortedCellAndSerialId[neighborSortedParticleId].y);
    const float equilibriumDistance =
        elasticConnectionsData[connectionBaseIdx + connectionSlot].y;

    float4 displacement = (sortedPosition[sortedParticleId] -
                           sortedPosition[neighborSortedParticleId]) *
                          simulationScale;
    displacement.w = 0.0f;

    const float distanceToNeighbor = length(displacement.xyz);
    if (distanceToNeighbor == 0.0f)
      continue;

    const float distanceDelta = distanceToNeighbor - equilibriumDistance;
    const float4 direction = displacement / distanceToNeighbor;

    // Worm body particles: both endpoints in (2.05, 2.25) → full elasticity.
    // Otherwise (agar/other): 25% elasticity.
    float stiffness = elasticityCoefficient;
    if (!(isWormType(originalPosition[elasticIndex].w) &&
          isWormType(originalPosition[neighborSerialId].w))) {
      stiffness *= 0.25f;
    }
    acceleration[sortedParticleId] -= direction * distanceDelta * stiffness;

    // Muscle activation: elasticConnectionsData[].z is 1-indexed muscle ID.
    const int muscleId = static_cast<int>(
        elasticConnectionsData[connectionBaseIdx + connectionSlot].z);
    if (muscleId > 0 && muscleId <= static_cast<int>(muscleCount)) {
      const float activation = muscleActivationSignal[muscleId - 1];
      if (activation > 0.0f) {
        acceleration[sortedParticleId] -=
            direction * activation * maxMuscleForce;
      }
    }
  }
}

// Computes per-particle density by summing poly6 neighbor contributions
// within hScaled, then applying the hScaled^6 minimum-density floor.
//
// Arguments:
// buffer(0) neighborMap: packed table of 32 float2 entries per particle.
//   x = neighbor sorted particle id (or -1 sentinel), y = neighbor distance.
// buffer(1) massMultWpoly6Coefficient: mass * 315 / (64 * pi * hScaled^9),
//   where hScaled = h * simulationScale. Multiplying poly6Sum by this converts
//   the raw kernel sum into density (mass/volume).
// buffer(2) hScaled2: (h * simulationScale)^2, matching the scale of distances
// in neighborMap.
// buffer(3) rho: output density array, indexed by sorted particle id.
// buffer(4) sortedParticleIdBySerialId: map from serial (original) particle id
// -> sorted index. buffer(5) particleCount: number of particles to process.
// thread_position_in_grid serialId: serial (original) particle index for this
// thread.
kernel void pcisph_computeDensity(
    const device float2 *neighborMap [[buffer(0)]],
    constant float &massMultWpoly6Coefficient [[buffer(1)]],
    constant float &hScaled2 [[buffer(2)]], device float *rho [[buffer(3)]],
    const device uint *sortedParticleIdBySerialId [[buffer(4)]],
    constant uint &particleCount [[buffer(5)]],
    uint serialId [[thread_position_in_grid]]) {
  if (serialId >= particleCount) {
    return;
  }

  const uint sortedParticleId = sortedParticleIdBySerialId[serialId];
  const uint idx = sortedParticleId * static_cast<uint>(kMaxNeighborCount);
  // poly6NeighborContributionSum = sum_j (hScaled2 - r_j^2)^3 for all
  // neighbors j with r_j <
  // hScaled. Floored to the Poly6 self-contribution at r=0:
  //   poly6SelfContribution = (hScaled2 - 0)^3 = hScaled^6.
  // This prevents the density kernel sum from becoming zero when a particle
  // has no
  // neighbors, which would otherwise produce rho=0 and a numerically
  // explosive pressure correction.
  float poly6NeighborContributionSum = 0.0f;
  const float poly6SelfContribution = poly6Core(0.0f, hScaled2);

  for (int neighborSlot = 0; neighborSlot < kMaxNeighborCount; ++neighborSlot) {
    const float2 neighbor = neighborMap[idx + static_cast<uint>(neighborSlot)];
    if (static_cast<int>(neighbor.x) != kNoParticleId) {
      float distanceToNeighbor2 = neighbor.y;
      distanceToNeighbor2 *= distanceToNeighbor2;
      if (distanceToNeighbor2 < hScaled2) {
        poly6NeighborContributionSum +=
            poly6Core(distanceToNeighbor2, hScaled2);
      }
    }
  }

  if (poly6NeighborContributionSum < poly6SelfContribution) {
    poly6NeighborContributionSum = poly6SelfContribution;
  }
  rho[sortedParticleId] =
      poly6NeighborContributionSum * massMultWpoly6Coefficient;
}

// Writes the non-pressure acceleration and zeros the pressure-force slot and
// pressure scalar for sortedParticleId.
//
// The acceleration buffer is a single allocation of 2*particleCount entries
// split into two logical halves:
//   [0 .. N-1]   non-pressure acceleration (viscosity + gravity + surface
//                tension) — written as nonPressureAccel (float4(0) for
//                boundary particles, total_accel for fluid particles).
//   [N .. 2N-1]  pressure acceleration accumulated by the PCISPH
//                predict/correct iterations that follow this kernel — zeroed
//                here so each timestep starts from a clean slate.
// Both are written unconditionally because sorted slots are reassigned every
// timestep and may have held a different particle's data in the prior step.
inline void writeAccelerationAndInitPressure(device float4 *acceleration,
                                             device float *pressure,
                                             uint particleCount,
                                             uint sortedParticleId,
                                             float4 nonPressureAccel) {
  acceleration[sortedParticleId] = nonPressureAccel;
  acceleration[particleCount + sortedParticleId] = float4(0.0f);
  pressure[sortedParticleId] = 0.0f;
}

/// PCISPH: Compute viscosity + surface tension + gravity forces.
/// Initializes pressure to 0 and acceleration[N..2N] to 0.
kernel void pcisph_computeForcesAndInitPressure(
    const device float2 *neighborMap [[buffer(0)]],
    const device float *rho [[buffer(1)]], device float *pressure [[buffer(2)]],
    const device float4 *sortedPosition [[buffer(3)]],
    const device float4 *sortedVelocity [[buffer(4)]],
    device float4 *acceleration [[buffer(5)]],
    const device uint *sortedParticleIdBySerialId [[buffer(6)]],
    constant float &surfTensCoeff [[buffer(7)]],
    constant float &massMultLaplacianWviscosityCoeff [[buffer(8)]],
    constant float &hScaled [[buffer(9)]], constant float &mu [[buffer(10)]],
    constant float &gravitationalAccelerationX [[buffer(11)]],
    constant float &gravitationalAccelerationY [[buffer(12)]],
    constant float &gravitationalAccelerationZ [[buffer(13)]],
    const device float4 *originalPosition [[buffer(14)]],
    const device uint2 *sortedCellAndSerialId [[buffer(15)]],
    constant uint &particleCount [[buffer(16)]],
    constant float &mass [[buffer(17)]],
    uint serialId [[thread_position_in_grid]]) {
  if (serialId >= particleCount)
    return;

  uint sortedParticleId = sortedParticleIdBySerialId[serialId];

  // Boundary particles don't move.
  // Note: originalPosition[].w still holds the original particle type;
  // sortedPosition[].w was overwritten with cell ID in sortPostPass. Match the
  // OpenCL behavior by classifying via integer cast (i.e. floor/truncate for
  // positive values) instead of exact float equality. Configuration
  // files (such as demo1) store subtype tags in the fractional part (e.g. 3.1),
  // so exact comparison with 3.0 would miss boundary-like values that should
  // map to boundary class 3.
  if (int(originalPosition[serialId].w) == kBoundaryParticle) {
    // Boundary particles never move, but sortedParticleId changes every
    // timestep as the sort reassigns sorted slots. The slot this particle
    // occupies now may have held a non-boundary particle last timestep, so
    // its acceleration entry could contain stale non-zero values. Write zeros
    // unconditionally to own the output slot cleanly.
    writeAccelerationAndInitPressure(acceleration, pressure, particleCount,
                                     sortedParticleId, float4(0.0f));
    return;
  }

  int neighborMapOffset = sortedParticleId * kMaxNeighborCount;
  float hScaled2 = hScaled * hScaled;

  float4 accelViscosity = float4(0.0f);
  float4 accelSurfaceTension = float4(0.0f);

  // Loop through neighbors
  for (int neighborSlot = 0; neighborSlot < kMaxNeighborCount; ++neighborSlot) {
    float2 neighborEntry = neighborMap[neighborMapOffset + neighborSlot];
    int neighborSortedParticleId = int(neighborEntry.x);
    if (neighborSortedParticleId == kNoParticleId)
      continue;

    float distanceToNeighbor = neighborEntry.y;
    // findNeighbors admits only neighbors with distanceSquared <= h*h, so
    // distanceToNeighbor <= hScaled is guaranteed. The >= check handles the
    // exact-equality boundary case: at distanceToNeighbor == hScaled both
    // the viscosity term (hScaled - distanceToNeighbor = 0) and the surface
    // tension kernel ((hScaled2 - distanceToNeighbor2) = 0) contribute
    // nothing, so skipping saves the work.
    if (distanceToNeighbor >= hScaled)
      continue;

    float distanceToNeighbor2 = distanceToNeighbor * distanceToNeighbor;

    uint neighborSerialId = sortedCellAndSerialId[neighborSortedParticleId].y;
    // 1.0 for fluid neighbors, 0.0 for boundary neighbors. Boundary particles
    // are stationary walls whose velocity contribution to viscosity is zeroed
    // out by multiplying vj by this mask.
    float neighborFluidMask =
        (int(originalPosition[neighborSerialId].w) != kBoundaryParticle) ? 1.0f
                                                                         : 0.0f;

    float4 particleVelocity = sortedVelocity[sortedParticleId];
    float4 neighborVelocity = sortedVelocity[neighborSortedParticleId];
    // originalPosition[].w holds the particle type (unchanged by sortPostPass).
    float particleType = originalPosition[serialId].w;
    float neighborType = originalPosition[neighborSerialId].w;
    // Müller et al. (2003), written as a neighbor sum:
    //   a_i += sum_j [ (mu / rho_j) * (m_j / rho_i) * (v_j - v_i)
    //                  * laplacian_W_visc(r_ij, h) ]
    // where laplacian_W_visc(r, h) = (45 / (pi * h^6)) * (h - r).
    // This loop body computes one j-term at a time; the sum_j is realized by
    // iterating over neighbor slots and accumulating into accelViscosity.
    // In this implementation, massMultLaplacianWviscosityCoeff packs
    // m_j * (45 / (pi * hScaled^6)) (with constant particle mass), while
    // 1 / rho_i is applied separately in the finalize step after the loop.
    //
    // Here viscCoeff replaces mu / rho_j: mu (the true dynamic viscosity,
    // Pa·s) is passed as buffer(10) but is not used. Instead, viscCoeff is
    // an empirical per-interaction-type tuning constant that lumps together
    // mu and the neighbor density rho_j into a single scalar.
    float viscCoeff =
        kViscosityCoeffDefault; // default (fluid–fluid or fluid–boundary)

    // Worm body ([kWormTypeMin, kWormTypeMax]) <->
    // fluid medium ([kFluidTypeMin, kFluidTypeMax]): use lower
    // interface viscosity so the worm can slide through the medium.
    // If this interface viscosity is too high, drag becomes too strong and
    // locomotion is over-damped.
    bool isWormParticle = isWormType(particleType);
    bool isWormNeighbor = isWormType(neighborType);
    bool isFluidParticle = isFluidType(particleType);
    bool isFluidNeighbor = isFluidType(neighborType);

    if ((isWormParticle && isFluidNeighbor) ||
        (isFluidParticle && isWormNeighbor)) {
      viscCoeff =
          kViscosityCoeffWormFluid; // worm-fluid interface: 10x lower viscosity
    }

    // This loop accumulates only the pairwise velocity-shape part:
    //   sum_j [ viscCoeff_ij * (v_j - v_i) * (hScaled - r_ij) ]
    // where viscCoeff_ij is an empirical stand-in for (mu / rho_j).
    //
    // Compared with Muller viscosity:
    //   a_i += sum_j [ (mu/rho_j) * (m_j/rho_i) * (v_j-v_i)
    //                  * (45/(pi*h^6)) * (h-r_ij) ]
    // this line intentionally omits ((m_j*45)/(pi*hScaled^6)) and (1/rho_i);
    // both are applied after the neighbor loop via:
    //   accelViscosity *= (1.5/1000) * massMultLaplacianWviscosityCoeff /
    //                     rho_i
    // so rho_i and kernel/mass scaling are deferred to the finalize step.
    accelViscosity +=
        viscCoeff * (neighborVelocity * neighborFluidMask - particleVelocity) *
        (hScaled - distanceToNeighbor);

    // Surface tension force (legacy Sibernetic form).
    // Discrete pairwise form used here:
    //   a_i^surf += C_st * s_ij * (x_i - x_j),
    //   s_ij = (hScaled^2 - r_ij^2)^3,
    // where C_st = (-1.7e-09) * surfTensCoeff / mass (the /mass is applied
    // in the finalize step below).
    // Literature attribution: this term is historically tagged in the OpenCL
    // implementation as "formula (16) [5]" from Becker and Teschner,
    // "Weakly compressible SPH for free surface flows" (SCA 2007,
    // DOI: 10.2312/SCA/SCA07/209-218). The -1.7e-09 scalar is a project-
    // specific empirical scaling, not a universal constant from the paper.
    const float surfKern = poly6Core(distanceToNeighbor2, hScaled2);
    // sortedPosition[].w stores cell ID (not position), so the .w part of this
    // subtraction is not physically meaningful; we zero accelSurfaceTension.w
    // after the loop.
    accelSurfaceTension += kSurfaceTensionScale * surfTensCoeff * surfKern *
                           (sortedPosition[sortedParticleId] -
                            sortedPosition[neighborSortedParticleId]);
  }

  // Finalize surface tension
  accelSurfaceTension.w = 0.0f;
  accelSurfaceTension /= mass;
  // Finalize viscosity
  // 1.5f/1000.0f is an empirical scale factor from this codebase,
  // not a coefficient from the Müller viscosity derivation.
  accelViscosity *= kViscosityScale * massMultLaplacianWviscosityCoeff /
                    rho[sortedParticleId];

  // Total acceleration = viscosity + gravitational acceleration +
  // surface tension
  float4 totalAccel = accelViscosity;
  totalAccel += float4(gravitationalAccelerationX, gravitationalAccelerationY,
                       gravitationalAccelerationZ, 0.0f);
  totalAccel += accelSurfaceTension;

  writeAccelerationAndInitPressure(acceleration, pressure, particleCount,
                                   sortedParticleId, totalAccel);
}

// Boundary handling helper for pcisph_predictPositions.
// Pushes the predicted position away from nearby boundary particles using
// weighted normal vectors (Ihmsen et al., 2010, Sec. 3.2).
// When correctVelocity is true, also removes the normal component of velocity
// and applies friction (eps = 0.99).
//
// Limitation: this is a soft repulsion field, not a hard constraint. The
// boundary weight w_c(i,b) = max(0, (r0-d)/r0) ramps linearly from 1 at d=0
// to 0 at d=r0. If a particle's velocity is large enough to jump past the
// entire r0 influence zone in one timestep, all weights clamp to zero and no
// correction is applied — the particle tunnels through the wall. The
// simulation relies on deltaTime being small enough that per-step displacement
// stays within the r0 zone.
inline void computeInteractionWithBoundaryParticles(
    int sortedParticleId, float r0, const device float2 *neighborMap,
    const device uint *sortedParticleIdBySerialId,
    const device uint2 *sortedCellAndSerialId,
    const device float4 *originalPosition, const device float4 *velocity,
    thread float4 *predictedPosition, bool correctVelocity,
    thread float4 *predictedVelocity, uint particleCount) {
  const int neighborMapOffset = sortedParticleId * kMaxNeighborCount;
  float4 n_c_i = float4(0.0f);
  float w_c_ib_sum = 0.0f;
  float w_c_ib_second_sum = 0.0f;

  for (int neighborSlot = 0; neighborSlot < kMaxNeighborCount; ++neighborSlot) {
    const int neighborSortedParticleId =
        static_cast<int>(neighborMap[neighborMapOffset + neighborSlot].x);
    if (neighborSortedParticleId == kNoParticleId)
      continue;

    const uint neighborSerialId =
        sortedCellAndSerialId[neighborSortedParticleId].y;
    if (static_cast<int>(originalPosition[neighborSerialId].w) !=
        kBoundaryParticle)
      continue;

    const float4 d = (*predictedPosition) - originalPosition[neighborSerialId];
    float distanceToNeighbor = sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
    const float w_c_ib =
        max(0.0f, (r0 - distanceToNeighbor) / r0); // Ihmsen (10)
    // Boundary particles store their outward normal in velocity[].
    const float4 n_b = velocity[neighborSerialId]; // Ihmsen (9)
    n_c_i += n_b * w_c_ib;
    w_c_ib_sum += w_c_ib; // Ihmsen (11) sum #1
    w_c_ib_second_sum +=
        w_c_ib * (r0 - distanceToNeighbor); // Ihmsen (11) sum #2
  }

  // Only .xyz of n_c_i is physically meaningful; .w accumulates from
  // velocity[].w of boundary particles (unused padding). Use .xyz explicitly
  // so correctness does not depend on .w being 0.
  float n_c_i_length = dot(n_c_i.xyz, n_c_i.xyz);
  if (n_c_i_length != 0.0f) {
    n_c_i_length = sqrt(n_c_i_length);
    // Ihmsen (11): position correction =
    //   normalized contact normal * weighted average penetration depth.
    //   n_c_i / n_c_i_length = unit contact normal direction
    //   w_c_ib_second_sum / w_c_ib_sum = weighted avg of (r0 - d) over
    //   boundary neighbors
    const float4 deltaPos =
        ((n_c_i / n_c_i_length) * w_c_ib_second_sum) / w_c_ib_sum;
    (*predictedPosition).xyz += deltaPos.xyz;

    if (correctVelocity) {
      // Projection of velocity onto the (unnormalized) contact normal.
      // Negative means the particle is moving into the wall; only then
      // do we strip the wall-normal component and apply friction.
      const float velocityNormalProjection =
          dot(n_c_i.xyz, (*predictedVelocity).xyz);
      if (velocityNormalProjection < 0.0f) {
        (*predictedVelocity).xyz -= n_c_i.xyz * velocityNormalProjection;
        (*predictedVelocity) *= kBoundaryFriction; // Ihmsen (12)
      }
    }
  }
}

// PCISPH position prediction: semi-implicit Euler step.
// Reads current sorted position/velocity and acceleration (non-pressure from
// the first half + pressure from the second half),
// writes predicted position into sortedPosition[particleCount +
// sortedParticleId]. The second half uses the same sorted index as the first
// half for 1:1 correspondence, but is not spatially sorted by cell ID —
// particles may have moved to different cells after integration. Re-sorting
// happens at the start of the next timestep (hashParticles → sort →
// sortPostPass). Boundary particles just copy their current position unchanged.
kernel void pcisph_predictPositions(
    device float4 *acceleration [[buffer(0)]],
    device float4 *sortedPosition [[buffer(1)]],
    const device float4 *sortedVelocity [[buffer(2)]],
    const device uint2 *sortedCellAndSerialId [[buffer(3)]],
    const device uint *sortedParticleIdBySerialId [[buffer(4)]],
    constant float &gravitationalAccelerationX [[buffer(5)]],
    constant float &gravitationalAccelerationY [[buffer(6)]],
    constant float &gravitationalAccelerationZ [[buffer(7)]],
    constant float &simulationScaleInv [[buffer(8)]],
    constant float &deltaTime [[buffer(9)]],
    // originalPosition: unsorted originalPosition array indexed by serialId.
    // .w holds particle type (e.g. kBoundaryParticle); used for type checks
    // and boundary interaction geometry.
    const device float4 *originalPosition [[buffer(10)]],
    // velocity: unsorted velocity array indexed by serialId.
    // For fluid particles: .xyz = velocity, .w = 0 (unused padding).
    // For boundary particles: .xyz = outward surface normal (used by
    // computeInteractionWithBoundaryParticles for Ihmsen wall repulsion).
    const device float4 *velocity [[buffer(11)]],
    // r0: boundary interaction radius (= h * 0.5). Cutoff distance for the
    // Ihmsen et al. (2010) boundary weight function w_c_ib = max(0, (r0-d)/r0).
    constant float &r0 [[buffer(12)]],
    // neighborMap: packed neighbor table, kMaxNeighborCount float2 entries per
    // particle. Each entry: (.x = neighbor sorted particle ID, .y = distance).
    const device float2 *neighborMap [[buffer(13)]],
    constant uint &particleCount [[buffer(14)]],
    uint serialId [[thread_position_in_grid]]) {
  if (serialId >= particleCount)
    return;

  const uint sortedParticleId = sortedParticleIdBySerialId[serialId];

  float4 currentPosition = sortedPosition[sortedParticleId];

  // Boundary particles are stationary; just copy current position.
  if (static_cast<int>(originalPosition[serialId].w) == kBoundaryParticle) {
    sortedPosition[particleCount + sortedParticleId] = currentPosition;
    return;
  }

  // Current non-pressure + pressure accelerations.
  // Only .xyz is physically meaningful; .w is not used.
  float4 currentAcceleration = acceleration[sortedParticleId] +
                               acceleration[particleCount + sortedParticleId];

  float4 currentVelocity = sortedVelocity[sortedParticleId];
  // Semi-implicit Euler integration.
  float4 predictedVelocity = currentVelocity + deltaTime * currentAcceleration;
  // Velocity is in simulation space, but position is in grid space.
  // Multiply deltaTime by simulationScaleInv (= 1/simulationScale) to convert
  // the simulation-space displacement into grid-space units.
  const float posDeltaTime = deltaTime * simulationScaleInv;
  float4 predictedPosition = currentPosition + posDeltaTime * predictedVelocity;

  // Position-only boundary correction (correctVelocity=false). Velocity is not
  // adjusted here because this kernel runs inside the iterative PCISPH
  // pressure loop — the prediction is tentative and may be recomputed.
  // The velocity friction correction (Ihmsen eq. 12) is applied once in
  // pcisph_correctPosition after the pressure loop converges
  // (correctVelocity=true).
  computeInteractionWithBoundaryParticles(
      static_cast<int>(sortedParticleId), r0, neighborMap,
      sortedParticleIdBySerialId, sortedCellAndSerialId, originalPosition,
      velocity, &predictedPosition, false, &predictedVelocity, particleCount);

  sortedPosition[particleCount + sortedParticleId] = predictedPosition;
}

// Computes predicted density from predicted positions for the PCISPH pressure
// correction loop. Reads predicted positions from sortedPosition[N..2N)
// (written by pcisph_predictPositions), accumulates Poly6 neighbor
// contributions scaled to simulation space, and writes the result to
// rho[N..2N).
kernel void pcisph_predictDensity(
    const device float2 *neighborMap [[buffer(0)]],
    const device uint *sortedParticleIdBySerialId [[buffer(1)]],
    constant float &massMultWpoly6Coefficient [[buffer(2)]],
    constant float &h [[buffer(3)]],
    // TODO(weiweng): Remove restDensity once the host no longer binds it.
    // It exists only to keep buffer indices consistent with the OpenCL
    // kernel signature; this kernel never reads it.
    constant float &restDensity [[buffer(4)]],
    constant float &simulationScale [[buffer(5)]],
    const device float4 *sortedPosition [[buffer(6)]],
    device float *rho [[buffer(7)]], constant uint &particleCount [[buffer(8)]],
    uint serialId [[thread_position_in_grid]]) {
  if (serialId >= particleCount)
    return;

  (void)restDensity;

  const uint sortedParticleId = sortedParticleIdBySerialId[serialId];
  const int neighborMapOffset =
      static_cast<int>(sortedParticleId) * kMaxNeighborCount;

  const float hScaled = h * simulationScale;
  const float hScaled2 = hScaled * hScaled;
  const float simulationScale2 = simulationScale * simulationScale;

  float densityAccum = 0.0f;
  const float poly6SelfContribution = poly6Core(0.0f, hScaled2);

  for (int neighborSlot = 0; neighborSlot < kMaxNeighborCount; ++neighborSlot) {
    const int neighborSortedParticleId =
        static_cast<int>(neighborMap[neighborMapOffset + neighborSlot].x);
    if (neighborSortedParticleId == kNoParticleId)
      continue;

    // Displacement between predicted positions in unscaled grid space.
    // sortedPosition[N..2N) holds predicted positions written by
    // pcisph_predictPositions.
    const float4 predictedDisplacement =
        sortedPosition[particleCount + sortedParticleId] -
        sortedPosition[particleCount +
                       static_cast<uint>(neighborSortedParticleId)];
    const float distSquaredScaled =
        dot(predictedDisplacement.xyz, predictedDisplacement.xyz) *
        simulationScale2;

    if (distSquaredScaled < hScaled2) {
      densityAccum += poly6Core(distSquaredScaled, hScaled2);
    }
  }

  // Floor to the Poly6 self-contribution at r=0 in simulation space.
  if (densityAccum < poly6SelfContribution) {
    densityAccum = poly6SelfContribution;
  }

  rho[particleCount + sortedParticleId] =
      densityAccum * massMultWpoly6Coefficient;
}

// PCISPH pressure correction: reads predicted density from rho[N..2N),
// computes the density error against the reference density rho0, scales by the
// precomputed correction factor delta, clamps to non-negative, and accumulates
// into pressure[0..N).
//
// This kernel runs once per PCISPH iteration, after pcisph_predictDensity has
// written the predicted density into rho[N..2N).
//
// OpenCL signature (sphFluid.cl, lines 929-952):
//   __kernel void pcisph_correctPressure(
//       __global uint *particleIndexBack,  // arg 0
//       float rho0,                        // arg 1
//       __global float *pressure,          // arg 2
//       __global float *rho,               // arg 3
//       float delta,                       // arg 4
//       uint PARTICLE_COUNT                // arg 5
//   )
kernel void pcisph_correctPressure(
    const device uint *sortedParticleIdBySerialId [[buffer(0)]],
    constant float &restDensity [[buffer(1)]],
    device float *pressure [[buffer(2)]], const device float *rho [[buffer(3)]],
    // delta: precomputed PCISPH pressure-correction scaling factor.
    // Solenthaler & Pajarola (2009) eq. 8 / Solenthaler dissertation eq. 3.6:
    //   delta = 1 / (beta * gradWspikyCoeff^2 * (|sum_j gradW_j|^2 + sum_j
    //   |gradW_j|^2))
    // where beta = 2 * dt^2 * mass^2 / rho0^2.
    // Computed once at startup by owConfigProperty::calcDelta() using a
    // synthetic 32-neighbor filled neighborhood.
    constant float &delta [[buffer(4)]],
    constant uint &particleCount [[buffer(5)]],
    uint serialId [[thread_position_in_grid]]) {
  if (serialId >= particleCount)
    return;

  const uint sortedParticleId = sortedParticleIdBySerialId[serialId];

  // Density error: predicted density minus reference density.
  const float rhoErr = rho[particleCount + sortedParticleId] - restDensity;

  // Pressure correction scaled by delta; clamp to non-negative.
  float pCorr = rhoErr * delta;
  if (pCorr < 0.0f)
    pCorr = 0.0f;

  // Accumulate correction into existing pressure.
  pressure[sortedParticleId] += pCorr;
}

// PCISPH pressure-gradient force acceleration.
//
// For each fluid particle, sums the pressure force over neighbors using the
// Spiky kernel gradient (Solenthaler & Pajarola 2009, formula 5):
//   a_i = -(m / rho_i) * sum_j (p_i + p_j) / (2 * rho_j)
//         * grad_W_spiky(r_ij)
// where grad_W_spiky ∝ -(hScaled - r_ij)^2 * (vec{x_i} - vec{x_j}) / r_ij.
//   r_ij = |vec{x_i} - vec{x_j}| is the scalar distance,
//   (vec{x_i} - vec{x_j}) is the displacement vector from j to i,
//   so (vec{x_i} - vec{x_j}) / r_ij is the unit direction from j toward i.
//
// Close-range correction: when r_ij < hScaled/4, the pressure term is replaced
// with a repulsive correction based on (restDensity * delta) to prevent
// particle overlap at very small separations.
//
// Reads predicted density from rho[N..2N).
// Writes pressure acceleration to acceleration[N..2N).
// Boundary particles (originalPosition[].w == 3) get zero acceleration.
kernel void pcisph_computePressureForceAcceleration(
    const device float2 *neighborMap [[buffer(0)]],
    const device float *pressure [[buffer(1)]],
    const device float *rho [[buffer(2)]],
    const device float4 *sortedPosition [[buffer(3)]],
    const device uint *sortedParticleIdBySerialId [[buffer(4)]],
    // delta: precomputed PCISPH pressure-correction scaling factor.
    // See pcisph_correctPressure for derivation.
    constant float &delta [[buffer(5)]],
    constant float &massMultGradWspikyCoefficient [[buffer(6)]],
    constant float &h [[buffer(7)]],
    constant float &simulationScale [[buffer(8)]],
    constant float &restDensity [[buffer(9)]],
    device float4 *acceleration [[buffer(10)]],
    const device float4 *originalPosition [[buffer(11)]],
    const device uint2 *sortedCellAndSerialId [[buffer(12)]],
    constant uint &particleCount [[buffer(13)]],
    uint serialId [[thread_position_in_grid]]) {
  if (serialId >= particleCount)
    return;

  const uint sortedParticleId = sortedParticleIdBySerialId[serialId];

  if (static_cast<int>(originalPosition[serialId].w) == kBoundaryParticle) {
    acceleration[particleCount + sortedParticleId] = float4(0.0f);
    return;
  }

  const int neighborMapOffset =
      static_cast<int>(sortedParticleId) * kMaxNeighborCount;
  const float hScaled = h * simulationScale;
  const float pressureI = pressure[sortedParticleId];
  const float rhoI = rho[particleCount + sortedParticleId];

  float4 result = float4(0.0f);

  for (int neighborSlot = 0; neighborSlot < kMaxNeighborCount; ++neighborSlot) {
    const float2 neighborEntry = neighborMap[neighborMapOffset + neighborSlot];
    const int neighborSortedParticleId = static_cast<int>(neighborEntry.x);
    if (neighborSortedParticleId == kNoParticleId)
      continue;

    float distanceToNeighbor = neighborEntry.y;
    if (distanceToNeighbor >= hScaled)
      continue;

    // Solenthaler (2009) formula (5): pressure force contribution.
    float value = -(hScaled - distanceToNeighbor) *
                  (hScaled - distanceToNeighbor) * 0.5f *
                  (pressureI + pressure[neighborSortedParticleId]) /
                  rho[particleCount + neighborSortedParticleId];

    // Direction vector in simulation space.
    float4 direction = (sortedPosition[sortedParticleId] -
                        sortedPosition[neighborSortedParticleId]) *
                       simulationScale;

    // Close-range correction: when particles are very close
    // (r < hScaled/4), replace the pressure term with a repulsive
    // correction to prevent particle overlap.
    if (distanceToNeighbor < hScaled * 0.25f) {
      value = -(hScaled * 0.25f - distanceToNeighbor) *
              (hScaled * 0.25f - distanceToNeighbor) * 0.5f *
              (restDensity * delta) /
              rho[particleCount + neighborSortedParticleId];
    }

    result += value * direction / distanceToNeighbor;
  }

  result *= massMultGradWspikyCoefficient / rhoI;
  acceleration[particleCount + sortedParticleId] = result;
}

// Leapfrog / semi-implicit Euler time integration.
//
// timestepIndex == 0: store combined acceleration for the first leapfrog step.
// mode == 0 (positions): leapfrog position update x(t+dt) = x(t) + v*dt +
// a*dt^2/2. mode == 1 (velocities): leapfrog velocity update v(t+dt) = v(t) +
// (a(t)+a(t+dt))*dt/2,
//   with boundary interaction correction.
// mode == 2 (semi-implicit Euler): v(t+dt) = v(t) + a(t+dt)*dt,
//   x(t+dt) = x(t) + v(t+dt)*dt, with boundary interaction correction.
kernel void pcisph_integrate(device float4 *acceleration [[buffer(0)]],
                             device float4 *sortedPosition [[buffer(1)]],
                             const device float4 *sortedVelocity [[buffer(2)]],
                             const device uint2 *sortedCellAndSerialId
                             [[buffer(3)]],
                             const device uint *sortedParticleIdBySerialId
                             [[buffer(4)]],
                             constant float &simulationScaleInv [[buffer(5)]],
                             constant float &deltaTime [[buffer(6)]],
                             device float4 *originalPosition [[buffer(7)]],
                             device float4 *velocity [[buffer(8)]],
                             constant float &r0 [[buffer(9)]],
                             const device float2 *neighborMap [[buffer(10)]],
                             constant uint &particleCount [[buffer(11)]],
                             constant int &timestepIndex [[buffer(12)]],
                             constant int &mode [[buffer(13)]],
                             uint serialId [[thread_position_in_grid]]) {
  if (serialId >= particleCount)
    return;

  const uint sortedParticleId = sortedParticleIdBySerialId[serialId];
  // sortedCellAndSerialId[sortedParticleId].y is the inverse of
  // sortedParticleIdBySerialId[serialId], so it always equals serialId.
  // The OpenCL kernel had an extra indirection because its thread index
  // was not the serial ID; in this Metal port threads are dispatched by
  // serialId directly, making the round-trip redundant.

  if (static_cast<int>(originalPosition[serialId].w) == kBoundaryParticle)
    return;

  if (timestepIndex == 0) {
    // Store this timestep's total acceleration (non-pressure + pressure)
    // into the [2N..3N) section, indexed by serialId so it survives
    // the next timestep's re-sort. The leapfrog velocity update (mode 1)
    // reads it back as accelerationPrev for:
    //   v(t+dt) = v(t) + (a(t) + a(t+dt)) * dt / 2
    acceleration[particleCount * 2 + serialId] =
        acceleration[sortedParticleId] +
        acceleration[particleCount + sortedParticleId];
    return;
  }

  float4 accelerationPrev = acceleration[particleCount * 2 + serialId];
  accelerationPrev.w = 0.0f;

  float4 accelerationCurrent = acceleration[sortedParticleId] +
                               acceleration[particleCount + sortedParticleId];
  accelerationCurrent.w = 0.0f;

  float4 velocityT = sortedVelocity[sortedParticleId];
  float particleType = originalPosition[serialId].w;

  if (mode == 0) {
    float4 positionT = sortedPosition[sortedParticleId];
    float4 positionNew =
        positionT + (velocityT * deltaTime +
                     accelerationPrev * deltaTime * deltaTime * 0.5f) *
                        simulationScaleInv;
    sortedPosition[sortedParticleId] = positionNew;
    sortedPosition[sortedParticleId].w = particleType;
  } else if (mode == 1) {
    float4 positionNew = sortedPosition[sortedParticleId];

    float4 velocityNew =
        velocityT + (accelerationPrev + accelerationCurrent) * deltaTime * 0.5f;

    computeInteractionWithBoundaryParticles(
        static_cast<int>(sortedParticleId), r0, neighborMap,
        sortedParticleIdBySerialId, sortedCellAndSerialId, originalPosition,
        velocity, &positionNew, true, &velocityNew, particleCount);

    velocity[serialId] = velocityNew;
    originalPosition[serialId] = positionNew;
    originalPosition[serialId].w = particleType;
    acceleration[particleCount * 2 + serialId] = accelerationCurrent;
  } else if (mode == 2) {
    float4 positionT = sortedPosition[sortedParticleId];
    float4 velocityNew = velocityT + accelerationCurrent * deltaTime;
    float4 positionNew =
        positionT + velocityNew * deltaTime * simulationScaleInv;

    computeInteractionWithBoundaryParticles(
        static_cast<int>(sortedParticleId), r0, neighborMap,
        sortedParticleIdBySerialId, sortedCellAndSerialId, originalPosition,
        velocity, &positionNew, true, &velocityNew, particleCount);

    velocity[serialId] = velocityNew;
    originalPosition[serialId] = positionNew;
    originalPosition[serialId].w = particleType;
    acceleration[particleCount * 2 + serialId] = accelerationCurrent;
  }
}

// ── Membrane interaction helpers
// ──────────────────────────────────────────────

// 3×3 determinant of column vectors (c1, c2, c3), using only .xyz components.
// Used by calculateProjectionOfPointToPlane (Cramer's rule).
inline float calcDeterminant3x3(float4 c1, float4 c2, float4 c3) {
  return c1.x * c2.y * c3.z + c1.y * c2.z * c3.x + c1.z * c2.x * c3.y -
         c1.z * c2.y * c3.x - c1.x * c2.z * c3.y - c1.y * c2.x * c3.z;
}

// Projects point ps onto the plane defined by triangle (pa, pb, pc) using
// Cramer's rule. Returns the projected point in .xyz; .w == -1.0f signals
// a degenerate triangle (zero-area, so the plane is undefined).
//
// Reference: http://ateist.spb.ru/mw/distpoint.htm
inline float4 calculateProjectionOfPointToPlane(float4 ps, float4 pa, float4 pb,
                                                float4 pc) {
  float4 pm = float4(0.0f);

  const float b_1 =
      pa.x * ((pb.y - pa.y) * (pc.z - pa.z) - (pb.z - pa.z) * (pc.y - pa.y)) +
      pa.y * ((pb.z - pa.z) * (pc.x - pa.x) - (pb.x - pa.x) * (pc.z - pa.z)) +
      pa.z * ((pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x));
  const float b_2 =
      ps.x * (pb.x - pa.x) + ps.y * (pb.y - pa.y) + ps.z * (pb.z - pa.z);
  const float b_3 =
      ps.x * (pc.x - pa.x) + ps.y * (pc.y - pa.y) + ps.z * (pc.z - pa.z);

  const float4 a_1 =
      float4((pb.y - pa.y) * (pc.z - pa.z) - (pb.z - pa.z) * (pc.y - pa.y),
             pb.x - pa.x, pc.x - pa.x, 0.0f);
  const float4 a_2 =
      float4((pb.z - pa.z) * (pc.x - pa.x) - (pb.x - pa.x) * (pc.z - pa.z),
             pb.y - pa.y, pc.y - pa.y, 0.0f);
  const float4 a_3 =
      float4((pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x),
             pb.z - pa.z, pc.z - pa.z, 0.0f);
  const float4 b = float4(b_1, b_2, b_3, 0.0f);

  const float denominator = calcDeterminant3x3(a_1, a_2, a_3);
  if (denominator != 0.0f) {
    pm.x = calcDeterminant3x3(b, a_2, a_3) / denominator;
    pm.y = calcDeterminant3x3(a_1, b, a_3) / denominator;
    pm.z = calcDeterminant3x3(a_1, a_2, b) / denominator;
  } else {
    pm.w = -1.0f; // degenerate triangle
  }
  return pm;
}

// ── Membrane interaction kernels ─────────────────────────────────────────────

// Handles particle-membrane collision detection and position correction.
// For each liquid particle, finds neighboring elastic particles that belong
// to membrane triangles, projects the particle onto each membrane plane to
// compute a surface normal, then accumulates a position correction (Ihmsen
// et al. 2010, Sec. 3.2) into the delta buffer at position[N + serialId].
//
// Thread count: particleCount (all particles). Boundary and non-liquid
// particles exit early.
//
// OpenCL signature (sphFluid.cl, lines 1120-1294):
//   __kernel void computeInteractionWithMembranes(
//       __global float4 *position,           // arg 0  (read + write delta)
//       __global float4 *velocity,           // arg 1
//       __global float4 *sortedPosition,     // arg 2  (unused)
//       __global uint2  *particleIndex,      // arg 3
//       __global uint   *particleIndexBack,  // arg 4
//       __global float2 *neighborMap,        // arg 5
//       __global int    *particleMembranesList, // arg 6
//       __global int    *membraneData,       // arg 7
//       int PARTICLE_COUNT,                  // arg 8
//       int numOfElasticP,                   // arg 9  (unused)
//       float r0                             // arg 10
//   )
kernel void computeInteractionWithMembranes(
    device float4 *position [[buffer(0)]],
    const device float4 *velocity [[buffer(1)]],
    const device uint2 *sortedCellAndSerialId [[buffer(2)]],
    const device uint *sortedParticleIdBySerialId [[buffer(3)]],
    const device float2 *neighborMap [[buffer(4)]],
    const device int *particleMembranesList [[buffer(5)]],
    const device int *membraneData [[buffer(6)]],
    constant uint &particleCount [[buffer(7)]],
    constant float &r0 [[buffer(8)]],
    uint serialId [[thread_position_in_grid]]) {
  if (serialId >= particleCount)
    return;

  // Only liquid particles interact with membranes.
  if (static_cast<int>(position[serialId].w) == kBoundaryParticle)
    return;
  if (static_cast<int>(position[serialId].w) != kLiquidParticle)
    return;

  const uint sortedParticleId = sortedParticleIdBySerialId[serialId];
  const int neighborMapOffset =
      static_cast<int>(sortedParticleId) * kMaxNeighborCount;

  // Per-neighbor membrane normal accumulator and distance.
  float4 membraneNeighborNormal[kMaxNeighborCount];
  float membraneNeighborDistance[kMaxNeighborCount];
  int membraneNeighborSerialId[kMaxNeighborCount];
  int membraneNeighborCount = 0;

  for (int i = 0; i < kMaxNeighborCount; ++i) {
    membraneNeighborNormal[i] = float4(0.0f);
  }

  // Pass 1: For each neighbor, check if it's an elastic (membrane) particle.
  // If so, project the current particle onto each membrane triangle containing
  // that neighbor and accumulate the surface normal.
  for (int neighborSlot = 0; neighborSlot < kMaxNeighborCount; ++neighborSlot) {
    const int neighborSortedParticleId =
        static_cast<int>(neighborMap[neighborMapOffset + neighborSlot].x);
    if (neighborSortedParticleId == kNoParticleId)
      break;

    const uint neighborSerialId =
        sortedCellAndSerialId[neighborSortedParticleId].y;

    if (static_cast<int>(position[neighborSerialId].w) != kElasticParticle)
      continue;

    int membraneTriangleCount = 0;

    // Distance from current particle to this elastic neighbor.
    // The OpenCL kernel zeros .z ("mv change from subscripting") — this is
    // preserved for exact behavioral compatibility even though it projects
    // the distance into the XY plane.
    float4 displacement = position[serialId] - position[neighborSerialId];
    displacement.z = 0.0f;
    const float distanceToNeighbor = sqrt(dot(displacement, displacement));

    // Check all membrane triangles that include this elastic neighbor.
    for (int membraneSlot = 0;
         membraneSlot < kMaxMembranesIncludingSameParticle; ++membraneSlot) {
      const int membraneIndex =
          particleMembranesList[static_cast<int>(neighborSerialId) *
                                    kMaxMembranesIncludingSameParticle +
                                membraneSlot];
      if (membraneIndex < 0)
        break;

      const int vi = membraneData[membraneIndex * 3 + 0];
      const int vj = membraneData[membraneIndex * 3 + 1];
      const int vk = membraneData[membraneIndex * 3 + 2];

      const float4 projection = calculateProjectionOfPointToPlane(
          position[serialId], position[vi], position[vj], position[vk]);
      if (projection.w == -1.0f)
        return; // degenerate triangle — bail out

      // Normal = particle position − its projection onto the membrane plane.
      float4 normalToPlane = position[serialId] - projection;
      const float normalLength = sqrt(normalToPlane.x * normalToPlane.x +
                                      normalToPlane.y * normalToPlane.y +
                                      normalToPlane.z * normalToPlane.z);

      if (normalLength <= 0.0f)
        return;

      normalToPlane /= normalLength;
      membraneNeighborNormal[membraneNeighborCount] += normalToPlane;
      membraneTriangleCount++;
    }

    if (membraneTriangleCount > 0) {
      membraneNeighborNormal[membraneNeighborCount] /=
          static_cast<float>(membraneTriangleCount);
      membraneNeighborDistance[membraneNeighborCount] = distanceToNeighbor;
      membraneNeighborSerialId[membraneNeighborCount] =
          static_cast<int>(neighborSerialId);
      membraneNeighborCount++;
    }
  }

  // Pass 2: Ihmsen boundary correction using the collected membrane normals.
  if (membraneNeighborCount > 0) {
    float4 n_c_i = float4(0.0f);
    float w_c_im_sum = 0.0f;
    float w_c_im_second_sum = 0.0f;

    for (int mi = 0; mi < membraneNeighborCount; ++mi) {
      const float x_im_dist = membraneNeighborDistance[mi];
      const float w_c_im = max(0.0f, (r0 - x_im_dist) / r0); // Ihmsen (10)
      const float4 n_m = membraneNeighborNormal[mi];
      n_c_i += n_m * w_c_im;                          // Ihmsen (9)
      w_c_im_sum += w_c_im;                           // Ihmsen (11) sum #1
      w_c_im_second_sum += w_c_im * (r0 - x_im_dist); // Ihmsen (11) sum #2
    }

    n_c_i.w = 0.0f;
    float n_c_i_length = dot(n_c_i.xyz, n_c_i.xyz);
    if (n_c_i_length != 0.0f) {
      n_c_i_length = sqrt(n_c_i_length);
      const float4 deltaPos =
          ((n_c_i / n_c_i_length) * w_c_im_second_sum) / w_c_im_sum;
      // Accumulate position correction into delta buffer [N..2N).
      position[particleCount + serialId].xyz += deltaPos.xyz; // Ihmsen (11)
    }
  }
}

// Applies accumulated membrane interaction deltas to particle positions.
// After computeInteractionWithMembranes has accumulated corrections into
// position[N..2N), this kernel adds them to the actual positions in [0..N).
//
// Thread count: particleCount (all particles).
// Boundary particles are skipped.
//
// OpenCL signature (sphFluid.cl, lines 1295-1328):
//   __kernel void computeInteractionWithMembranes_finalize(
//       __global float4 *position,           // arg 0
//       __global float4 *velocity,           // arg 1  (unused)
//       __global uint2  *particleIndex,      // arg 2
//       __global uint   *particleIndexBack,  // arg 3
//       int PARTICLE_COUNT                   // arg 4
//   )
kernel void computeInteractionWithMembranes_finalize(
    device float4 *position [[buffer(0)]],
    const device uint *sortedParticleIdBySerialId [[buffer(1)]],
    constant uint &particleCount [[buffer(2)]],
    uint serialId [[thread_position_in_grid]]) {
  if (serialId >= particleCount)
    return;

  if (static_cast<int>(position[serialId].w) == kBoundaryParticle)
    return;

  position[serialId] += position[particleCount + serialId];
}
