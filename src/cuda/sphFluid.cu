// CUDA kernel scaffolding — HISTORICAL placeholder.
//
// Preserved as part of the original OpenCL-mirroring work plan
// (see ./README.md). The as-built native CUDA differentiable substrate
// took a different architectural approach (XPBD constraints instead of
// PCISPH; standalone sib_cuda binary instead of an owSolver-derived
// backend) and lives in ../cuda_diff/. This file is NOT part of any
// build target and contains no executable kernel bodies; it is kept
// for context only.
//
// The original plan was to literally port each PCISPH kernel in
// src/sphFluid.cl. That plan was superseded; see ./README.md "Status"
// header for the reasoning.

#ifndef SIBERNETIC_CUDA_SPHFLUID_CU
#define SIBERNETIC_CUDA_SPHFLUID_CU

#include <cuda_runtime.h>

// -----------------------------------------------------------------------
// Neighbor search
// -----------------------------------------------------------------------

__global__ void clearBuffers(/* TODO: signature mirroring sphFluid.cl */) {
  // TODO: port from src/sphFluid.cl ::clearBuffers
}

__global__ void hashParticles(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::hashParticles
  // OpenCL: global_id → particle index, compute cell, write to gridCellIndex
}

// Sort: in CUDA we'd use thrust::sort_by_key on (cellId, particleIndex).
// No __global__ kernel needed; called from the host in CudaBackend.cpp.

__global__ void sortPostPass(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::sortPostPass
}

__global__ void indexx(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::indexx
}

__global__ void indexPostPass(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::indexPostPass
}

__global__ void findNeighbors(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::findNeighbors
  // 27-cell neighborhood walk; build per-particle neighbor list.
}

// -----------------------------------------------------------------------
// PCISPH physics
// -----------------------------------------------------------------------

__global__ void pcisphComputeDensity(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::pcisph_computeDensity
  // Wpoly6 kernel, density estimation per particle.
}

__global__ void pcisphComputeForcesAndInitPressure(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::pcisph_computeForcesAndInitPressure
  // Viscosity, surface tension, gravity, body forces.
}

__global__ void pcisphComputeElasticForces(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::pcisph_computeElasticForces
  // Spring forces between elastic-bonded particles.
  // CRITICAL: this is the kernel where Taichi's coordinate-scale bug
  // manifests. CUDA port should match OpenCL's coordinate handling
  // exactly — keep elastic forces in WORLD coordinates, not scaled.
}

__global__ void pcisphPredictPositions(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::pcisph_predictPositions
  // PCISPH iteration step 1: predict positions under current acceleration.
}

__global__ void pcisphPredictDensity(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::pcisph_predictDensity
  // PCISPH iteration step 2: re-evaluate density at predicted positions.
}

__global__ void pcisphCorrectPressure(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::pcisph_correctPressure
  // PCISPH iteration step 3: update pressure field to enforce
  // incompressibility (density deviation < 1%).
}

__global__ void pcisphComputePressureForceAcceleration(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::pcisph_computePressureForceAcceleration
  // gradWspiky; symmetric pressure force.
}

__global__ void pcisphIntegrate(/* TODO */) {
  // TODO: port from src/sphFluid.cl ::pcisph_integrate
  // Leapfrog time integration; mode 0 = position update,
  // mode 1 = velocity update.
}

// -----------------------------------------------------------------------
// Membrane interaction
// -----------------------------------------------------------------------

__global__ void clearMembraneBuffers(/* TODO */) {
  // TODO
}

__global__ void computeInteractionWithMembranes(/* TODO */) {
  // TODO
}

#endif // SIBERNETIC_CUDA_SPHFLUID_CU
