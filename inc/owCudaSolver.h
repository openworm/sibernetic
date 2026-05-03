// owCudaSolver — native CUDA backend for Sibernetic's PCISPH solver.
//
// STATUS: skeleton only. Public interface defined; method bodies are
// declared but not implemented. See ../src/cuda/README.md for the
// implementation plan.
//
// Mirrors owOpenCLSolver.h's interface exactly so owPhysicsFluidSimulator
// can dispatch between OpenCL, Metal (PR #222), and CUDA via the same
// base class. Once PR #222's owSolver.h lands, this should inherit from
// owSolver and override the virtual methods rather than duplicating the
// signatures.

#ifndef OWCUDASOLVER_H_
#define OWCUDASOLVER_H_

#include "owConfigProperty.h"

class owCudaSolver {
 public:
  // TODO: constructors mirroring owOpenCLSolver
  // owCudaSolver(float *positions, float *velocities, owConfigProperty *config,
  //              float *connections = nullptr, float *membranes = nullptr,
  //              int *particleMembranes = nullptr);
  // ~owCudaSolver();

  // ── Neighbor search ──
  // unsigned int _runClearBuffers(owConfigProperty *config);
  // unsigned int _runHashParticles(owConfigProperty *config);
  // void _runSort(owConfigProperty *config);
  // unsigned int _runSortPostPass(owConfigProperty *config);
  // unsigned int _runIndexx(owConfigProperty *config);
  // void _runIndexPostPass(owConfigProperty *config);
  // unsigned int _runFindNeighbors(owConfigProperty *config);

  // ── PCISPH physics ──
  // unsigned int _run_pcisph_computeDensity(owConfigProperty *config);
  // unsigned int _run_pcisph_computeForcesAndInitPressure(owConfigProperty *config);
  // unsigned int _run_pcisph_computeElasticForces(owConfigProperty *config);
  // unsigned int _run_pcisph_predictPositions(owConfigProperty *config);
  // unsigned int _run_pcisph_predictDensity(owConfigProperty *config);
  // unsigned int _run_pcisph_correctPressure(owConfigProperty *config);
  // unsigned int _run_pcisph_computePressureForceAcceleration(owConfigProperty *config);
  // unsigned int _run_pcisph_integrate(int iterationCount,
  //                                     int pcisph_integrate_mode,
  //                                     owConfigProperty *config);

  // ── Membrane interaction ──
  // unsigned int _run_clearMembraneBuffers(owConfigProperty *config);
  // unsigned int _run_computeInteractionWithMembranes(owConfigProperty *config);

  // ── Buffer access (mirrors owOpenCLSolver pattern) ──
  // void read_position_buffer(float *positions, owConfigProperty *config);
  // void read_velocity_buffer(float *velocities, owConfigProperty *config);
  // void read_density_buffer(float *density, owConfigProperty *config);

 private:
  // TODO: CUDA stream, device pointers for each buffer (positions,
  // velocities, accelerations, density, pressure, neighbor map, grid
  // index, etc.). See PR #222's owMetalSolver.cpp for the equivalent
  // Metal implementation as a structural template.
};

#endif  // OWCUDASOLVER_H_
