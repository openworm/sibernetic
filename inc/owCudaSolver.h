// owCudaSolver — HISTORICAL placeholder for the originally-planned
// integration of a native CUDA backend into the owSolver virtual
// interface. NOT part of any current build target.
//
// The as-built native CUDA differentiable substrate took a different
// architectural approach and lives in src/cuda_diff/ as a standalone
// binary (sib_cuda), parallel to src/metal_diff/sib_metal. Integrating
// either substrate into the main Release/Sibernetic binary via an
// owSolver-derived class is a natural follow-up — the top-level README
// "Next steps for closing the parity gap" #5 proposes this for Metal;
// the CUDA analog would slot in the same way. Neither has been done;
// this header is preserved as scaffolding for that future work.
//
// See ../src/cuda/README.md for the original work plan and
// ../src/cuda_diff/README.md for the as-built substrate.

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
