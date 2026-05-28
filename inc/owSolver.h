#pragma once

#include "owConfigProperty.h"

/// Abstract solver interface.
///
/// Both owOpenCLSolver and owMetalSolver implement this interface so that
/// owPhysicsFluidSimulator can dispatch to either backend.
class owSolver {
public:
  virtual ~owSolver() = default;

  // ── Neighbor search ──
  virtual unsigned int _runClearBuffers(owConfigProperty *config) = 0;
  virtual unsigned int _runHashParticles(owConfigProperty *config) = 0;
  virtual void _runSort(owConfigProperty *config) = 0;
  virtual unsigned int _runSortPostPass(owConfigProperty *config) = 0;
  virtual unsigned int _runIndexx(owConfigProperty *config) = 0;
  virtual void _runIndexPostPass(owConfigProperty *config) = 0;
  virtual unsigned int _runFindNeighbors(owConfigProperty *config) = 0;

  // ── PCISPH physics ──
  virtual unsigned int _run_pcisph_computeDensity(owConfigProperty *config) = 0;
  virtual unsigned int
  _run_pcisph_computeForcesAndInitPressure(owConfigProperty *config) = 0;
  virtual unsigned int
  _run_pcisph_computeElasticForces(owConfigProperty *config) = 0;
  virtual unsigned int
  _run_pcisph_predictPositions(owConfigProperty *config) = 0;
  virtual unsigned int
  _run_pcisph_predictDensity(owConfigProperty *config) = 0;
  virtual unsigned int
  _run_pcisph_correctPressure(owConfigProperty *config) = 0;
  virtual unsigned int
  _run_pcisph_computePressureForceAcceleration(owConfigProperty *config) = 0;
  virtual unsigned int _run_pcisph_integrate(int iterationCount,
                                             int pcisph_integrate_mode,
                                             owConfigProperty *config) = 0;

  // ── Membrane interaction ──
  virtual unsigned int
  _run_clearMembraneBuffers(owConfigProperty *config) = 0;
  virtual unsigned int
  _run_computeInteractionWithMembranes(owConfigProperty *config) = 0;
  virtual unsigned int
  _run_computeInteractionWithMembranes_finalize(owConfigProperty *config) = 0;

  // ── Buffer readback ──
  virtual void read_position_buffer(float *position_cpp,
                                    owConfigProperty *config) = 0;
  virtual void read_velocity_buffer(float *velocity_cpp,
                                    owConfigProperty *config) = 0;
  virtual void read_density_buffer(float *density_cpp,
                                   owConfigProperty *config) = 0;
  virtual void read_particleIndex_buffer(unsigned int *particleIndex_cpp,
                                         owConfigProperty *config) = 0;
  virtual void read_pressure_buffer(float *pressure_cpp,
                                    owConfigProperty *config) = 0;

  // ── Misc ──
  virtual void updateMuscleActivityData(float *muscle_activation_signal_cpp,
                                        owConfigProperty *config) = 0;
  virtual void reset(const float *position_cpp, const float *velocity_cpp,
                     owConfigProperty *config,
                     const float *elasticConnectionsData_cpp = nullptr,
                     const int *membraneData_cpp = nullptr,
                     const int *particleMembranesList_cpp = nullptr) = 0;
};
