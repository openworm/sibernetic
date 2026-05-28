#pragma once

#ifdef SIBERNETIC_USE_METAL

#include "owConfigProperty.h"
#include "owSolver.h"

#include "Metal/Metal/MTLBuffer.hpp"
#include "../src/backend/MetalBackend.h"

#include <memory>

/// Metal GPU solver – mirrors owOpenCLSolver's public interface using the
/// Metal compute backend and the kernel-argument structs from src/kernels/.
class owMetalSolver : public owSolver {
public:
  owMetalSolver(const float *position_cpp, const float *velocity_cpp,
                owConfigProperty *config,
                const float *elasticConnectionsData_cpp = nullptr,
                const int *membraneData_cpp = nullptr,
                const int *particleMembranesList_cpp = nullptr);
  ~owMetalSolver() override;

  // ── Neighbor search ──
  unsigned int _runClearBuffers(owConfigProperty *config) override;
  unsigned int _runHashParticles(owConfigProperty *config) override;
  void _runSort(owConfigProperty *config) override;
  unsigned int _runSortPostPass(owConfigProperty *config) override;
  unsigned int _runIndexx(owConfigProperty *config) override;
  void _runIndexPostPass(owConfigProperty *config) override;
  unsigned int _runFindNeighbors(owConfigProperty *config) override;

  // ── PCISPH physics ──
  unsigned int _run_pcisph_computeDensity(owConfigProperty *config) override;
  unsigned int
  _run_pcisph_computeForcesAndInitPressure(owConfigProperty *config) override;
  unsigned int
  _run_pcisph_computeElasticForces(owConfigProperty *config) override;
  unsigned int
  _run_pcisph_predictPositions(owConfigProperty *config) override;
  unsigned int _run_pcisph_predictDensity(owConfigProperty *config) override;
  unsigned int _run_pcisph_correctPressure(owConfigProperty *config) override;
  unsigned int
  _run_pcisph_computePressureForceAcceleration(
      owConfigProperty *config) override;
  unsigned int _run_pcisph_integrate(int iterationCount,
                                     int pcisph_integrate_mode,
                                     owConfigProperty *config) override;

  // ── Membrane interaction ──
  unsigned int _run_clearMembraneBuffers(owConfigProperty *config) override;
  unsigned int
  _run_computeInteractionWithMembranes(owConfigProperty *config) override;
  unsigned int _run_computeInteractionWithMembranes_finalize(
      owConfigProperty *config) override;

  // ── Buffer readback ──
  void read_position_buffer(float *position_cpp,
                            owConfigProperty *config) override;
  void read_velocity_buffer(float *velocity_cpp,
                            owConfigProperty *config) override;
  void read_density_buffer(float *density_cpp,
                           owConfigProperty *config) override;
  void read_particleIndex_buffer(unsigned int *particleIndex_cpp,
                                 owConfigProperty *config) override;
  void read_pressure_buffer(float *pressure_cpp,
                            owConfigProperty *config) override;

  // ── Misc ──
  void updateMuscleActivityData(float *muscle_activation_signal_cpp,
                                owConfigProperty *config) override;
  void reset(const float *position_cpp, const float *velocity_cpp,
             owConfigProperty *config,
             const float *elasticConnectionsData_cpp = nullptr,
             const int *membraneData_cpp = nullptr,
             const int *particleMembranesList_cpp = nullptr) override;

private:
  void initializeBuffers(const float *position_cpp,
                         const float *velocity_cpp, owConfigProperty *config,
                         const float *elasticConnectionsData_cpp,
                         const int *membraneData_cpp,
                         const int *particleMembranesList_cpp);
  void releaseBuffers();

  std::unique_ptr<Sibernetic::MetalBackend> backend_;

  // ── GPU buffers (MTL::Buffer*, shared storage mode) ──
  MTL::Buffer *accelerationBuf_ = nullptr; // 3×N float4
  MTL::Buffer *gridCellIndexBuf_ = nullptr;
  MTL::Buffer *gridCellIndexFixedUpBuf_ = nullptr;
  MTL::Buffer *neighborMapBuf_ = nullptr; // N×32 float2
  MTL::Buffer *particleIndexBuf_ = nullptr; // N uint2
  MTL::Buffer *particleIndexBackBuf_ = nullptr; // N uint
  MTL::Buffer *positionBuf_ = nullptr; // 2×N float4
  MTL::Buffer *pressureBuf_ = nullptr; // N float
  MTL::Buffer *rhoBuf_ = nullptr; // 2×N float
  MTL::Buffer *sortedPositionBuf_ = nullptr; // 2×N float4
  MTL::Buffer *sortedVelocityBuf_ = nullptr; // N float4
  MTL::Buffer *velocityBuf_ = nullptr; // 2×N float4
  MTL::Buffer *muscleActivationSignalBuf_ = nullptr;
  MTL::Buffer *membraneDataBuf_ = nullptr;
  MTL::Buffer *particleMembranesListBuf_ = nullptr;
  MTL::Buffer *elasticConnectionsDataBuf_ = nullptr;
};

#endif // SIBERNETIC_USE_METAL
