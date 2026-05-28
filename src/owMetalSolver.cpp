#include "owMetalSolver.h"

#ifdef SIBERNETIC_USE_METAL

#include "owOpenCLConstant.h"

#include "kernels/ClearBuffersKernel.h"
#include "kernels/ClearMembraneBuffersKernel.h"
#include "kernels/ComputeDensityKernel.h"
#include "kernels/ComputeInteractionWithMembranesKernel.h"
#include "kernels/FindNeighborsKernel.h"
#include "kernels/HashParticlesKernel.h"
#include "kernels/IndexxKernel.h"
#include "kernels/PcisphComputeElasticForcesKernel.h"
#include "kernels/PcisphComputeForcesKernel.h"
#include "kernels/PcisphComputePressureForceAccelerationKernel.h"
#include "kernels/PcisphCorrectPressureKernel.h"
#include "kernels/PcisphIntegrateKernel.h"
#include "kernels/PcisphPredictDensityKernel.h"
#include "kernels/PcisphPredictPositionsKernel.h"
#include "kernels/SortPostPassKernel.h"

#include <cstring>

using namespace Sibernetic;

// ── qsort comparator (same as owOpenCLSolver) ──

static int comparator(const void *v1, const void *v2) {
  const int *f1 = static_cast<const int *>(v1);
  const int *f2 = static_cast<const int *>(v2);
  if (f1[0] < f2[0])
    return -1;
  if (f1[0] > f2[0])
    return +1;
  return 0;
}

// ── Helpers ──

static MTL::Buffer *newSharedBuffer(MTL::Device *device, size_t bytes) {
  return device->newBuffer(bytes, MTL::ResourceStorageModeShared);
}

static MTL::Buffer *newSharedBuffer(MTL::Device *device, const void *data,
                                    size_t bytes) {
  return device->newBuffer(data, bytes, MTL::ResourceStorageModeShared);
}

// ── Constructor / Destructor ──

owMetalSolver::owMetalSolver(const float *position_cpp,
                             const float *velocity_cpp,
                             owConfigProperty *config,
                             const float *elasticConnectionsData_cpp,
                             const int *membraneData_cpp,
                             const int *particleMembranesList_cpp) {
  backend_ = std::make_unique<MetalBackend>();
  initializeBuffers(position_cpp, velocity_cpp, config,
                    elasticConnectionsData_cpp, membraneData_cpp,
                    particleMembranesList_cpp);
}

owMetalSolver::~owMetalSolver() { releaseBuffers(); }

// ── Buffer management ──

void owMetalSolver::initializeBuffers(
    const float *position_cpp, const float *velocity_cpp,
    owConfigProperty *config, const float *elasticConnectionsData_cpp,
    const int *membraneData_cpp, const int *particleMembranesList_cpp) {

  MTL::Device *dev = backend_->device();
  const int N = config->getParticleCount();

  accelerationBuf_ = newSharedBuffer(dev, N * sizeof(float) * 4 * 3);
  memset(accelerationBuf_->contents(), 0, accelerationBuf_->length());

  gridCellIndexBuf_ =
      newSharedBuffer(dev, (config->gridCellCount + 1) * sizeof(unsigned int));
  gridCellIndexFixedUpBuf_ =
      newSharedBuffer(dev, (config->gridCellCount + 1) * sizeof(unsigned int));

  neighborMapBuf_ =
      newSharedBuffer(dev, MAX_NEIGHBOR_COUNT * N * sizeof(float) * 2);
  particleIndexBuf_ = newSharedBuffer(dev, N * sizeof(unsigned int) * 2);
  particleIndexBackBuf_ = newSharedBuffer(dev, N * sizeof(unsigned int));

  positionBuf_ = newSharedBuffer(dev, N * sizeof(float) * 4 * 2);
  memset(positionBuf_->contents(), 0, positionBuf_->length());
  memcpy(positionBuf_->contents(), position_cpp, N * sizeof(float) * 4);

  pressureBuf_ = newSharedBuffer(dev, N * sizeof(float));
  rhoBuf_ = newSharedBuffer(dev, N * sizeof(float) * 2);
  sortedPositionBuf_ = newSharedBuffer(dev, N * sizeof(float) * 4 * 2);
  sortedVelocityBuf_ = newSharedBuffer(dev, N * sizeof(float) * 4);

  velocityBuf_ = newSharedBuffer(dev, N * sizeof(float) * 4 * 2);
  memset(velocityBuf_->contents(), 0, velocityBuf_->length());
  memcpy(velocityBuf_->contents(), velocity_cpp, N * sizeof(float) * 4);

  muscleActivationSignalBuf_ =
      newSharedBuffer(dev, config->MUSCLE_COUNT * sizeof(float));
  memset(muscleActivationSignalBuf_->contents(), 0,
         muscleActivationSignalBuf_->length());

  if (membraneData_cpp && particleMembranesList_cpp) {
    membraneDataBuf_ = newSharedBuffer(dev, membraneData_cpp,
                                       config->numOfMembranes * 3 * sizeof(int));
    particleMembranesListBuf_ = newSharedBuffer(
        dev, particleMembranesList_cpp,
        config->numOfElasticP * MAX_MEMBRANES_INCLUDING_SAME_PARTICLE *
            sizeof(int));
  }

  if (elasticConnectionsData_cpp) {
    elasticConnectionsDataBuf_ = newSharedBuffer(
        dev, elasticConnectionsData_cpp,
        config->numOfElasticP * MAX_NEIGHBOR_COUNT * sizeof(float) * 4);
  }
}

void owMetalSolver::releaseBuffers() {
  auto release = [](MTL::Buffer *&buf) {
    if (buf) {
      buf->release();
      buf = nullptr;
    }
  };
  release(accelerationBuf_);
  release(gridCellIndexBuf_);
  release(gridCellIndexFixedUpBuf_);
  release(neighborMapBuf_);
  release(particleIndexBuf_);
  release(particleIndexBackBuf_);
  release(positionBuf_);
  release(pressureBuf_);
  release(rhoBuf_);
  release(sortedPositionBuf_);
  release(sortedVelocityBuf_);
  release(velocityBuf_);
  release(muscleActivationSignalBuf_);
  release(membraneDataBuf_);
  release(particleMembranesListBuf_);
  release(elasticConnectionsDataBuf_);
}

void owMetalSolver::reset(const float *position_cpp,
                          const float *velocity_cpp,
                          owConfigProperty *config,
                          const float *elasticConnectionsData_cpp,
                          const int *membraneData_cpp,
                          const int *particleMembranesList_cpp) {
  releaseBuffers();
  initializeBuffers(position_cpp, velocity_cpp, config,
                    elasticConnectionsData_cpp, membraneData_cpp,
                    particleMembranesList_cpp);
}

// ── Kernel dispatches ──

unsigned int owMetalSolver::_runClearBuffers(owConfigProperty *config) {
  ClearBuffersMetalArgs args{};
  args.neighborMap = neighborMapBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(kClearBuffersKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

unsigned int owMetalSolver::_runHashParticles(owConfigProperty *config) {
  HashParticlesMetalArgs args{};
  args.position = positionBuf_;
  args.gridCellsX = config->gridCellsX;
  args.gridCellsY = config->gridCellsY;
  args.gridCellsZ = config->gridCellsZ;
  args.hashGridCellSizeInv = config->getConst("hashGridCellSizeInv");
  args.xmin = config->xmin;
  args.ymin = config->ymin;
  args.zmin = config->zmin;
  args.particleIndex = particleIndexBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(kHashParticlesKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

void owMetalSolver::_runSort(owConfigProperty *config) {
  // Must wait for hashParticles to finish before CPU reads particleIndexBuf_.
  backend_->finish();
  // CPU-side qsort on shared Metal buffer (same approach as OpenCL solver).
  int *data = static_cast<int *>(particleIndexBuf_->contents());
  qsort(data, config->getParticleCount(), 2 * sizeof(int), comparator);
}

unsigned int owMetalSolver::_runSortPostPass(owConfigProperty *config) {
  SortPostPassMetalArgs args{};
  args.particleIndex = particleIndexBuf_;
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.position = positionBuf_;
  args.velocity = velocityBuf_;
  args.sortedPosition = sortedPositionBuf_;
  args.sortedVelocity = sortedVelocityBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(kSortPostPassKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

unsigned int owMetalSolver::_runIndexx(owConfigProperty *config) {
  IndexxMetalArgs args{};
  args.particleIndex = particleIndexBuf_;
  args.gridCellCount = config->gridCellCount;
  args.gridCellIndex = gridCellIndexBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  // indexx dispatches over grid cells, not particles.
  backend_->dispatch(kIndexxKernelName, config->gridCellCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

void owMetalSolver::_runIndexPostPass(owConfigProperty *config) {
  // Must wait for indexx to finish before CPU reads gridCellIndexBuf_.
  backend_->finish();
  // CPU-side fixup on shared Metal buffers (same approach as OpenCL solver).
  int *src = static_cast<int *>(gridCellIndexBuf_->contents());
  int *dst = static_cast<int *>(gridCellIndexFixedUpBuf_->contents());
  const int cellCount = static_cast<int>(config->gridCellCount);

  memcpy(dst, src, (cellCount + 1) * sizeof(int));

  int recentNonEmptyCell = cellCount;
  for (int i = cellCount; i >= 0; i--) {
    if (dst[i] == NO_CELL_ID)
      dst[i] = recentNonEmptyCell;
    else
      recentNonEmptyCell = dst[i];
  }
}

unsigned int owMetalSolver::_runFindNeighbors(owConfigProperty *config) {
  FindNeighborsMetalArgs args{};
  args.gridCellIndexFixedUp = gridCellIndexFixedUpBuf_;
  args.sortedPosition = sortedPositionBuf_;
  args.gridCellCount = config->gridCellCount;
  args.gridCellsX = config->gridCellsX;
  args.gridCellsY = config->gridCellsY;
  args.gridCellsZ = config->gridCellsZ;
  args.h = config->getConst("h");
  args.hashGridCellSize = config->getConst("hashGridCellSize");
  args.hashGridCellSizeInv = config->getConst("hashGridCellSizeInv");
  args.simulationScale = config->getConst("simulationScale");
  args.xmin = config->xmin;
  args.ymin = config->ymin;
  args.zmin = config->zmin;
  args.neighborMap = neighborMapBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(kFindNeighborsKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

// ── PCISPH kernels ──

unsigned int
owMetalSolver::_run_pcisph_computeDensity(owConfigProperty *config) {
  ComputeDensityMetalArgs args{};
  args.neighborMap = neighborMapBuf_;
  args.massMultWpoly6Coefficient =
      config->getConst("mass_mult_Wpoly6Coefficient");
  args.hScaled2 = config->getConst("_hScaled2");
  args.rho = rhoBuf_;
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(kComputeDensityKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

unsigned int owMetalSolver::_run_pcisph_computeForcesAndInitPressure(
    owConfigProperty *config) {
  PcisphComputeForcesMetalArgs args{};
  args.neighborMap = neighborMapBuf_;
  args.rho = rhoBuf_;
  args.pressure = pressureBuf_;
  args.sortedPosition = sortedPositionBuf_;
  args.sortedVelocity = sortedVelocityBuf_;
  args.acceleration = accelerationBuf_;
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.surfTensCoeff = config->getConst("surfTensCoeff");
  args.massMultLaplacianWviscosityCoeff =
      config->getConst("mass_mult_divgradWviscosityCoefficient");
  args.hScaled = config->getConst("_hScaled");
  args.mu = config->getConst("viscosity");
  args.gravity_x = config->getConst("gravity_x");
  args.gravity_y = config->getConst("gravity_y");
  args.gravity_z = config->getConst("gravity_z");
  args.position = positionBuf_;
  args.particleIndex = particleIndexBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());
  args.mass = config->getConst("mass");

  backend_->dispatch(kPcisphComputeForcesKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

unsigned int
owMetalSolver::_run_pcisph_computeElasticForces(owConfigProperty *config) {
  if (config->numOfElasticP == 0)
    return 0;

  PcisphComputeElasticForcesMetalArgs args{};
  args.sortedPosition = sortedPositionBuf_;
  args.acceleration = accelerationBuf_;
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.sortedCellAndSerialId = particleIndexBuf_;
  args.maxMuscleForce = config->getConst("max_muscle_force");
  args.simulationScale = config->getConst("simulationScale");
  args.numOfElasticP = config->numOfElasticP;
  args.elasticConnectionsData = elasticConnectionsDataBuf_;
  args.muscleCount = config->MUSCLE_COUNT;
  args.muscleActivationSignal = muscleActivationSignalBuf_;
  args.originalPosition = positionBuf_;
  args.elasticityCoefficient = config->getConst("elasticityCoefficient");

  backend_->dispatch(kPcisphComputeElasticForcesKernelName,
                     config->numOfElasticP,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

unsigned int
owMetalSolver::_run_pcisph_predictPositions(owConfigProperty *config) {
  PcisphPredictPositionsMetalArgs args{};
  args.acceleration = accelerationBuf_;
  args.sortedPosition = sortedPositionBuf_;
  args.sortedVelocity = sortedVelocityBuf_;
  args.sortedCellAndSerialId = particleIndexBuf_;
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.gravitationalAccelerationX = config->getConst("gravity_x");
  args.gravitationalAccelerationY = config->getConst("gravity_y");
  args.gravitationalAccelerationZ = config->getConst("gravity_z");
  args.simulationScaleInv = config->getConst("simulationScaleInv");
  args.deltaTime = config->getTimeStep();
  args.originalPosition = positionBuf_;
  args.velocity = velocityBuf_;
  args.r0 = config->getConst("r0");
  args.neighborMap = neighborMapBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(kPcisphPredictPositionsKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

unsigned int
owMetalSolver::_run_pcisph_predictDensity(owConfigProperty *config) {
  PcisphPredictDensityMetalArgs args{};
  args.neighborMap = neighborMapBuf_;
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.massMultWpoly6Coefficient =
      config->getConst("mass_mult_Wpoly6Coefficient");
  args.h = config->getConst("h");
  args.restDensity = config->getConst("rho0");
  args.simulationScale = config->getConst("simulationScale");
  args.sortedPosition = sortedPositionBuf_;
  args.rho = rhoBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(kPcisphPredictDensityKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

unsigned int
owMetalSolver::_run_pcisph_correctPressure(owConfigProperty *config) {
  PcisphCorrectPressureMetalArgs args{};
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.restDensity = config->getConst("rho0");
  args.pressure = pressureBuf_;
  args.rho = rhoBuf_;
  args.delta = config->getDelta();
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(kPcisphCorrectPressureKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

unsigned int owMetalSolver::_run_pcisph_computePressureForceAcceleration(
    owConfigProperty *config) {
  PcisphComputePressureForceAccelerationMetalArgs args{};
  args.neighborMap = neighborMapBuf_;
  args.pressure = pressureBuf_;
  args.rho = rhoBuf_;
  args.sortedPosition = sortedPositionBuf_;
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.delta = config->getDelta();
  args.massMultGradWspikyCoefficient =
      config->getConst("mass_mult_gradWspikyCoefficient");
  args.h = config->getConst("h");
  args.simulationScale = config->getConst("simulationScale");
  args.restDensity = config->getConst("rho0");
  args.acceleration = accelerationBuf_;
  args.originalPosition = positionBuf_;
  args.sortedCellAndSerialId = particleIndexBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(
      kPcisphComputePressureForceAccelerationKernelName, args.particleCount,
      [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
      /*waitForCompletion=*/false);
  return 0;
}

unsigned int owMetalSolver::_run_pcisph_integrate(int iterationCount,
                                                  int pcisph_integrate_mode,
                                                  owConfigProperty *config) {
  PcisphIntegrateMetalArgs args{};
  args.acceleration = accelerationBuf_;
  args.sortedPosition = sortedPositionBuf_;
  args.sortedVelocity = sortedVelocityBuf_;
  args.sortedCellAndSerialId = particleIndexBuf_;
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.simulationScaleInv = config->getConst("simulationScaleInv");
  args.deltaTime = config->getTimeStep();
  args.originalPosition = positionBuf_;
  args.velocity = velocityBuf_;
  args.r0 = config->getConst("r0");
  args.neighborMap = neighborMapBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());
  args.timestepIndex = iterationCount;
  args.mode = pcisph_integrate_mode;

  backend_->dispatch(kPcisphIntegrateKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

// ── Membrane kernels ──

unsigned int
owMetalSolver::_run_clearMembraneBuffers(owConfigProperty *config) {
  ClearMembraneBuffersMetalArgs args{};
  args.position = positionBuf_;
  args.velocity = velocityBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(kClearMembraneBuffersKernelName, args.particleCount,
                     [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
                     /*waitForCompletion=*/false);
  return 0;
}

unsigned int owMetalSolver::_run_computeInteractionWithMembranes(
    owConfigProperty *config) {
  ComputeInteractionWithMembranesMetalArgs args{};
  args.position = positionBuf_;
  args.velocity = velocityBuf_;
  args.sortedCellAndSerialId = particleIndexBuf_;
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.neighborMap = neighborMapBuf_;
  args.particleMembranesList = particleMembranesListBuf_;
  args.membraneData = membraneDataBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());
  args.r0 = config->getConst("r0");

  backend_->dispatch(
      kComputeInteractionWithMembranesKernelName, args.particleCount,
      [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
      /*waitForCompletion=*/false);
  return 0;
}

unsigned int owMetalSolver::_run_computeInteractionWithMembranes_finalize(
    owConfigProperty *config) {
  ComputeInteractionWithMembranesFinalizeMetalArgs args{};
  args.position = positionBuf_;
  args.sortedParticleIdBySerialId = particleIndexBackBuf_;
  args.particleCount = static_cast<uint32_t>(config->getParticleCount());

  backend_->dispatch(
      kComputeInteractionWithMembranesFinalizeKernelName, args.particleCount,
      [&](MTL::ComputeCommandEncoder *enc) { args.bind(enc); },
      /*waitForCompletion=*/false);
  return 0;
}

// ── Buffer readback (shared memory – just memcpy) ──

void owMetalSolver::read_position_buffer(float *position_cpp,
                                         owConfigProperty *config) {
  backend_->finish();
  memcpy(position_cpp, positionBuf_->contents(),
         config->getParticleCount() * sizeof(float) * 4);
}

void owMetalSolver::read_velocity_buffer(float *velocity_cpp,
                                         owConfigProperty *config) {
  backend_->finish();
  memcpy(velocity_cpp, velocityBuf_->contents(),
         config->getParticleCount() * sizeof(float) * 4);
}

void owMetalSolver::read_density_buffer(float *density_cpp,
                                        owConfigProperty *config) {
  backend_->finish();
  memcpy(density_cpp, rhoBuf_->contents(),
         config->getParticleCount() * sizeof(float));
}

void owMetalSolver::read_particleIndex_buffer(unsigned int *particleIndex_cpp,
                                              owConfigProperty *config) {
  backend_->finish();
  memcpy(particleIndex_cpp, particleIndexBuf_->contents(),
         config->getParticleCount() * sizeof(unsigned int) * 2);
}

void owMetalSolver::read_pressure_buffer(float *pressure_cpp,
                                         owConfigProperty *config) {
  backend_->finish();
  memcpy(pressure_cpp, pressureBuf_->contents(),
         config->getParticleCount() * sizeof(float));
}

void owMetalSolver::updateMuscleActivityData(
    float *muscle_activation_signal_cpp, owConfigProperty *config) {
  memcpy(muscleActivationSignalBuf_->contents(), muscle_activation_signal_cpp,
         config->MUSCLE_COUNT * sizeof(float));
}

#endif // SIBERNETIC_USE_METAL
