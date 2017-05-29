/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "owPhysicsFluidSimulator.h"
#include "owSignalSimulator.h"
#include "owVtkExport.h"

/** Constructor method for owPhysicsFluidSimulator.
 *
 *  @param helper
 *  pointer to owHelper object with helper function.
 *  @param dev_type
 *  defines preferable device type for current configuration
 */
owPhysicsFluidSimulator::owPhysicsFluidSimulator(owHelper *helper, int argc,
                                                 char **argv) {
  // int generateInitialConfiguration = 1;//1 to generate initial configuration,
  // 0 - load from file

  try {
    iterationCount = 0;
    config = new owConfigProperty(argc, argv);
    // LOAD FROM FILE
    owHelper::preLoadConfiguration(config);
    config->initGridCells();
    position_cpp = new float[4 * config->getParticleCount()];
    velocity_cpp = new float[4 * config->getParticleCount()];
    muscle_activation_signal_cpp = new float[config->MUSCLE_COUNT];
    if (config->numOfElasticP != 0)
      elasticConnectionsData_cpp =
          new float[4 * config->numOfElasticP * MAX_NEIGHBOR_COUNT];
    if (config->numOfMembranes <= 0)
      membraneData_cpp = NULL;
    else
      membraneData_cpp = new int[config->numOfMembranes * 3];
    if (config->numOfElasticP <= 0)
      particleMembranesList_cpp = NULL;
    else
      particleMembranesList_cpp =
          new int[config->numOfElasticP *
                  MAX_MEMBRANES_INCLUDING_SAME_PARTICLE];
    for (unsigned int i = 0; i < config->MUSCLE_COUNT; ++i) {
      muscle_activation_signal_cpp[i] = 0.f;
    }

    // The buffers listed below are only for usability and debug
    density_cpp = new float[1 * config->getParticleCount()];
    particleIndex_cpp = new unsigned int[config->getParticleCount() * 2];

    // LOAD FROM FILE
    owHelper::loadConfiguration(
        position_cpp, velocity_cpp, elasticConnectionsData_cpp,
        membraneData_cpp, particleMembranesList_cpp,
        config); // Load configuration from file to buffer
    this->helper = helper;
    if (config->numOfElasticP != 0) {
      ocl_solver = new owOpenCLSolver(
          position_cpp, velocity_cpp, config, elasticConnectionsData_cpp,
          membraneData_cpp,
          particleMembranesList_cpp); // Create new openCLsolver instance
    } else
      ocl_solver =
          new owOpenCLSolver(position_cpp, velocity_cpp,
                             config); // Create new openCLsolver instance
  } catch (std::runtime_error &ex) {
    /* Clearing all allocated buffers and created object only not ocl_solver
     * case it wont be created yet only if exception is throwing from its
     * constructor
     * but in this case ocl_solver wont be created
     * */
    destroy();
    delete config;
    throw;
  }
}
/** Reset simulation
 *
 *  Restart simulation with new or current simulation configuration.
 *  It redefines all required data buffers and restart owOpenCLSolver
 *  by run owOpenCLSolver::reset(...).
 */
void owPhysicsFluidSimulator::reset() {
  // Free all buffers
  destroy();
  config->resetNeuronSimulation();
  iterationCount = 0;
  config->numOfBoundaryP = 0;
  config->numOfElasticP = 0;
  config->numOfLiquidP = 0;
  config->numOfMembranes = 0;
  // LOAD FROM FILE
  owHelper::preLoadConfiguration(config);
  config->initGridCells();
  position_cpp = new float[4 * config->getParticleCount()];
  velocity_cpp = new float[4 * config->getParticleCount()];
  muscle_activation_signal_cpp = new float[config->MUSCLE_COUNT];
  if (config->numOfElasticP != 0)
    elasticConnectionsData_cpp =
        new float[4 * config->numOfElasticP * MAX_NEIGHBOR_COUNT];
  if (config->numOfMembranes <= 0)
    membraneData_cpp = NULL;
  else
    membraneData_cpp = new int[config->numOfMembranes * 3];
  if (config->numOfElasticP <= 0)
    particleMembranesList_cpp = NULL;
  else
    particleMembranesList_cpp =
        new int[config->numOfElasticP * MAX_MEMBRANES_INCLUDING_SAME_PARTICLE];
  for (unsigned int i = 0; i < config->MUSCLE_COUNT; ++i) {
    muscle_activation_signal_cpp[i] = 0.f;
  }
  // The buffers listed below are only for usability and debug
  density_cpp = new float[1 * config->getParticleCount()];
  particleIndex_cpp = new unsigned int[config->getParticleCount() * 2];
  // LOAD FROM FILE
  owHelper::loadConfiguration(position_cpp, velocity_cpp,
                              elasticConnectionsData_cpp, membraneData_cpp,
                              particleMembranesList_cpp,
                              config); // Load configuration from file to buffer
  if (config->numOfElasticP != 0) {
    ocl_solver->reset(
        position_cpp, velocity_cpp, config, elasticConnectionsData_cpp,
        membraneData_cpp,
        particleMembranesList_cpp); // Create new openCLsolver instance
  } else
    ocl_solver->reset(position_cpp, velocity_cpp,
                      config); // Create new openCLsolver instance
}
/** Run one simulation step
 *
 *  Run simulation step in pipeline manner.
 *  It starts with neighbor search algorithm than
 *  physic simulation algorithms: PCI SPH [1],
 *  elastic matter simulation, boundary handling [2],
 *  membranes handling and finally numerical integration.
 *  [1] http://www.ifi.uzh.ch/vmml/publications/pcisph/pcisph.pdf
 *  [2] M. Ihmsen, N. Akinci, M. Gissler, M. Teschner,
 *      Boundary Handling and Adaptive Time-stepping for PCISPH
 *      Proc. VRIPHYS, Copenhagen, Denmark, pp. 79-88, Nov 11-12, 2010
 *
 *  @param looad_to
 *  If it's true than Sibernetic works "load simulation data in file" mode.
 */
double owPhysicsFluidSimulator::simulationStep(const bool load_to) {
  int iter = 0; // PCISPH prediction-correction iterations counter
  //
  // now we will implement sensory system of the c. elegans worm, mechanosensory
  // one
  // here we plan to implement the part of openworm sensory system, which is
  // still
  // one of the grand challenges of this project

  // if(iterationCount==0) return 0.0;//uncomment this line to stop movement of
  // the scene

  helper->refreshTime();
  std::cout << "\n[[ Step " << iterationCount << ", sim. phys. time: " << iterationCount*timeStep << " s ]]\n";
  // SEARCH FOR NEIGHBOURS PART
  // ocl_solver->_runClearBuffers();
  // helper->watch_report("_runClearBuffers: \t%9.3f ms\n");
  ocl_solver->_runHashParticles(config);
  helper->watch_report("_runHashParticles: \t%9.3f ms\n");
  ocl_solver->_runSort(config);
  helper->watch_report("_runSort: \t\t%9.3f ms\n");
  ocl_solver->_runSortPostPass(config);
  helper->watch_report("_runSortPostPass: \t%9.3f ms\n");
  ocl_solver->_runIndexx(config);
  helper->watch_report("_runIndexx: \t\t%9.3f ms\n");
  ocl_solver->_runIndexPostPass(config);
  helper->watch_report("_runIndexPostPass: \t%9.3f ms\n");
  ocl_solver->_runFindNeighbors(config);
  helper->watch_report("_runFindNeighbors: \t%9.3f ms\n");
  // PCISPH PART
  if (config->getIntegrationMethod() == LEAPFROG) { // in this case we should
                                                    // remmember value of
                                                    // position on stem i - 1
    // Calc next time (t+dt) positions x(t+dt)
    ocl_solver->_run_pcisph_integrate(iterationCount, 0 /*=positions_mode*/,
                                      config);
  }
  ocl_solver->_run_pcisph_computeDensity(config);
  ocl_solver->_run_pcisph_computeForcesAndInitPressure(config, iterationCount);
  ocl_solver->_run_pcisph_computeElasticForces(config);
  do {
    // printf("\n^^^^ iter %d ^^^^\n",iter);
    ocl_solver->_run_pcisph_predictPositions(config);
    ocl_solver->_run_pcisph_predictDensity(config);
    ocl_solver->_run_pcisph_correctPressure(config);
    ocl_solver->_run_pcisph_computePressureForceAcceleration(config);
    iter++;
  } while (iter < maxIteration);

  // and finally calculate v(t+dt)
  if (config->getIntegrationMethod() == LEAPFROG) {
    ocl_solver->_run_pcisph_integrate(iterationCount, 1 /*=velocities_mode*/,
                                      config);
    helper->watch_report("_runPCISPH: \t\t%9.3f ms\t3 iteration(s)\n");
  } else {
    ocl_solver->_run_pcisph_integrate(iterationCount, 2, config);
    helper->watch_report("_runPCISPH: \t\t%9.3f ms\t3 iteration(s)\n");
  }
  // Handling of Interaction with membranes
  if (config->numOfMembranes > 0) {
    ocl_solver->_run_clearMembraneBuffers(config);
    ocl_solver->_run_computeInteractionWithMembranes(config);
    // compute change of coordinates due to interactions with membranes
    ocl_solver->_run_computeInteractionWithMembranes_finalize(config);
    helper->watch_report("membraneHadling: \t%9.3f ms\n");
  }
  // END
  ocl_solver->read_position_buffer(position_cpp, config);
  helper->watch_report("_readBuffer: \t\t%9.3f ms\n");

  // END PCISPH algorithm
  printf("------------------------------------\n");
  printf("_Total_step_time:\t%9.3f ms\n", helper->getElapsedTime());
  printf("------------------------------------\n");
  if (load_to) {
    if (iterationCount == 0) {
      owHelper::loadConfigurationToFile(position_cpp, config,
                                        elasticConnectionsData_cpp,
                                        membraneData_cpp, true);
    } else {
      if (iterationCount % config->getLogStep() == 0) {
        owHelper::loadConfigurationToFile(position_cpp, config, NULL, NULL,
                                          false);
      }
    }
  }
  if (owVtkExport::isActive) {
    if (iterationCount % config->getLogStep() == 0) {
      getvelocity_cpp();
      owVtkExport::exportState(iterationCount, config, position_cpp,
                               elasticConnectionsData_cpp, velocity_cpp,
                               membraneData_cpp, muscle_activation_signal_cpp);
    }
  }

  //==================================================================
  // writing locomotion log
  //==================================================================

  if (iterationCount % 1000 == 0) {
    char f_motion_log_filename[256];
    char f_motion_log_mode[10] = "wt";
    sprintf(f_motion_log_filename, "visc_coeff=%.3e m.m.f.=%4d motion log.txt",
            viscosity, (int)max_muscle_force); // everywhere here log is for
                                               // log-file, not logarithm
    float log_x[200], log_y[200], log_z[200], log_n[200];

    if (iterationCount > 0)
      sprintf(f_motion_log_mode, "a+");

    FILE *f_motion_log = fopen(f_motion_log_filename, f_motion_log_mode);

    fprintf(f_motion_log, "%e\tX:\t",
            (float)iterationCount * timeStep); // timeStep);

    float *ec_cpp = getElasticConnectionsData_cpp();
    float *p_cpp = getPosition_cpp();
    int L_index_i, i, ecc = 0; // elastic connections counter;

    // for(int i_ec=0; i_ec < config->numOfElasticP * MAX_NEIGHBOR_COUNT;
    // i_ec++)
    //{i = (i_ec / MAX_NEIGHBOR_COUNT);

    for (i = 0; i < 200; i++) {
      log_x[i] = log_y[i] = log_z[i] = 0.f;
      log_n[i] = 0;
    }

    for (i = 0; i < config->numOfElasticP; i++) {
      if ((p_cpp[i * 4 + 3] > 2.05f) && (p_cpp[i * 4 + 3] < 2.25f)) {
        if ((p_cpp[i * 4 + 3] > 2.05f) && (p_cpp[i * 4 + 3] < 2.15f))
          L_index_i = (int)((p_cpp[i * 4 + 3] - 2.1f) * 10000.f) + 1.f; //-100;
        else
          L_index_i = (int)((p_cpp[i * 4 + 3] - 2.2f) * 10000.f) + 1.f; //-100;

        if ((L_index_i >= 1) && (L_index_i < 200 /*100*/)) {
          log_x[L_index_i] += p_cpp[i * 4 + 0];
          log_y[L_index_i] += p_cpp[i * 4 + 1];
          log_z[L_index_i] += p_cpp[i * 4 + 2];
          log_n[L_index_i]++;
        }
      }
    }

    for (i = 0; i < 200; i++) {
      if (log_n[i] > 0) {
        log_x[i] /= (float)log_n[i];
        log_y[i] /= (float)log_n[i];
        log_z[i] /= (float)log_n[i];
      }
    }

    for (i = 1; i < 100; i++) {
      fprintf(f_motion_log, "%e\t",
              log_z[i + 50 /*2*/] * simulationScale * 1000.f);
    }
    fprintf(f_motion_log, "\tY:\t");

    for (i = 1; i < 100; i++) {
      fprintf(f_motion_log, "%e\t",
              log_x[i + 50 /*2*/] * simulationScale * 1000.f);
    }
    fprintf(f_motion_log, "\tZ:\t");

    for (i = 1; i < 100; i++) {
      fprintf(f_motion_log, "%e\t",
              log_y[i + 50 /*2*/] * simulationScale * 1000.f);
    }
    fprintf(f_motion_log, "\n");

    fclose(f_motion_log);
  }

  //==================================================================

  float correction_coeff;

  /*for (unsigned int i = 0; i < config->MUSCLE_COUNT; ++i) {
    correction_coeff = sqrt(
        1.f - ((1 + i % 24 - 12.5f) / 12.5f) * ((1 + i % 24 - 12.5f) / 12.5f));
    // printf("\n%d\t%d\t%f\n",i,1+i%24,correction_coeff);
    muscle_activation_signal_cpp[i] *= correction_coeff;
  }*/

  config->updateNeuronSimulation(muscle_activation_signal_cpp);

  int half_i;

  for (int i = 0; i < config->MUSCLE_COUNT; i++) {
    half_i = (i % 24) / 2;
    muscle_activation_signal_cpp[i] *= muscle_activation_signal_cpp[i];
    // if(muscle_activation_signal_cpp[i]<0.5) muscle_activation_signal_cpp[i] =
    // 0;

    // correction_coeff = sqrt( 1.f -
    // ((1+i%24-12.5f)/12.5f)*((1+i%24-12.5f)/12.5f) );

    // correction_coeff = 1.0f;

    /*
    correction_coeff = 1.0f*sqrt( 1.f -
    ((1+half_i-6.5f)/6.5f)*((1+half_i-6.5f)/6.5f)) - 0.0;

    if(half_i>=6)
    {
            correction_coeff = (1.5f*sqrt( 1.f -
    ((1+half_i-6.5f)/6.5f)*((1+half_i-6.5f)/6.5f)) - 0.5f)*0.8f ;
            //if( ((i>=0)&&(i<=3)) || ((i>=20)&&(i<=23)) )
    correction_coeff*=0.5;
            //printf("\n%d\t%d\t%f\n",i,half_i,correction_coeff);
            muscle_activation_signal_cpp[i] *=
    correction_coeff;//*correction_coeff;
    }*/

    // if(half_i>=6)
    /**/
    {
      // correction_coeff = 1.05*sin((half_i+1)/5.5f+0.7f);
      // correction_coeff = 1.05*sin((half_i+1)/5.5f+0.7f);

      if (iterationCount > 1200000)
        correction_coeff = 1.2f * sin((half_i + 0) / 5.4f + 0.95f); // swimming
      else
        correction_coeff = (2500 / 2000) * 1.2f *
                           (1 + sin((half_i + 0) / 5.4f + 0.95f)) /
                           2.f; // crawling
      // printf("\n%d\t%d\t%f\n",i,half_i,correction_coeff);
      muscle_activation_signal_cpp[i] *= correction_coeff;
    } /**/

    /*
    if(iterationCount<5000)
    muscle_activation_signal_cpp[i] = 0.0f;//1.f;
    else
    muscle_activation_signal_cpp[i] = 0.75f;//1.f;
    /**/
  }

  /**/
  if (iterationCount < 10000) {
    for (int i = 0; i < config->MUSCLE_COUNT; i++) {
      muscle_activation_signal_cpp[i] *= (float)iterationCount / 10000.f;
    }
  } /**/

  ocl_solver->updateMuscleActivityData(muscle_activation_signal_cpp, config);
  iterationCount++;
  return helper->getElapsedTime();
}
/** Prepare data and log it to special configuration
 *  file you can run your simulation from place you snapshoted it
 *
 *  @param fileName - name of file where saved configuration will be stored
 */
void owPhysicsFluidSimulator::makeSnapshot() {
  getvelocity_cpp();
  std::string fileName = config->getSnapshotFileName();
  owHelper::loadConfigurationToFile(
      position_cpp, velocity_cpp, elasticConnectionsData_cpp, membraneData_cpp,
      particleMembranesList_cpp, fileName.c_str(), config);
}

// Destructor
owPhysicsFluidSimulator::~owPhysicsFluidSimulator(void) {
  destroy();
  delete config;
  delete ocl_solver;
}

void owPhysicsFluidSimulator::destroy() {
  delete[] position_cpp;
  delete[] velocity_cpp;
  delete[] density_cpp;
  delete[] particleIndex_cpp;
  if (elasticConnectionsData_cpp != NULL)
    delete[] elasticConnectionsData_cpp;
  delete[] muscle_activation_signal_cpp;
  if (membraneData_cpp != NULL) {
    delete[] membraneData_cpp;
    delete[] particleMembranesList_cpp;
  }
}
