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
#include <iomanip>

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
    gettimeofday(&simulation_start, NULL);
    iterationCount = 0;
    config = new owConfigProperty(argc, argv);
    // LOAD FROM FILE
    config->initGridCells();
    position_cpp = new float[4 * config->getParticleCount()];
    velocity_cpp = new float[4 * config->getParticleCount()];
    pressure_cpp = new float[1 * config->getParticleCount()];
    muscle_activation_signal_cpp = new float[config->MUSCLE_COUNT];
    if (config->numOfElasticP != 0)
      elasticConnectionsData_cpp =
          new float[4 * config->numOfElasticP * MAX_NEIGHBOR_COUNT];
    if (config->numOfMembranes <= 0)
      membraneData_cpp = nullptr;
    else
      membraneData_cpp = new int[config->numOfMembranes * 3];
    if (config->numOfElasticP <= 0)
      particleMembranesList_cpp = nullptr;
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
    this->genShellPaticlesList();
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
  config->initGridCells();
  position_cpp = new float[4 * config->getParticleCount()];
  velocity_cpp = new float[4 * config->getParticleCount()];
  pressure_cpp = new float[1 * config->getParticleCount()];
  muscle_activation_signal_cpp = new float[config->MUSCLE_COUNT];
  if (config->numOfElasticP != 0)
    elasticConnectionsData_cpp =
        new float[4 * config->numOfElasticP * MAX_NEIGHBOR_COUNT];
  if (config->numOfMembranes <= 0)
    membraneData_cpp = nullptr;
  else
    membraneData_cpp = new int[config->numOfMembranes * 3];
  if (config->numOfElasticP <= 0)
    particleMembranesList_cpp = nullptr;
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
  this->genShellPaticlesList();
}
int update_muscle_activity_signals_log_file(int iterationCount,
                                            float *muscle_activation_signal_cpp,
                                            owConfigProperty *config) {
  /*	char muscle_log_file_mode [10] = "wt";

            if(iterationCount > 0 ) sprintf(muscle_log_file_mode,"a+");

            FILE *f_muscle_log =
     fopen(muscle_log_file_name,muscle_log_file_mode);

            if(!f_muscle_log) return -1;// Can't open file

            fprintf(f_muscle_log,"%10d\t%e\t",iterationCount,(float)iterationCount*timeStep);

            for(int i_m=0;i_m<=95;i_m++)
            {
                    fprintf(f_muscle_log,"%.3f\t",muscle_activation_signal_cpp[i_m]);
            }

            fprintf(f_muscle_log,"\n");

            fclose(f_muscle_log);*/

  // write log file with muscular activity data

  std::ofstream musclesActivityFile;
  std::string musclesActivityFileName =
      config->getLoadPath() + std::string("/muscles_activity_buffer.txt");

  if (iterationCount == 0) {
    musclesActivityFile.open(musclesActivityFileName.c_str(),
                             std::ofstream::trunc);
    if (!musclesActivityFile)
      throw std::runtime_error("There was a problem with creation of muscles "
                               "activity file for logging. Please check the "
                               "path.");
  } else {
    musclesActivityFile.open(musclesActivityFileName.c_str(),
                             std::ofstream::app);
    if (!musclesActivityFile)
      throw std::runtime_error("There was a problem with creation of muscles "
                               "activity file for logging. Please check the "
                               "path.");
  }

  musclesActivityFile << std::scientific;

  if ((muscle_activation_signal_cpp != nullptr)) {
    for (unsigned int i = 0; i < config->MUSCLE_COUNT; i++) {
      if (i == config->MUSCLE_COUNT - 1)
        musclesActivityFile << muscle_activation_signal_cpp[i] << "\n";
      else
        musclesActivityFile << muscle_activation_signal_cpp[i] << "\t";
    }
  }

  musclesActivityFile.close();

  return 0;
}

int update_worm_motion_log_file(
    int iterationCount, float *ec_cpp /*getElasticConnectionsData_cpp()*/,
    float *p_cpp /*getPosition_cpp()*/, owConfigProperty *config) {
  //	char motion_log_file_mode [10] = "wt";
  float log_x[200], log_y[200], log_z[200], log_n[200];
  // float * ec_cpp = getElasticConnectionsData_cpp();
  // float * p_cpp = getPosition_cpp();
  int L_index_i, i; // elastic connections counter;

  std::ofstream wormMotionLogFile;
  std::string wormMotionLogFileName =
      config->getLoadPath() + std::string("/worm_motion_log.txt");

  if (iterationCount == 0) {
    wormMotionLogFile.open(wormMotionLogFileName.c_str(), std::ofstream::trunc);
    if (!wormMotionLogFile)
      throw std::runtime_error("There was a problem with creation of worm "
                               "motion log file. Check the path.");
  } else {
    wormMotionLogFile.open(wormMotionLogFileName.c_str(), std::ofstream::app);
    if (!wormMotionLogFile)
      throw std::runtime_error("There was a problem with creation of worm "
                               "motion log file. Check the path.");
  }

  for (i = 0; i < 200; i++) {
    log_x[i] = log_y[i] = log_z[i] = 0.f;
    log_n[i] = 0;
  }

  // fprintf(f_motion_log,"%e\tX:\t",(float)iterationCount*timeStep);
  wormMotionLogFile << std::scientific;
  wormMotionLogFile << (float)iterationCount * config->getTimeStep()
                    << "\tX:\t";
  // wormMotionLogFile << (float)iterationCount*timeStep << "\tX:\t";

  for (i = 0; (unsigned)i < config->numOfElasticP; i++) {
    if ((p_cpp[i * 4 + 3] > 2.05f) && (p_cpp[i * 4 + 3] < 2.25f)) {
      if ((p_cpp[i * 4 + 3] > 2.05f) && (p_cpp[i * 4 + 3] < 2.15f)) {
        L_index_i = (int)((p_cpp[i * 4 + 3] + 0.0000003f - 2.1000f) * 10000.f) +
                    1; //-100;
      } else {
        L_index_i = (int)((p_cpp[i * 4 + 3] + 0.0000003f - 2.2000f) * 10000.f) +
                    1; //-100;
      }

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
    // fprintf(f_motion_log,"%e\t",log_z[i+50/*2*/]*simulationScale*1000.f);
    wormMotionLogFile << log_z[i + 50 /*2*/] * simulationScale * 1000.f << "\t";
  } // fprintf(f_motion_log,"\tY:\t");
  wormMotionLogFile << "\tY:\t";

  for (i = 1; i < 100; i++) {
    // fprintf(f_motion_log,"%e\t",log_x[i+50/*2*/]*simulationScale*1000.f);
    wormMotionLogFile << log_x[i + 50 /*2*/] * simulationScale * 1000.f << "\t";
  } // fprintf(f_motion_log,"\tZ:\t");
  wormMotionLogFile << "\tZ:\t";

  for (i = 1; i < 100; i++) {
    // fprintf(f_motion_log,"%e\t",log_y[i+50/*2*/]*simulationScale*1000.f);
    wormMotionLogFile << log_y[i + 50 /*2*/] * simulationScale * 1000.f << "\t";
  } // fprintf(f_motion_log,"\n");
  wormMotionLogFile << "\n";

  // fclose(f_motion_log);
  wormMotionLogFile.close();

  return 0;
}

/** Gen list of particles in shell
 */
void owPhysicsFluidSimulator::genShellPaticlesList() {
  for (size_t i = 0; i < config->numOfElasticP; ++i) {
    for (size_t j = 0; j < MAX_MEMBRANES_INCLUDING_SAME_PARTICLE; ++j) {
      if (particleMembranesList_cpp[i * MAX_MEMBRANES_INCLUDING_SAME_PARTICLE +
                                    j] != -1) {
        shellIndexes.push_back(i);
        break;
      }
    }
  }
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
  std::cout << "\n[[ Step " << iterationCount << " (total steps: ";
  if (config->getNumberOfIterations() == 0)
    std::cout << "unlimited";
  else
    std::cout << config->getNumberOfIterations();

  std::cout << ", t in sim: " << iterationCount * config->getTimeStep()
            << "s) dt: " << config->getTimeStep() << " (in s)";

  struct timeval current_time;
  gettimeofday(&current_time, NULL);
  float elapsed_seconds = (float)(current_time.tv_sec - simulation_start.tv_sec);
  float time_elapsed = elapsed_seconds / 60.0;
  std::string time_elapsed_unit = "(in min)";
  if (time_elapsed > 60.0) {
    time_elapsed /= 60.0;
    time_elapsed_unit = "(in h)";
  }
  printf(", time elapsed: %.2f %s", time_elapsed, time_elapsed_unit.c_str());
  if (config->getNumberOfIterations() > 0) {
    int steps_left = config->getNumberOfIterations() - iterationCount;
    float time_left = ((elapsed_seconds/iterationCount)*steps_left/60.0);
    std::string time_left_unit = "(in min)";
    if (time_left > 60) {
        time_left /= 60.0;
        time_left_unit = "(in h)";
    }
    printf(", time left: %.2f %s", time_left, time_left_unit.c_str());
  }
  std::cout << " ]]\n";

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
  ocl_solver->_run_pcisph_computeForcesAndInitPressure(config);
  ocl_solver->_run_pcisph_computeElasticForces(config);
  do {
    // printf("\n^^^^ iter %d ^^^^\n",iter);
    ocl_solver->_run_pcisph_predictPositions(config);
    ocl_solver->_run_pcisph_predictDensity(config);
    ocl_solver->_run_pcisph_correctPressure(config);
    ocl_solver->_run_pcisph_computePressureForceAcceleration(config);
    iter++;
  } while (iter < config->getConst("maxIteration"));

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
    helper->watch_report("membraneHandling: \t%9.3f ms\n");
  }
  // END
  ocl_solver->read_position_buffer(position_cpp, config);
  ocl_solver->read_pressure_buffer(pressure_cpp, config);
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
      owHelper::loadPressureToFile(pressure_cpp, shellIndexes, position_cpp, iterationCount, config);
    } else {
      if (iterationCount % config->getLogStep() == 0) {
        owHelper::loadConfigurationToFile(position_cpp, config, nullptr,
                                          nullptr, false);
        owHelper::loadPressureToFile(pressure_cpp, shellIndexes, position_cpp, iterationCount, config);
      }
    }
  }
  if (owVtkExport::isActive) {
    if (iterationCount % config->getLogStep() == 0) {
      getVelocity_cpp();
      owVtkExport::exportState(iterationCount, config, position_cpp,
                               elasticConnectionsData_cpp, velocity_cpp,
                               membraneData_cpp, muscle_activation_signal_cpp);
    }
  }

  config->updateNeuronSimulation(muscle_activation_signal_cpp);
  /* // signal correction switched off
    float correction_coeff;

    for (unsigned int i = 0; i < config->MUSCLE_COUNT; ++i) {
      correction_coeff = sqrt(
          1.f - ((1 + i % 24 - 12.5f) / 12.5f) * ((1 + i % 24 - 12.5f)
  / 12.5f));
      // printf("\n%d\t%d\t%f\n",i,1+i%24,correction_coeff);
      //muscle_activation_signal_cpp[i] *= correction_coeff;
            muscle_activation_signal_cpp[i] *= muscle_activation_signal_cpp[i];
            muscle_activation_signal_cpp[i] *= 1.0f*(1.f-0.4f*(i%24)/24.f);
          //if(muscle_activation_signal_cpp[i]<0.25f*0.9f)
  muscle_activation_signal_cpp[i] = 0.f;
    }



   //smooth start switched off
        if(iterationCount<5000)
        {
                for(int i=0;i<config->MUSCLE_COUNT;i++)
                {
                        muscle_activation_signal_cpp[i] *=
     (float)iterationCount/5000.f;
                }
        }*/

  if (iterationCount % config->getLogStep() == 0) {
    update_muscle_activity_signals_log_file(
        iterationCount, muscle_activation_signal_cpp, config);
    update_worm_motion_log_file(iterationCount, elasticConnectionsData_cpp,
                                position_cpp,
                                config); // format: time, x_1...x_n, y_1...y_n,
                                         // z_1...z_n (worm body central line
                                         // coordinates, from head to tail)
  }

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
  getVelocity_cpp();
  std::string fileName = config->getSnapshotFileName();
  owHelper::loadConfigurationToFile(
      position_cpp, velocity_cpp, elasticConnectionsData_cpp, membraneData_cpp,
      particleMembranesList_cpp, fileName.c_str(),
      config /*, nullptr/muscle_activation_signal_cpp*/);
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
  if (elasticConnectionsData_cpp != nullptr)
    delete[] elasticConnectionsData_cpp;
  delete[] muscle_activation_signal_cpp;
  if (membraneData_cpp != nullptr) {
    delete[] membraneData_cpp;
    delete[] particleMembranesList_cpp;
  }
}
