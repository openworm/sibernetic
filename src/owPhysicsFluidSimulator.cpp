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
#include <cstring>  // for memcpy
#include <Python.h>

// NumPy C API for fast array access
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "owPhysicsFluidSimulator.h"
#include "owSignalSimulator.h"
#include "owVtkExport.h"

// Make pytorch_solver.py / taichi_solver.py and packages installed in the
// project-local .venv importable from the embedded interpreter without the
// user having to export PYTHONPATH first. Assumes the binary is launched from
// the repo root, which matches the existing OpenCL kernel-path convention.
static void setupEmbeddedPythonPath() {
  PyRun_SimpleString(
      "import sys, os, glob\n"
      "_cwd = os.getcwd()\n"
      "if _cwd not in sys.path:\n"
      "    sys.path.insert(0, _cwd)\n"
      "for _p in glob.glob(os.path.join(_cwd, '.venv', 'lib', 'python*', 'site-packages')):\n"
      "    if _p not in sys.path:\n"
      "        sys.path.append(_p)\n"
  );
}

// Phase 1.4: Helper to eliminate duplicate init code between constructor and reset()
void owPhysicsFluidSimulator::initTorchSolver(bool isReset) {
  // Clean up existing solver if this is a reset
  if (isReset && torchSolver) {
    Py_DECREF(torchSolver);
    torchSolver = nullptr;
  }

  // Load module and class
  PyObject *pName = PyUnicode_FromString("pytorch_solver");
  PyObject *pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  if (!pModule) {
    PyErr_Print();
    throw std::runtime_error("Failed to load pytorch_solver module");
  }
  PyObject *pClass = PyObject_GetAttrString(pModule, "PytorchSolver");
  Py_DECREF(pModule);
  if (!pClass) {
    PyErr_Print();
    throw std::runtime_error("Failed to get PytorchSolver class");
  }

  // Build position and velocity lists
  int count = config->getParticleCount();
  PyObject *pos_list = PyList_New(count);
  PyObject *vel_list = PyList_New(count);
  for (int i = 0; i < count; ++i) {
    PyObject *row_p = PyList_New(4);
    PyObject *row_v = PyList_New(4);
    for (int j = 0; j < 4; ++j) {
      PyList_SetItem(row_p, j, PyFloat_FromDouble(position_cpp[i * 4 + j]));
      PyList_SetItem(row_v, j, PyFloat_FromDouble(velocity_cpp[i * 4 + j]));
    }
    PyList_SetItem(pos_list, i, row_p);
    PyList_SetItem(vel_list, i, row_v);
  }

  // Build config dictionary with proper reference counting
  PyObject *cfg = PyDict_New();
  PyObject *val;

  val = PyFloat_FromDouble(config->xmin);
  PyDict_SetItemString(cfg, "xmin", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->ymin);
  PyDict_SetItemString(cfg, "ymin", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->zmin);
  PyDict_SetItemString(cfg, "zmin", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("hashGridCellSizeInv"));
  PyDict_SetItemString(cfg, "hash_grid_cell_size_inv", val); Py_DECREF(val);
  val = PyLong_FromUnsignedLong(config->gridCellsX);
  PyDict_SetItemString(cfg, "grid_cells_x", val); Py_DECREF(val);
  val = PyLong_FromUnsignedLong(config->gridCellsY);
  PyDict_SetItemString(cfg, "grid_cells_y", val); Py_DECREF(val);
  val = PyLong_FromUnsignedLong(config->gridCellsZ);
  PyDict_SetItemString(cfg, "grid_cells_z", val); Py_DECREF(val);
  val = PyLong_FromUnsignedLong(config->gridCellCount);
  PyDict_SetItemString(cfg, "grid_cell_count", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("h"));
  PyDict_SetItemString(cfg, "h", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("mass_mult_Wpoly6Coefficient"));
  PyDict_SetItemString(cfg, "mass_mult_Wpoly6Coefficient", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("mass_mult_gradWspikyCoefficient"));
  PyDict_SetItemString(cfg, "mass_mult_gradWspikyCoefficient", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("rho0"));
  PyDict_SetItemString(cfg, "rho0", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getDelta());
  PyDict_SetItemString(cfg, "delta", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getTimeStep());
  PyDict_SetItemString(cfg, "time_step", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("gravity_x"));
  PyDict_SetItemString(cfg, "gravity_x", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("gravity_y"));
  PyDict_SetItemString(cfg, "gravity_y", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("gravity_z"));
  PyDict_SetItemString(cfg, "gravity_z", val); Py_DECREF(val);
  // Phase 1.2: Additional config parameters for PyTorch solver parity
  val = PyLong_FromLong(MAX_NEIGHBOR_COUNT);
  PyDict_SetItemString(cfg, "max_neighbor_count", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("simulationScale"));
  PyDict_SetItemString(cfg, "simulation_scale", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("simulationScaleInv"));
  PyDict_SetItemString(cfg, "simulation_scale_inv", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("r0"));
  PyDict_SetItemString(cfg, "r0", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("viscosity"));
  PyDict_SetItemString(cfg, "viscosity", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("mass_mult_divgradWviscosityCoefficient"));
  PyDict_SetItemString(cfg, "mass_mult_divgradWviscosityCoefficient", val); Py_DECREF(val);
  val = PyLong_FromLong(3);  // PCISPH uses 3 iterations
  PyDict_SetItemString(cfg, "max_iteration", val); Py_DECREF(val);
  val = PyUnicode_FromString(config->getTorchDevice().c_str());
  PyDict_SetItemString(cfg, "device", val); Py_DECREF(val);

  // Create solver instance
  PyObject *args = PyTuple_Pack(3, pos_list, vel_list, cfg);
  torchSolver = PyObject_CallObject(pClass, args);
  Py_DECREF(args);
  Py_DECREF(pos_list);
  Py_DECREF(vel_list);
  Py_DECREF(cfg);
  Py_DECREF(pClass);

  if (!torchSolver) {
    PyErr_Print();
    throw std::runtime_error("Failed to create PytorchSolver instance");
  }
}

// Initialize Taichi GPU solver
void owPhysicsFluidSimulator::initTaichiSolver(bool isReset) {
  // Clean up existing solver if this is a reset
  if (isReset && taichiSolver) {
    Py_DECREF(taichiSolver);
    taichiSolver = nullptr;
  }

  // Load module and class
  PyObject *pName = PyUnicode_FromString("taichi_solver");
  PyObject *pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  if (!pModule) {
    PyErr_Print();
    throw std::runtime_error("Failed to load taichi_solver module");
  }
  PyObject *pClass = PyObject_GetAttrString(pModule, "TaichiSolver");
  Py_DECREF(pModule);
  if (!pClass) {
    PyErr_Print();
    throw std::runtime_error("Failed to get TaichiSolver class");
  }

  // Build position and velocity lists
  int count = config->getParticleCount();
  PyObject *pos_list = PyList_New(count);
  PyObject *vel_list = PyList_New(count);
  for (int i = 0; i < count; ++i) {
    PyObject *row_p = PyList_New(4);
    PyObject *row_v = PyList_New(4);
    for (int j = 0; j < 4; ++j) {
      PyList_SetItem(row_p, j, PyFloat_FromDouble(position_cpp[i * 4 + j]));
      PyList_SetItem(row_v, j, PyFloat_FromDouble(velocity_cpp[i * 4 + j]));
    }
    PyList_SetItem(pos_list, i, row_p);
    PyList_SetItem(vel_list, i, row_v);
  }

  // Build config dictionary
  PyObject *cfg = PyDict_New();
  PyObject *val;

  // Grid boundaries
  val = PyFloat_FromDouble(config->xmin);
  PyDict_SetItemString(cfg, "xmin", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->xmax);
  PyDict_SetItemString(cfg, "xmax", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->ymin);
  PyDict_SetItemString(cfg, "ymin", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->ymax);
  PyDict_SetItemString(cfg, "ymax", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->zmin);
  PyDict_SetItemString(cfg, "zmin", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->zmax);
  PyDict_SetItemString(cfg, "zmax", val); Py_DECREF(val);

  // SPH parameters
  val = PyFloat_FromDouble(config->getConst("h"));
  PyDict_SetItemString(cfg, "h", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("rho0"));
  PyDict_SetItemString(cfg, "rho0", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getDelta());
  PyDict_SetItemString(cfg, "delta", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getTimeStep());
  PyDict_SetItemString(cfg, "time_step", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("mass"));
  PyDict_SetItemString(cfg, "mass", val); Py_DECREF(val);

  // Gravity
  val = PyFloat_FromDouble(config->getConst("gravity_x"));
  PyDict_SetItemString(cfg, "gravity_x", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("gravity_y"));
  PyDict_SetItemString(cfg, "gravity_y", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("gravity_z"));
  PyDict_SetItemString(cfg, "gravity_z", val); Py_DECREF(val);

  // Additional parameters
  val = PyFloat_FromDouble(config->getConst("r0"));
  PyDict_SetItemString(cfg, "r0", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("viscosity"));
  PyDict_SetItemString(cfg, "mu", val); Py_DECREF(val);
  val = PyLong_FromLong(MAX_NEIGHBOR_COUNT);
  PyDict_SetItemString(cfg, "max_neighbor_count", val); Py_DECREF(val);

  // Floor constraint
  val = PyFloat_FromDouble(0.0);  // floor_y
  PyDict_SetItemString(cfg, "floor_y", val); Py_DECREF(val);
  val = PyFloat_FromDouble(0.3);  // floor_restitution
  PyDict_SetItemString(cfg, "floor_restitution", val); Py_DECREF(val);

  // Device (metal, cuda, cpu)
  val = PyUnicode_FromString(config->getTaichiDevice().c_str());
  PyDict_SetItemString(cfg, "device", val); Py_DECREF(val);

  // Simulation scale (CRITICAL for elastic force calculations!)
  val = PyFloat_FromDouble(config->getConst("simulationScale"));
  PyDict_SetItemString(cfg, "simulation_scale", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("simulationScaleInv"));
  PyDict_SetItemString(cfg, "simulation_scale_inv", val); Py_DECREF(val);

  // Elastic parameters
  val = PyLong_FromUnsignedLong(config->numOfElasticP);
  PyDict_SetItemString(cfg, "num_elastic_particles", val); Py_DECREF(val);
  val = PyFloat_FromDouble(config->getConst("elasticityCoefficient"));
  PyDict_SetItemString(cfg, "elasticity_coefficient", val); Py_DECREF(val);

  // Build elastic connections list (if we have elastic particles)
  PyObject *elastic_list = Py_None;
  Py_INCREF(Py_None);
  if (config->numOfElasticP > 0 && elasticConnectionsData_cpp != nullptr) {
    Py_DECREF(elastic_list);
    // Each elastic particle has MAX_NEIGHBOR_COUNT connections
    // Each connection is 4 floats: [particle_id, rest_length, muscle_id, unused]
    int num_connections = config->numOfElasticP * MAX_NEIGHBOR_COUNT;
    elastic_list = PyList_New(num_connections);
    for (int i = 0; i < num_connections; ++i) {
      PyObject *conn = PyList_New(4);
      for (int j = 0; j < 4; ++j) {
        PyList_SetItem(conn, j, PyFloat_FromDouble(elasticConnectionsData_cpp[i * 4 + j]));
      }
      PyList_SetItem(elastic_list, i, conn);
    }
  }

  // Create solver instance with elastic connections
  PyObject *args = PyTuple_Pack(4, pos_list, vel_list, cfg, elastic_list);
  taichiSolver = PyObject_CallObject(pClass, args);
  Py_DECREF(args);
  Py_DECREF(pos_list);
  Py_DECREF(vel_list);
  Py_DECREF(cfg);
  Py_DECREF(elastic_list);
  Py_DECREF(pClass);

  if (!taichiSolver) {
    PyErr_Print();
    throw std::runtime_error("Failed to create TaichiSolver instance");
  }
}

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
    torchSolver = nullptr;
    taichiSolver = nullptr;
    useTorchBackend = config->torchEnabled();
    useTaichiBackend = config->taichiEnabled();

    if (useTorchBackend) {
      Py_Initialize();
      setupEmbeddedPythonPath();  // Add cwd + .venv site-packages to sys.path
      _import_array();  // Initialize NumPy C API for fast array access (no return)
      initTorchSolver(false);  // Phase 1.4: Use helper (not a reset)
#ifndef OW_NO_OPENCL
      ocl_solver = nullptr;
#endif
    } else if (useTaichiBackend) {
      Py_Initialize();
      setupEmbeddedPythonPath();  // Add cwd + .venv site-packages to sys.path
      _import_array();  // Initialize NumPy C API for fast array access (no return)
      initTaichiSolver(false);  // Initialize Taichi GPU solver
#ifndef OW_NO_OPENCL
      ocl_solver = nullptr;
#endif
    }
#ifndef OW_NO_OPENCL
    else {
      if (config->numOfElasticP != 0) {
        ocl_solver = new owOpenCLSolver(
            position_cpp, velocity_cpp, config, elasticConnectionsData_cpp,
            membraneData_cpp,
            particleMembranesList_cpp); // Create new openCLsolver instance
      } else
        ocl_solver = new owOpenCLSolver(position_cpp, velocity_cpp,
                                        config); // Create new openCLsolver instance
    }
#endif
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
  if (useTorchBackend) {
    initTorchSolver(true);  // Phase 1.4: Use helper (is a reset)
  } else if (useTaichiBackend) {
    initTaichiSolver(true);  // Reset Taichi solver
  }
#ifndef OW_NO_OPENCL
  else {
    if (config->numOfElasticP != 0) {
      ocl_solver->reset(
          position_cpp, velocity_cpp, config, elasticConnectionsData_cpp,
          membraneData_cpp,
          particleMembranesList_cpp); // Create new openCLsolver instance
    } else
      ocl_solver->reset(position_cpp, velocity_cpp,
                        config); // Create new openCLsolver instance
  }
#endif
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

  if (useTorchBackend) {
    // Helper lambda for error-checked method calls
    auto callMethod = [this](const char* method) {
      PyObject *result = PyObject_CallMethod(torchSolver, method, nullptr);
      if (!result) {
        PyErr_Print();
        throw std::runtime_error(std::string("PyTorch solver method failed: ") + method);
      }
      Py_DECREF(result);
    };

    callMethod("run_hash_particles");
    callMethod("run_sort");
    callMethod("run_index");
    callMethod("run_index_post_pass");
    callMethod("run_find_neighbors");
    callMethod("run_compute_density");
    callMethod("run_compute_pressure");
    callMethod("run_compute_pressure_force_acceleration");
    callMethod("run_integrate");
    // Get state as numpy arrays (optimized - no Python list conversion)
    PyObject *state = PyObject_CallMethod(torchSolver, "get_state", nullptr);
    if (state && PyTuple_Check(state)) {
      PyArrayObject *pos_arr = (PyArrayObject*)PyTuple_GetItem(state, 0);
      PyArrayObject *vel_arr = (PyArrayObject*)PyTuple_GetItem(state, 1);

      // Fast path: direct memcpy if arrays are contiguous float32
      if (PyArray_Check(pos_arr) && PyArray_Check(vel_arr) &&
          PyArray_IS_C_CONTIGUOUS(pos_arr) && PyArray_IS_C_CONTIGUOUS(vel_arr) &&
          PyArray_TYPE(pos_arr) == NPY_FLOAT32 && PyArray_TYPE(vel_arr) == NPY_FLOAT32) {
        float *pos_data = (float*)PyArray_DATA(pos_arr);
        float *vel_data = (float*)PyArray_DATA(vel_arr);
        size_t bytes = config->getParticleCount() * 4 * sizeof(float);
        std::memcpy(position_cpp, pos_data, bytes);
        std::memcpy(velocity_cpp, vel_data, bytes);
      } else {
        // Fallback: slow path for non-contiguous or non-float32 arrays
        for (Py_ssize_t i = 0; i < PyArray_DIM(pos_arr, 0); ++i) {
          for (int j = 0; j < 4; ++j) {
            position_cpp[i * 4 + j] = *(float*)PyArray_GETPTR2(pos_arr, i, j);
            velocity_cpp[i * 4 + j] = *(float*)PyArray_GETPTR2(vel_arr, i, j);
          }
        }
      }
    }
    Py_XDECREF(state);
    if (load_to) {
      if (iterationCount == 0) {
        owHelper::loadConfigurationToFile(position_cpp, config,
                                          elasticConnectionsData_cpp,
                                          membraneData_cpp, true);
        owHelper::loadVelocityToFile(velocity_cpp, iterationCount, config);
      } else if (iterationCount % config->getLogStep() == 0) {
        owHelper::loadConfigurationToFile(position_cpp, config, nullptr,
                                          nullptr, false);
        owHelper::loadVelocityToFile(velocity_cpp, iterationCount, config);
      }
    }
    iterationCount++;
    return helper->getElapsedTime();
  }

  if (useTaichiBackend) {
    // Taichi GPU solver path
    auto callMethod = [this](const char* method) {
      PyObject *result = PyObject_CallMethod(taichiSolver, method, nullptr);
      if (!result) {
        PyErr_Print();
        throw std::runtime_error(std::string("Taichi solver method failed: ") + method);
      }
      Py_DECREF(result);
    };

    // Run full simulation step on GPU
    callMethod("run_step");

    // Get updated state back to CPU (optimized - no Python list conversion)
    PyObject *state = PyObject_CallMethod(taichiSolver, "get_state", nullptr);
    if (state && PyTuple_Check(state)) {
      PyArrayObject *pos_arr = (PyArrayObject*)PyTuple_GetItem(state, 0);
      PyArrayObject *vel_arr = (PyArrayObject*)PyTuple_GetItem(state, 1);

      // Fast path: direct memcpy if arrays are contiguous float32
      if (PyArray_Check(pos_arr) && PyArray_Check(vel_arr) &&
          PyArray_IS_C_CONTIGUOUS(pos_arr) && PyArray_IS_C_CONTIGUOUS(vel_arr) &&
          PyArray_TYPE(pos_arr) == NPY_FLOAT32 && PyArray_TYPE(vel_arr) == NPY_FLOAT32) {
        float *pos_data = (float*)PyArray_DATA(pos_arr);
        float *vel_data = (float*)PyArray_DATA(vel_arr);
        size_t bytes = config->getParticleCount() * 4 * sizeof(float);
        std::memcpy(position_cpp, pos_data, bytes);
        std::memcpy(velocity_cpp, vel_data, bytes);
      } else {
        // Fallback: slow path for non-contiguous or non-float32 arrays
        for (Py_ssize_t i = 0; i < PyArray_DIM(pos_arr, 0); ++i) {
          for (int j = 0; j < 4; ++j) {
            position_cpp[i * 4 + j] = *(float*)PyArray_GETPTR2(pos_arr, i, j);
            velocity_cpp[i * 4 + j] = *(float*)PyArray_GETPTR2(vel_arr, i, j);
          }
        }
      }
    }
    Py_XDECREF(state);

    // Handle logging
    if (load_to) {
      if (iterationCount == 0) {
        owHelper::loadConfigurationToFile(position_cpp, config,
                                          elasticConnectionsData_cpp,
                                          membraneData_cpp, true);
        owHelper::loadVelocityToFile(velocity_cpp, iterationCount, config);
      } else if (iterationCount % config->getLogStep() == 0) {
        owHelper::loadConfigurationToFile(position_cpp, config, nullptr,
                                          nullptr, false);
        owHelper::loadVelocityToFile(velocity_cpp, iterationCount, config);
      }
    }
    iterationCount++;
    return helper->getElapsedTime();
  }

#ifndef OW_NO_OPENCL
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
#else
  (void)load_to;
  throw std::runtime_error("OpenCL backend disabled");
#endif
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
#ifndef OW_NO_OPENCL
  if (ocl_solver)
    delete ocl_solver;
#endif
  if (torchSolver) {
    Py_DECREF(torchSolver);
  }
  if (taichiSolver) {
    Py_DECREF(taichiSolver);
  }
  if (torchSolver || taichiSolver) {
    Py_Finalize();
  }
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
