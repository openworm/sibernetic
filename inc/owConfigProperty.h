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
/*
 * This file contains definition for struct configuration
 */
#ifndef OWCONFIGURATION_H_
#define OWCONFIGURATION_H_

#include "owNeuronSimulator.h"
#include "owOpenCLConstant.h"
#include "owPhysicsConstant.h"
#include "owSignalSimulator.h"
#include <algorithm>
#include <ctime>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class param_error : public std::runtime_error {
public:
  param_error(const std::string &err_msg) : std::runtime_error(err_msg) {}
};

struct owConfigProperty {
  // This value defines boundary of box in which simulation is
  // Sizes of the box containing simulated 'world'
  // Sizes choice is realized this way because it should be proportional to
  // smoothing radius h
public:
  typedef unsigned int uint;
  const int getParticleCount() { return PARTICLE_COUNT; }
  void setParticleCount(int value) {
    PARTICLE_COUNT = value;
    PARTICLE_COUNT_RoundedUp =
        (((PARTICLE_COUNT - 1) / local_NDRange_size) + 1) * local_NDRange_size;
  }
  float getConst(const std::string &name) /* throw(param_error) */ {
    if (constMap.find(name) == constMap.end())
      throw param_error(std::string("No param with name...") + name);
    return constMap[name];
  }
  void setConst(std::pair<std::string, float> new_const) {
    constMap.insert(new_const);
  }
  void setDeviceType(DEVICE type) { prefDeviceType = type; }
  const int getParticleCount_RoundUp() { return PARTICLE_COUNT_RoundedUp; }
  const int getDeviceType() const { return prefDeviceType; };
  const int getNumberOfIterations() const { return totalNumberOfIterations; }
  const char *getDeviceName() const { return devFullName.c_str(); }
  const std::string &getSourceFileName() const { return sourceFileName; }
  void setDeviceName(const char *name) { devFullName = name; }
  INTEGRATOR getIntegrationMethod() const { return integrationMethod; }
  const std::string &getConfigFileName() const { return configFileName; }
  const std::string &getConfigPath() const { return path; }
  const std::string &getLoadPath() const { return loadPath; }
  // SignalSimulator & getPyramidalSimulation() { return simulation; }
  owINeuronSimulator *getPyramidalSimulation() { return simulation; }
  void updateNeuronSimulation(float *muscleActivationSignal) {
    if (isWormConfig() || nrnSimRun) {
      std::vector<float> muscle_vector = simulation->run();
      for (unsigned int index = 0; index < muscle_vector.size(); index++) {
        muscleActivationSignal[index] = muscle_vector[index];
      }
    }
  }
  bool isWormConfig() {
    return (configFileName.find("worm") != std::string::npos);
  } //   == "worm" || configFileName == "worm_no_water")? true:false; }
  bool isC302() {
    // printf("\nisC302? nc302=%d\n",c302==true);
    return c302;
  }
  void setConfigFileName(const char *name) { configFileName = name; }
  void resetNeuronSimulation() {
    if (isWormConfig() || nrnSimRun) {
      delete simulation;
      if (!isWormConfig()) {
        simulation = new SignalSimulator();
      } else
        simulation =
            new owNeuronSimulator(1, this->timeStep, nrnSimulationFileName);
    }
  }
  void setTimeStep(float value) { this->timeStep = value; }
  void setLogStep(int value) { logStep = value; }
  int getLogStep() { return logStep; }
  float getTimeStep() const { return this->timeStep; }
  float getDelta() const { return delta; }
  std::string getSnapshotFileName() {
    std::string fileName = "./configuration/snapshot/" + configFileName + "_";
    std::stringstream ss;
    time_t t = time(0); // get time now
    struct tm *now = localtime(&t);
    ss << now->tm_hour;
    ss << "-";
    ss << now->tm_min;
    ss << "-";
    ss << now->tm_sec;
    ss << "_";
    ss << now->tm_mday;
    ss << ".";
    ss << (now->tm_mon + 1);
    ss << ".";
    ss << (now->tm_year + 1900);
    fileName += ss.str();
    return fileName;
  }
  // Constructor
  owConfigProperty(int argc, char **argv);
  void initGridCells() {
    // TODO move initialization to configuration class
    gridCellsX = static_cast<uint>((xmax - xmin) / h) + 1;
    gridCellsY = static_cast<uint>((ymax - ymin) / h) + 1;
    gridCellsZ = static_cast<uint>((zmax - zmin) / h) + 1;
    gridCellCount = gridCellsX * gridCellsY * gridCellsZ;
  }
  ~owConfigProperty() {
    if (simulation != nullptr)
      delete simulation;
  }
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  float zmin;
  float zmax;
  uint gridCellsX;
  uint gridCellsY;
  uint gridCellsZ;
  uint gridCellCount;
  uint numOfElasticP;
  uint numOfLiquidP;
  uint numOfBoundaryP;
  uint numOfMembranes;
  uint MUSCLE_COUNT;

private:
  /** Calculating delta parameter.
   *
   *  "In these situations,
   *  the SPH equations result in falsified values. To circumvent that problem,
   * we pre- compute a single scaling factor δ according to the following
   * formula [1, eq. 8] which is evaluated for a prototype particle with a
   * filled neighborhood. The resulting value is then used for all particles.
   * Finally, we end up with the following equations which are used in the
   * PCISPH method" [1]. [1]
   * http://www.ifi.uzh.ch/vmml/publications/pcisph/pcisph.pdf
   */
  inline void calcDelta(double gradWspikyCoefficient) {

    float x[] = {1, 1, 0, -1, -1, -1, 0, 1, 1, 1,  0, -1, -1, -1, 0, 1,
                 1, 1, 0, -1, -1, -1, 0, 1, 2, -2, 0, 0,  0,  0,  0, 0};
    float y[] = {0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1,  0, -1, -1, -1,
                 0, 1, 1, 1, 0, -1, -1, -1, 0, 0, 2, -2, 0, 0,  0,  0};
    float z[] = {0,  0,  0,  0,  0,  0,  0,  0,  1, 1, 1, 1, 1, 1,  1, 1,
                 -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 2, -2, 1, -1};
    float sum1_x = 0.f;
    float sum1_y = 0.f;
    float sum1_z = 0.f;
    double sum1 = 0.0, sum2 = 0.0;
    float v_x = 0.f;
    float v_y = 0.f;
    float v_z = 0.f;
    float dist;
    float particleRadius =
        pow(constMap["mass"] / constMap["rho0"],
            1.f / 3.f); // It's equal to simulationScale
                        // TODO: replace it with simulation scale
    float h_r_2;

    for (int i = 0; i < MAX_NEIGHBOR_COUNT; i++) {
      v_x = x[i] * 0.8f /*1.f*/ *
            particleRadius; // return it back to 0.8 it's more stable
      v_y = y[i] * 0.8f /*1.f*/ *
            particleRadius; // return it back to 0.8 it's more stable
      v_z = z[i] * 0.8f /*1.f*/ *
            particleRadius; // return it back to 0.8 it's more stable

      dist = sqrt(v_x * v_x + v_y * v_y + v_z * v_z); // scaled, right?

      if (dist <= constMap["h"] * constMap["simulationScale"]) {
        h_r_2 = pow((constMap["h"] * constMap["simulationScale"] - dist),
                    2); // scaled

        sum1_x += h_r_2 * v_x / dist;
        sum1_y += h_r_2 * v_y / dist;
        sum1_z += h_r_2 * v_z / dist;

        sum2 += h_r_2 * h_r_2;
      }
    }
    sum1 = sum1_x * sum1_x + sum1_y * sum1_y + sum1_z * sum1_z;
    double result = 1.0 / (constMap["beta"] * gradWspikyCoefficient *
                           gradWspikyCoefficient * (sum1 + sum2));
    // return  1.0f / (beta * gradWspikyCoefficient * gradWspikyCoefficient *
    // (sum1 + sum2));
    delta = (float)result;
  }
  void fillConstMap();
  int PARTICLE_COUNT;
  int PARTICLE_COUNT_RoundedUp;
  int totalNumberOfIterations;
  int logStep;
  float timeStep;
  float timeLim;
  float beta;
  float delta;
  DEVICE prefDeviceType;        // 0-CPU, 1-GPU
  INTEGRATOR integrationMethod; // DEFAULT is EULER
  std::string configFileName;
  std::string path;     // PATH to configuration files
  std::string loadPath; // PATH to load buffer files
  // SignalSimulator simulation;
  owINeuronSimulator *simulation;
  std::string devFullName;
  std::string sourceFileName;
  bool nrnSimRun; // indicates if we also run NEURON simulation
  bool c302;      // indicates if we also run NEURON simulation
  std::string nrnSimulationFileName;
  std::map<std::string, float> constMap;
};

#endif /* OWCONFIGURATION_H_ */
