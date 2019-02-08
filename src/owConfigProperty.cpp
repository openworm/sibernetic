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
#include "owConfigProperty.h"
#include "owHelper.h"
#include <stdio.h>

owConfigProperty::owConfigProperty(int argc, char **argv)
    : numOfElasticP(0), numOfLiquidP(0), numOfBoundaryP(0), numOfMembranes(0),
      MUSCLE_COUNT(96), logStep(10), path("./configuration/"),
      loadPath("./buffers/"), sourceFileName(OPENCL_PROGRAM_PATH) {
  prefDeviceType = ALL;
  this->timeStep = ::timeStep;
  timeLim = 0.f;
  beta = ::beta;
  integrationMethod = EULER;
  std::string strTemp;
  configFileName = "demo1"; // by default
  std::string simName = "";
  nrnSimRun = false;
  nrnSimulationFileName = "";
  simulation = nullptr;
    
  fillConstMap(); // map must be filled before parsing arguments, otherwise beta will be NaN because of division by zero

  for (int i = 1; i < argc; i++) {
    strTemp = argv[i];
    if (strTemp.find("device=") == 0) {
      std::transform(strTemp.begin(), strTemp.end(), strTemp.begin(),
                     ::tolower);
      if (strTemp.find("gpu") != std::string::npos)
        prefDeviceType = GPU;
      if (strTemp.find("cpu") != std::string::npos)
        prefDeviceType = CPU;
    }
    if (strTemp.find("timestep=") == 0) {

      this->timeStep = ::atof(strTemp.substr(strTemp.find('=') + 1).c_str());
      if (this->timeStep < 0.0f)
        std::cout << "timeStep < 0 using default value " << ::timeStep
                  << std::endl;
      this->timeStep = (this->timeStep > 0) ? this->timeStep : ::timeStep;

      // also we shoisSimulationRun = true;uld recalculate beta if time_step is
      // different from default value of timeStep in owPhysicsConstant
      beta = this->timeStep * this->timeStep * constMap["mass"] *
             constMap["mass"] * 2 / (constMap["rho0"] * constMap["rho0"]);
    }

    if (strTemp.find("timelimit=") == 0) {
      timeLim = ::atof(strTemp.substr(strTemp.find('=') + 1).c_str());
      if (timeLim < 0.0)
        throw std::runtime_error(
            "timelimit could not be less than 0 check input parameters");
    }
    if (strTemp.find("LEAPFROG") != std::string::npos ||
        strTemp.find("leapfrog") != std::string::npos) {
      integrationMethod = LEAPFROG;
    }
    if (strTemp.find("logstep=") == 0) {
      logStep = ::atoi(strTemp.substr(strTemp.find('=') + 1).c_str());
      if (logStep < 1)
        throw std::runtime_error(
            "logStep could not be less than 1 check input parameters");
    }
    if (strTemp.find("lpath=") != std::string::npos) {
      loadPath = strTemp.substr(strTemp.find('=') + 1).c_str();
    }
    if (strTemp.find("oclsourcepath=") != std::string::npos) {
      sourceFileName = strTemp.substr(strTemp.find('=') + 1).c_str();
    }
    if (strTemp == "-f") {
      if (i + 1 < argc) {
        configFileName = argv[i + 1];
        if (configFileName.find("\\") != std::string::npos ||
            configFileName.find("/") != std::string::npos) {
          std::size_t found = configFileName.find_last_of("/\\");
          path = configFileName.substr(0, found + 1);
          configFileName = configFileName.substr(found + 1);
        }
      } else
        throw std::runtime_error("You forget add configuration file name. "
                                 "Please add it and try again");
    }
    if (strTemp.find("sigsim=") == 0) {
      simName = strTemp.substr(strTemp.find('=') + 1).c_str();
    }
    if (strTemp.find("-nrn") == 0) {
      nrnSimRun = true;
      // Next step load custom model file
      if (i + 1 < argc) {
        nrnSimulationFileName = argv[i + 1];
      } else
        throw std::runtime_error("You forget add NEURON model file name. "
                                 "Please add it and try again");
    }

    if (strTemp.find("-c302") == 0) {
      // int result = strTemp.find("-c302");
      // printf("\n%s: %d c302=%d\n",strTemp.c_str(),result,c302==true);
      // if(result == 0){
      c302 = true;
    }
  }
  owHelper::preLoadConfiguration(this);

  totalNumberOfIterations = timeLim / this->timeStep; // if it equals to 0 it
                                                      // means that simulation
                                                      // will work infinitely
  if (isWormConfig() ||
      nrnSimRun) { // in case if we run worm configuration TODO make it optional

    if (isWormConfig() && !nrnSimRun) {
      std::string pythonClass = "MuscleSimulation";
      // printf("\nc302=%d\n",(int)isC302());
      if (isC302() == true) {
        pythonClass = "C302NRNSimulation";
      }

      if (simName.compare("") == 0)
        simulation = new SignalSimulator("main_sim", pythonClass, this->timeStep);
      else
        simulation = new SignalSimulator(simName, pythonClass, this->timeStep);
    } else {
      simulation =
          new owNeuronSimulator(1, this->timeStep, nrnSimulationFileName);
    }
  }
}
void showAllConst(std::map<std::string, float> &consts);
void owConfigProperty::fillConstMap() {
  if (constMap.find("mass") == constMap.end())
    constMap["mass"] = mass;

  if (constMap.find("h") == constMap.end())
    constMap["h"] = h;

  if (constMap.find("rho0") == constMap.end())
    constMap["rho0"] = rho0;

  if (constMap.find("timeStep") == constMap.end())
    constMap["timeStep"] = timeStep;
  else
    timeStep = constMap["timeStep"];

  if (constMap.find("simulationScale") == constMap.end())
    constMap["simulationScale"] = simulationScale;

  if (constMap.find("hashGridCellSize") == constMap.end())
    constMap["hashGridCellSize"] = 2.0f * constMap["h"];

  if (constMap.find("r0") == constMap.end())
    constMap["r0"] = constMap["h"] * 0.5f;

  if (constMap.find("viscosity") == constMap.end())
    constMap["viscosity"] = viscosity;

  if (constMap.find("beta") == constMap.end())
    constMap["beta"] = beta;

  if (constMap.find("gravity_x") == constMap.end())
    constMap["gravity_x"] = gravity_x;

  if (constMap.find("gravity_y") == constMap.end())
    constMap["gravity_y"] = gravity_y;

  if (constMap.find("gravity_z") == constMap.end())
    constMap["gravity_z"] = gravity_z;

  if (constMap.find("maxIteration") == constMap.end())
    constMap["maxIteration"] = maxIteration;

  if (constMap.find("hashGridCellSizeInv") == constMap.end())
    constMap["hashGridCellSizeInv"] = 1.0f / constMap["hashGridCellSize"];

  if (constMap.find("simulationScaleInv") == constMap.end())
    constMap["simulationScaleInv"] = 1.0f / constMap["simulationScale"];

  if (constMap.find("_hScaled") == constMap.end())
    constMap["_hScaled"] = constMap["h"] * constMap["simulationScale"];

  if (constMap.find("_hScaled2") == constMap.end())
    constMap["_hScaled2"] = constMap["_hScaled"] * constMap["_hScaled"];

  if (constMap.find("simulationScaleInv") == constMap.end())
    constMap["simulationScaleInv"] = 1.0f / constMap["simulationScale"];

  if (constMap.find("surfTensCoeff") == constMap.end())
    constMap["surfTensCoeff"] = surfTensCoeff;

  if (constMap.find("elasticityCoefficient") == constMap.end())
    constMap["elasticityCoefficient"] = elasticityCoefficient;

  if (constMap.find("max_muscle_force") == constMap.end())
    constMap["max_muscle_force"] = max_muscle_force;

  double Wpoly6Coefficient =
      315.0 / (64.0 * M_PI *
               pow((double)(constMap["h"] * constMap["simulationScale"]),
                   9.0)); // Wpoly6Coefficient for kernel Wpoly6 [1]
                          // [1] Solenthaler (Dissertation) page 17 eq. (2.20)

  double gradWspikyCoefficient =
      -45.0 /
      (M_PI * pow((double)(constMap["h"] * constMap["simulationScale"]),
                  6.0)); // gradWspikyCoefficient for kernel gradWspiky [1]
                         // [1] Solenthaler (Dissertation) page 18 eq. (2.21)

  double divgradWviscosityCoefficient = -gradWspikyCoefficient;

  constMap["mass_mult_Wpoly6Coefficient"] =
      (float)((double)constMap["mass"] * Wpoly6Coefficient);

  constMap["mass_mult_gradWspikyCoefficient"] =
      (float)((double)constMap["mass"] * gradWspikyCoefficient);

  constMap["mass_mult_divgradWviscosityCoefficient"] =
      (float)((double)constMap["mass"] * divgradWviscosityCoefficient);

  constMap["beta"] = this->timeStep * this->timeStep * constMap["mass"] *
                     constMap["mass"] * 2 /
                     (constMap["rho0"] * constMap["rho0"]);
  calcDelta(gradWspikyCoefficient);
  // showAllConst(constMap);
  // exit(0);
}

void showAllConst(std::map<std::string, float> &consts) {
  std::cout << "rho0: " << consts["rho0"] << std::endl;
  std::cout << "mass: " << consts["mass"] << std::endl;
  std::cout << "timeStep: " << consts["timeStep"] << std::endl;
  std::cout << "simulationScale: " << consts["simulationScale"] << std::endl;
  std::cout << "h: " << consts["h"] << std::endl;
  std::cout << "hashGridCellSize: " << consts["hashGridCellSize"] << std::endl;
  std::cout << "r0: " << consts["r0"] << std::endl;
  std::cout << "viscosity: " << consts["viscosity"] << std::endl;
  std::cout << "beta: " << consts["beta"] << std::endl;
  std::cout << "gravity_x: " << consts["gravity_x"] << std::endl;
  std::cout << "gravity_y: " << consts["gravity_y"] << std::endl;
  std::cout << "gravity_z: " << consts["gravity_z"] << std::endl;
  std::cout << "maxIteration: " << consts["maxIteration"] << std::endl;
  std::cout << "mass_mult_Wpoly6Coefficient: "
            << consts["mass_mult_Wpoly6Coefficient"] << std::endl;
  std::cout << "mass_mult_gradWspikyCoefficient: "
            << consts["mass_mult_gradWspikyCoefficient"] << std::endl;
  std::cout << "mass_mult_divgradWviscosityCoefficient: "
            << consts["mass_mult_divgradWviscosityCoefficient"] << std::endl;
  std::cout << "hashGridCellSizeInv: " << consts["hashGridCellSizeInv"]
            << std::endl;
  std::cout << "simulationScaleInv: " << consts["simulationScaleInv"]
            << std::endl;
  std::cout << "_hScaled: " << consts["_hScaled"] << std::endl;
  std::cout << "_hScaled2: " << consts["_hScaled2"] << std::endl;
  std::cout << "surfTensCoeff: " << consts["surfTensCoeff"] << std::endl;
  std::cout << "elasticityCoefficient: " << consts["elasticityCoefficient"]
            << std::endl;
}
