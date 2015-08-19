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

#include <vector>
#include <string>

#include "owOpenCLConstant.h"
#include "owPhysicsConstant.h"
#include "PyramidalSimulation.h"


struct owConfigProrerty{
	//This value defines boundary of box in which simulation is
	//Sizes of the box containing simulated 'world'
	//Sizes choice is realized this way because it should be proportional to smoothing radius h
public:
	const int getParticleCount(){ return PARTICLE_COUNT; }
	void setParticleCount(int value){
		PARTICLE_COUNT = value;
		PARTICLE_COUNT_RoundedUp = ((( PARTICLE_COUNT - 1 ) / local_NDRange_size ) + 1 ) * local_NDRange_size;
	}
	void setDeviceType(DEVICE type) { preferable_device_type = type; }
	const int getParticleCount_RoundUp(){ return PARTICLE_COUNT_RoundedUp; }
	const int getDeviceType() const { return preferable_device_type; };
	const int getNumberOfIteration() const { return totalNumberOfIteration ;}
	const char * getDeviceName() const { return device_full_name.c_str(); }
	void setDeviceName(const char * name) {
		device_full_name = name;
	}
	INTEGRATOR getIntegrationMethod() const { return integration_method; }
	const std::string & getCofigFileName() const { return configFileName; }
	const std::string & getCofigPath() const { return path; }
	PyramidalSimulation & getPyramidalSimulation() { return simulation; }
	void updatePyramidalSimulation(float * muscleActivationSignal){
		if(isWormConfig()){
			std::vector<float> muscle_vector = simulation.run();
			for (unsigned int index = 0; index < muscle_vector.size(); index++){
				muscleActivationSignal[index] = muscle_vector[index];
			}
		}
	}
	bool isWormConfig(){ return (configFileName == "worm" || configFileName == "worm_no_water")? true:false; }
	void setCofigFileName( const char * name ) { configFileName = name; }
	// Constructor
	owConfigProrerty(int argc, char** argv):numOfElasticP(0), numOfLiquidP(0), numOfBoundaryP(0), numOfMembranes(0), MUSCLE_COUNT(100), path("./configuration/"){
		preferable_device_type = ALL;
		time_step = timeStep;
		time_limit = 0.f;
		beta = ::beta;
		integration_method = EULER;
		std::string s_temp;
		configFileName = "demo1";
		for(int i = 1; i<argc; i++){
			s_temp = argv[i];
			if(s_temp.find("device=") == 0){
				if(s_temp.find("GPU") != std::string::npos || s_temp.find("gpu") != std::string::npos)
					preferable_device_type = GPU;
				if(s_temp.find("CPU") != std::string::npos || s_temp.find("cpu") != std::string::npos)
					preferable_device_type = CPU;
			}
			if(s_temp.find("timestep=") == 0){

				time_step = ::atof( s_temp.substr(s_temp.find('=')+1).c_str());
				time_step = (time_step > 0) ? time_step : timeStep;
				//also we shoisSimulationRun = true;uld recalculate beta if time_step is different from default value of timeStep in owPhysicsConstant
				beta = time_step*time_step*mass*mass*2/(rho0*rho0);
			}
			if(s_temp.find("timelimit=") == 0){
				time_limit = ::atof( s_temp.substr(s_temp.find('=')+1).c_str());
			}
			if(s_temp.find("LEAPFROG") != std::string::npos || s_temp.find("leapfrog") != std::string::npos){
				integration_method = LEAPFROG;
			}
			if(s_temp == "-f"){
				if(i + 1 < argc){
					configFileName = argv[i+1];
					if(configFileName.find("\\") != std::string::npos || configFileName.find("/") != std::string::npos) // If it find \ or / it means that
																														// you indicate full path to configuration in other case it will try to find configuration in configurations folder
						path = "";
				}
				else
					throw std::runtime_error("You forget add configuration file name. Please add it and try again");
			}
		}
		totalNumberOfIteration = time_limit/time_step; // if it equals to 0 it means that simulation will work infinitely
		calcDelta();
		if(isWormConfig()){ // in case if we run worm configuration TODO make it optional
			try{
				simulation.setup();
			}catch(const char * ex){
				std::cerr << ex << std::endl;// TODO make more careful exception catching
				exit(-1);
			}
		}
	};
	float getTimeStep() const { return timeStep; };
	float getDelta() const { return delta; };
	float xmin;
	float xmax;
	float ymin;
	float ymax;
	float zmin;
	float zmax;
	int gridCellsX;
	int gridCellsY;
	int gridCellsZ;
	int gridCellCount;
	int numOfElasticP;
	int numOfLiquidP;
	int numOfBoundaryP;
	int numOfMembranes;
	int MUSCLE_COUNT;
private:
	/** Calculating delta parameter.
	 *
	 *  "In these situations,
	 *  the SPH equations result in falsified values. To circumvent that problem, we pre-
	 *  compute a single scaling factor δ according to the following formula [1, eq. 8] which is
	 *  evaluated for a prototype particle with a filled neighborhood. The resulting value
	 *  is then used for all particles. Finally, we end up with the following equations
	 *  which are used in the PCISPH method" [1].
	 *  [1] http://www.ifi.uzh.ch/vmml/publications/pcisph/pcisph.pdf
	 */
	inline void calcDelta(){

		float x[] = { 1, 1, 0,-1,-1,-1, 0, 1, 1, 1, 0,-1,-1,-1, 0, 1, 1, 1, 0,-1,-1,-1, 0, 1, 2,-2, 0, 0, 0, 0, 0, 0 };
	    float y[] = { 0, 1, 1, 1, 0,-1,-1,-1, 0, 1, 1, 1, 0,-1,-1,-1, 0, 1, 1, 1, 0,-1,-1,-1, 0, 0, 2,-2, 0, 0, 0, 0 };
	    float z[] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 2,-2, 1,-1 };
	    float sum1_x = 0.f;
		float sum1_y = 0.f;
		float sum1_z = 0.f;
	    double sum1 = 0.0, sum2 = 0.0;
		float v_x = 0.f;
		float v_y = 0.f;
		float v_z = 0.f;
		float dist;
		float particleRadius = pow(mass/rho0,1.f/3.f);  // It's equal to simulationScale
														// TODO: replace it with simulation scale
		float h_r_2;

	    for (int i = 0; i < MAX_NEIGHBOR_COUNT; i++)
	    {
			v_x = x[i] * 0.8f/*1.f*/ * particleRadius; // return it back to 0.8 it's more stable
			v_y = y[i] * 0.8f/*1.f*/ * particleRadius; // return it back to 0.8 it's more stable
			v_z = z[i] * 0.8f/*1.f*/ * particleRadius; // return it back to 0.8 it's more stable

	        dist = sqrt(v_x*v_x+v_y*v_y+v_z*v_z);//scaled, right?

	        if (dist <= h*simulationScale)
	        {
				h_r_2 = pow((h*simulationScale - dist),2);//scaled

	            sum1_x += h_r_2 * v_x / dist;
				sum1_y += h_r_2 * v_y / dist;
				sum1_z += h_r_2 * v_z / dist;

	            sum2 += h_r_2 * h_r_2;
	        }
	    }
		sum1 = sum1_x*sum1_x + sum1_y*sum1_y + sum1_z*sum1_z;
		double result = 1.0 / (beta * gradWspikyCoefficient * gradWspikyCoefficient * (sum1 + sum2));
		//return  1.0f / (beta * gradWspikyCoefficient * gradWspikyCoefficient * (sum1 + sum2));
		delta = (float)result;
	}

	int PARTICLE_COUNT;
	int PARTICLE_COUNT_RoundedUp;
	int totalNumberOfIteration;
	float time_step;
	float time_limit;
	float beta;
	float delta;
	DEVICE preferable_device_type;// 0-CPU, 1-GPU
	INTEGRATOR integration_method; //DEFAULT is EULER
	std::string configFileName;
	std::string path; // PATH to configuration files
	PyramidalSimulation simulation;
	std::string device_full_name;
};

#endif /* OWCONFIGURATION_H_ */
