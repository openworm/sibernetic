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

#include <stdexcept>
#include <iostream>
#include <fstream>

#include "PyramidalSimulation.h"
#include "owPhysicsFluidSimulator.h"

float calcDelta();
extern const float delta = calcDelta();
int iterationCount = 0;
int numOfElasticConnections = 0;
int numOfLiquidP = 0;
int numOfElasticP = 0;
int numOfBoundaryP = 0;
int numOfMembranes = 0;
int * _particleIndex;
unsigned int * gridNextNonEmptyCellBuffer;
extern int gridCellCount;
extern float * muscle_activation_signal_cpp;
int iter_step = 10;

//mv
//need to find a more elegant design for this - at the moment the use of a global
//is pretty ugly:
#ifdef PY_NETWORK_SIMULATION
PyramidalSimulation simulation;
#endif
owPhysicsFluidSimulator::owPhysicsFluidSimulator(owHelper * helper)
{
	//int generateInitialConfiguration = 1;//1 to generate initial configuration, 0 - load from file

	try{
		if(generateWormBodyConfiguration)
		// GENERATE THE SCENE
		owHelper::generateConfiguration(0, position_cpp, velocity_cpp, elasticConnectionsData_cpp, membraneData_cpp, numOfLiquidP, numOfElasticP, numOfBoundaryP, numOfElasticConnections, numOfMembranes, particleMembranesList_cpp);	
		else								
		// LOAD FROM FILE
		owHelper::preLoadConfiguration(numOfMembranes);
#ifdef PY_NETWORK_SIMULATION
        //mv
		simulation.setup();
#endif
											//=======================

		position_cpp = new float[ 4 * PARTICLE_COUNT ];
		velocity_cpp = new float[ 4 * PARTICLE_COUNT ];
		_particleIndex = new   int[ 2 * PARTICLE_COUNT ];
		gridNextNonEmptyCellBuffer = new unsigned int[gridCellCount+1];
		muscle_activation_signal_cpp = new float [MUSCLE_COUNT];
		if(numOfMembranes<=0) membraneData_cpp = NULL; else membraneData_cpp = new int [numOfMembranes*3];
		if(numOfElasticP<=0) particleMembranesList_cpp = NULL; 
			else particleMembranesList_cpp = new int [numOfElasticP*MAX_MEMBRANES_INCLUDING_SAME_PARTICLE];

		for(int i=0;i<MUSCLE_COUNT;i++)
		{
			muscle_activation_signal_cpp[i] = 0.f;
		}

		//The buffers listed below are only for usability and debug
		density_cpp = new float[ 1 * PARTICLE_COUNT ];
		particleIndex_cpp = new unsigned int[PARTICLE_COUNT * 2];
		acceleration_cpp = new float[PARTICLE_COUNT * 4];//TODO REMOVE IT AFTER FIXING
		
		if(generateWormBodyConfiguration)
		// GENERATE THE SCENE
		owHelper::generateConfiguration(1,position_cpp, velocity_cpp, elasticConnectionsData_cpp, membraneData_cpp, numOfLiquidP, numOfElasticP, numOfBoundaryP, numOfElasticConnections, numOfMembranes, particleMembranesList_cpp );	
		else 
		// LOAD FROM FILE	
		owHelper::loadConfiguration( position_cpp, velocity_cpp, elasticConnectionsData_cpp, numOfLiquidP, numOfElasticP, numOfBoundaryP, numOfElasticConnections, numOfMembranes,membraneData_cpp, particleMembranesList_cpp );		//Load configuration from file to buffer
											
		if(numOfElasticP != 0){
			ocl_solver = new owOpenCLSolver(position_cpp, velocity_cpp, elasticConnectionsData_cpp, membraneData_cpp, particleMembranesList_cpp);	//Create new openCLsolver instance
		}else
			ocl_solver = new owOpenCLSolver(position_cpp,velocity_cpp);	//Create new openCLsolver instance
		this->helper = helper;
	}catch( std::exception &e ){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}

double owPhysicsFluidSimulator::simulationStep(const bool load_to)
{
	//PCISPH algorithm
	int iter = 0;//PCISPH prediction-correction iterations conter
	//if(iterationCount!=0) return 0.0;//uncomment this line to stop movement of the scene
	helper->refreshTime();
	printf("\n[[ Step %d ]]\n",iterationCount);
	try{
		//SEARCH FOR NEIGHBOURS PART
//		ocl_solver->_runClearBuffers();								helper->watch_report("_runClearBuffers: \t%9.3f ms\n");
		ocl_solver->_runHashParticles();							helper->watch_report("_runHashParticles: \t%9.3f ms\n");
		ocl_solver->_runSort();										helper->watch_report("_runSort: \t\t%9.3f ms\n");
		ocl_solver->_runSortPostPass();								helper->watch_report("_runSortPostPass: \t%9.3f ms\n");
		ocl_solver->_runIndexx();									helper->watch_report("_runIndexx: \t\t%9.3f ms\n");
		ocl_solver->_runIndexPostPass();							helper->watch_report("_runIndexPostPass: \t%9.3f ms\n");
		ocl_solver->_runFindNeighbors();							helper->watch_report("_runFindNeighbors: \t%9.3f ms\n");
		//PCISPH PART
		ocl_solver->_run_pcisph_computeDensity();					
		ocl_solver->_run_pcisph_computeForcesAndInitPressure();		
		ocl_solver->_run_pcisph_computeElasticForces();
		do{
			//printf("\n^^^^ iter %d ^^^^\n",iter);
			ocl_solver->_run_pcisph_predictPositions();				
			ocl_solver->_run_pcisph_predictDensity();				
			ocl_solver->_run_pcisph_correctPressure();				
			ocl_solver->_run_pcisph_computePressureForceAcceleration();
			iter++;
		}while( iter < maxIteration );

		ocl_solver->_run_pcisph_integrate(iterationCount);			helper->watch_report("_runPCISPH: \t\t%9.3f ms\t3 iteration(s)\n");
		//Handling of Interaction with membranes
		ocl_solver->_run_clearMembraneBuffers();
		ocl_solver->_run_computeInteractionWithMembranes();
		// compute change of coordinates due to interactions with membranes
		ocl_solver->_run_computeInteractionWithMembranes_finalize();
		//END
		ocl_solver->read_position_buffer(position_cpp);				helper->watch_report("_readBuffer: \t\t%9.3f ms\n"); 

		//END PCISPH algorithm
		printf("------------------------------------\n");
		printf("_Total_step_time:\t%9.3f ms\n",helper->get_elapsedTime());
		printf("------------------------------------\n");
		if(load_to){
			if(iterationCount == 0){
				owHelper::loadConfigurationToFile(position_cpp,elasticConnectionsData_cpp,membraneData_cpp);
			}else{
				if(iterationCount % iter_step == 0){
					owHelper::loadConfigurationToFile(position_cpp, NULL, NULL, false);
				}
			}
		}
		iterationCount++;
		//for(int i=0;i<MUSCLE_COUNT;i++) { muscle_activation_signal_cpp[i] *= 0.9f; }
#ifdef PY_NETWORK_SIMULATION
        //mv
        vector<float> muscle_vector = simulation.run();
        for(int i=0; i<MUSCLE_COUNT; i++){
        	for (long index = 0; index < muscle_vector.size(); index++){
        		muscle_activation_signal_cpp[index] = muscle_vector[index];
        	}
        }
#endif
		ocl_solver->updateMuscleActivityData(muscle_activation_signal_cpp);
		return helper->get_elapsedTime();
	}
	catch(std::exception &e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}

//Destructor
owPhysicsFluidSimulator::~owPhysicsFluidSimulator(void)
{
	delete [] position_cpp;
	delete [] velocity_cpp;
	delete [] density_cpp;
	delete [] particleIndex_cpp;
	delete [] muscle_activation_signal_cpp;
	if(membraneData_cpp) delete [] membraneData_cpp;
	//if(particleMembranesList_cpp) delete [] particleMembranesList_cpp;
	ocl_solver->~owOpenCLSolver();
}

float calcDelta()
{
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
	float particleRadius = pow(mass/rho0,1.f/3.f);  // the value is about 0.01 instead of 
	float h_r_2;									// my previous estimate = simulationScale*h/2 = 0.0066

    for (int i = 0; i < 32; i++)
    {
		v_x = x[i] * 0.8f * particleRadius;
		v_y = y[i] * 0.8f * particleRadius;
		v_z = z[i] * 0.8f * particleRadius;

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
	return (float)result;
}
