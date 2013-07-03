#include "owPhysicsFluidSimulator.h"
#include <stdexcept>
#include <iostream>
#include <fstream>

float calcDelta();
extern const float delta = calcDelta();
int iterationCount = 0;
/*extern*/ int numOfElasticConnections = 0;
/*extern*/ int numOfLiquidP = 0;
/*extern*/ int numOfElasticP = 0;
/*extern*/ int numOfBoundaryP = 0;
int numOfMembranes = 0;
int * _particleIndex;
unsigned int * gridNextNonEmptyCellBuffer;
extern int gridCellCount;
extern float * muscle_activation_signal_cpp;

owPhysicsFluidSimulator::owPhysicsFluidSimulator(owHelper * helper)
{
	//int generateInitialConfiguration = 1;//1 to generate initial configuration, 0 - load from file

	try{
		if(generateInitialConfiguration)
		// GENERATE THE SCENE
		owHelper::generateConfiguration(0, position_cpp, velocity_cpp, elasticConnectionsData_cpp, membraneData_cpp, numOfLiquidP, numOfElasticP, numOfBoundaryP, numOfElasticConnections, numOfMembranes, particleMembranesList_cpp);	
		else								
		// LOAD FROM FILE
		owHelper::preLoadConfiguration();	
											//=======================

		position_cpp = new float[ 8 * PARTICLE_COUNT ];
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
		
		if(generateInitialConfiguration)	
		// GENERATE THE SCENE
		owHelper::generateConfiguration(1,position_cpp, velocity_cpp, elasticConnectionsData_cpp, membraneData_cpp, numOfLiquidP, numOfElasticP, numOfBoundaryP, numOfElasticConnections, numOfMembranes, particleMembranesList_cpp );	
		else 
		// LOAD FROM FILE	
		owHelper::loadConfiguration( position_cpp, velocity_cpp, elasticConnectionsData_cpp, numOfLiquidP, numOfElasticP, numOfBoundaryP, numOfElasticConnections );		//Load configuration from file to buffer
											
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

double owPhysicsFluidSimulator::simulationStep()
{
	//PCISPH algorithm
	int iter = 0;
	helper->refreshTime();
	printf("\n[[ Step %d ]]\n",iterationCount);
	try{
		//SEARCH FOR NEIGHBOURS PART
		ocl_solver->_runClearBuffers();								helper->watch_report("_runClearBuffers: \t%9.3f ms\n");
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
			ocl_solver->_run_pcisph_predictPositions();				
			ocl_solver->_run_pcisph_predictDensity();				
			ocl_solver->_run_pcisph_correctPressure();				
			ocl_solver->_run_pcisph_computePressureForceAcceleration();
			iter++;
		}while( iter < maxIteration );

		ocl_solver->_run_clearMembraneBuffers();
		ocl_solver->_run_computeInteractionWithMembranes();
		ocl_solver->_run_pcisph_integrate();						helper->watch_report("_runPCISPH: \t\t%9.3f ms\t3 iteration(s)\n");
		ocl_solver->read_position_buffer(position_cpp);				helper->watch_report("_readBuffer: \t\t%9.3f ms\n"); 
		//END PCISPH algorithm
		printf("------------------------------------\n");
		printf("_Total_step_time:\t%9.3f ms\n",helper->get_elapsedTime());
		printf("------------------------------------\n");
		iterationCount++;
		//for(int i=0;i<MUSCLE_COUNT;i++) { muscle_activation_signal_cpp[i] *= 0.9f; }
		ocl_solver->updateMuscleActivityData(muscle_activation_signal_cpp);
		return helper->get_elapsedTime();
	}catch(std::exception &e){
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