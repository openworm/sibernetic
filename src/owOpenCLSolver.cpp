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

#include "owOpenCLSolver.h"

extern int numOfElasticP;
extern int numOfBoundaryP;
extern int numOfMembranes;
extern int * _particleIndex;
extern unsigned int * gridNextNonEmptyCellBuffer;

int myCompare( const void * v1, const void * v2 ); 

owOpenCLSolver::owOpenCLSolver(const float * position_cpp, const float * velocity_cpp, owConfigProrerty * config, const float * elasticConnectionsData_cpp, const int * membraneData_cpp, const int * particleMembranesList_cpp)
{
	try{
		initializeOpenCL();
		// Create OpenCL buffers
		create_ocl_buffer( "acceleration", acceleration, CL_MEM_READ_WRITE, ( config->getParticleCount() * sizeof( float ) * 4 * 3 ) );
		create_ocl_buffer( "gridCellIndex", gridCellIndex, CL_MEM_READ_WRITE, ( ( config->gridCellCount + 1 ) * sizeof( unsigned int ) * 1 ) );
		create_ocl_buffer( "gridCellIndexFixedUp", gridCellIndexFixedUp, CL_MEM_READ_WRITE, ( ( config->gridCellCount + 1 ) * sizeof( unsigned int ) * 1 ) );
		create_ocl_buffer( "neighborMap", neighborMap, CL_MEM_READ_WRITE, ( MAX_NEIGHBOR_COUNT * config->getParticleCount() * sizeof( float ) * 2 ) );
		create_ocl_buffer( "particleIndex", particleIndex, CL_MEM_READ_WRITE, ( config->getParticleCount() * sizeof( unsigned int ) * 2 ) );
		create_ocl_buffer( "particleIndexBack", particleIndexBack, CL_MEM_READ_WRITE, ( config->getParticleCount() * sizeof( unsigned int ) ) );
		create_ocl_buffer( "position", position, CL_MEM_READ_WRITE, ( config->getParticleCount() * sizeof( float ) * 4 * (1 + 1/*1 extra, for membrane handling*/)) );
		create_ocl_buffer( "pressure", pressure, CL_MEM_READ_WRITE, ( config->getParticleCount() * sizeof( float ) * 1 ) );
		create_ocl_buffer( "rho", rho, CL_MEM_READ_WRITE, ( config->getParticleCount() * sizeof( float ) * 2 ) );
		create_ocl_buffer( "sortedPosition", sortedPosition, CL_MEM_READ_WRITE, ( config->getParticleCount() * sizeof( float ) * 4 * 2 ) );
		create_ocl_buffer( "sortedVelocity", sortedVelocity, CL_MEM_READ_WRITE, ( config->getParticleCount() * sizeof( float ) * 4 ) );
		create_ocl_buffer( "velocity", velocity, CL_MEM_READ_WRITE, ( config->getParticleCount() * sizeof( float ) * 4 * (1 + 1/*1 extra, for membrane handling*/) ) );
		create_ocl_buffer( "muscle_activation_signal", muscle_activation_signal, CL_MEM_READ_WRITE, ( MUSCLE_COUNT * sizeof( float ) ) );
		// Create OpenCL kernels
		create_ocl_kernel("clearBuffers", clearBuffers);
		create_ocl_kernel("findNeighbors", findNeighbors);
		create_ocl_kernel("hashParticles", hashParticles);
		create_ocl_kernel("indexx", indexx);
		create_ocl_kernel("sortPostPass", sortPostPass);
		// Additional PCISPH-related kernels 
		create_ocl_kernel("pcisph_computeForcesAndInitPressure", pcisph_computeForcesAndInitPressure);
		create_ocl_kernel("pcisph_integrate", pcisph_integrate);
		create_ocl_kernel("pcisph_predictPositions", pcisph_predictPositions);
		create_ocl_kernel("pcisph_predictDensity", pcisph_predictDensity);
		create_ocl_kernel("pcisph_correctPressure", pcisph_correctPressure);
		create_ocl_kernel("pcisph_computePressureForceAcceleration", pcisph_computePressureForceAcceleration);
		create_ocl_kernel("pcisph_computeDensity", pcisph_computeDensity);
		create_ocl_kernel("pcisph_computeElasticForces", pcisph_computeElasticForces);
		// membrane handling kernels
		create_ocl_kernel("clearMembraneBuffers",clearMembraneBuffers);
		create_ocl_kernel("computeInteractionWithMembranes",computeInteractionWithMembranes);
		create_ocl_kernel("computeInteractionWithMembranes_finalize",computeInteractionWithMembranes_finalize);
		//Copy position_cpp and velocity_cpp to the OpenCL Device
		copy_buffer_to_device( position_cpp, position, config->getParticleCount() * sizeof( float ) * 4 );
		copy_buffer_to_device( velocity_cpp, velocity, config->getParticleCount() * sizeof( float ) * 4 );
		//membranes
		if(membraneData_cpp != NULL )
		{
			create_ocl_buffer( "membraneData", membraneData, CL_MEM_READ_WRITE, ( numOfMembranes * sizeof( int ) * 3 ) );
			copy_buffer_to_device( membraneData_cpp, membraneData, numOfMembranes * sizeof( int ) * 3 );

			if(particleMembranesList_cpp != NULL) //in actual version I'm going to support only membrance built upon elastic matter particles (interconnected with springs -- highly recommended)
			{
				create_ocl_buffer("particleMembranesList", particleMembranesList,CL_MEM_READ_WRITE, numOfElasticP * MAX_MEMBRANES_INCLUDING_SAME_PARTICLE * sizeof(int) );
				copy_buffer_to_device( particleMembranesList_cpp, particleMembranesList, numOfElasticP * MAX_MEMBRANES_INCLUDING_SAME_PARTICLE * sizeof( int ) );
			}

			if(particleMembranesList_cpp) delete [] particleMembranesList_cpp;
		}
		//elastic connections
		if(elasticConnectionsData_cpp != NULL){
			create_ocl_buffer("elasticConnectionsData", elasticConnectionsData,CL_MEM_READ_WRITE, numOfElasticP * MAX_NEIGHBOR_COUNT * sizeof(float) * 4);
			copy_buffer_to_device(elasticConnectionsData_cpp, elasticConnectionsData, numOfElasticP * MAX_NEIGHBOR_COUNT * sizeof(float) * 4);
		}

	}catch( std::exception &e ){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}

extern char device_full_name [1000];
void owOpenCLSolver::initializeOpenCL()
{
	cl_int err;
	std::vector< cl::Platform > platformList;
	err = cl::Platform::get( &platformList );
	if( platformList.size() < 1 ){
		throw std::runtime_error( "No OpenCL platforms found" );
	}
	char cBuffer[1024];
	cl_platform_id cl_pl_id[10];
	cl_uint n_pl;
	clGetPlatformIDs(10,cl_pl_id,&n_pl);
	cl_int ciErrNum;
	int sz;
	for(int i=0;i<(int)n_pl;i++)
	{
		// Get OpenCL platform name and version
		ciErrNum = clGetPlatformInfo (cl_pl_id[i], CL_PLATFORM_VERSION, sz = sizeof(cBuffer), cBuffer, NULL);
		if (ciErrNum == CL_SUCCESS)
		{
			printf(" CL_PLATFORM_VERSION [%d]: \t%s\n", i, cBuffer);
		} 
		else
		{
			printf(" Error %i in clGetPlatformInfo Call !!!\n\n", ciErrNum);
		}
	}
	//0-CPU, 1-GPU // depends on the time order of system OpenCL drivers installation on your local machine
	// CL_DEVICE_TYPE
    cl_device_type type;
	const int device_type [] = {CL_DEVICE_TYPE_CPU,CL_DEVICE_TYPE_GPU};
	int preferable_device_type = 0;// 0-CPU, 1-GPU
	
	unsigned int plList = 0;//selected platform index in platformList array [choose CPU by default]
	//added autodetection of device number corresonding to preferrable device type (CPU|GPU) | otherwise the choice will be made from list of existing devices
	cl_uint ciDeviceCount;
	cl_device_id * devices_t;
	bool bPassed = true, findDevice = false;
	cl_int result;
	for(int clSelectedPlatformID = 0;clSelectedPlatformID < (int)n_pl;clSelectedPlatformID++){
		if(findDevice)
			break;
		clGetDeviceIDs (cl_pl_id[clSelectedPlatformID], CL_DEVICE_TYPE_ALL, 0, NULL, &ciDeviceCount);
		if ((devices_t = (cl_device_id*)malloc(sizeof(cl_device_id) * ciDeviceCount)) == NULL){
		   bPassed = false;
		}
		if(bPassed){
			result= clGetDeviceIDs (cl_pl_id[clSelectedPlatformID], CL_DEVICE_TYPE_ALL, ciDeviceCount, devices_t, &ciDeviceCount);
			if( result == CL_SUCCESS){
				for( cl_uint i =0; i < ciDeviceCount; ++i ){
					clGetDeviceInfo(devices_t[i], CL_DEVICE_TYPE, sizeof(type), &type, NULL);		
					if( type & device_type[preferable_device_type]){
						plList = clSelectedPlatformID;
						findDevice = true;
						break;
					}
				}
			}
		}
	}
	if(!findDevice) plList = 0;
	cl_context_properties cprops[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties) (platformList[plList])(), 0 };
	context = cl::Context( device_type[preferable_device_type], cprops, NULL, NULL, &err );
	devices = context.getInfo< CL_CONTEXT_DEVICES >();
	if( devices.size() < 1 ){
		throw std::runtime_error( "No OpenCL devices found" );
	}
	//Print some information about chosen platform
	int value;
    unsigned long val2;
    size_t val3;
	//uint deviceNum = 0;// causes "error C2065: 'uint' : undeclared identifier"
    unsigned int deviceNum = 0;
	result = devices[deviceNum].getInfo(CL_DEVICE_NAME,&cBuffer);// CL_INVALID_VALUE = -30;
	if(result == CL_SUCCESS) std::cout << "CL_CONTEXT_PLATFORM ["<< plList << "]: CL_DEVICE_NAME [" << deviceNum << "]:\t" << cBuffer << "\n" << std::endl;
	if(strlen(cBuffer)<1000) strcpy(device_full_name,cBuffer);
	result = devices[deviceNum].getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE,&val3);
	if(result == CL_SUCCESS) std::cout << "CL_CONTEXT_PLATFORM ["<< plList << "]: CL_DEVICE_MAX_WORK_GROUP_SIZE [" <<  deviceNum <<"]: \t" << val3 <<std::endl;
	result = devices[deviceNum].getInfo(CL_DEVICE_MAX_COMPUTE_UNITS,&value);
	if(result == CL_SUCCESS) std::cout<<"CL_CONTEXT_PLATFORM [" << plList << "]: CL_DEVICE_MAX_COMPUTE_UNITS [" << deviceNum << "]: \t" << value  << std::endl;
	result = devices[deviceNum].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE,&val2);
	if(result == CL_SUCCESS) std::cout<<"CL_CONTEXT_PLATFORM [" << plList <<"]: CL_DEVICE_GLOBAL_MEM_SIZE ["<< deviceNum <<"]: \t" << deviceNum <<std::endl;
	result = devices[deviceNum].getInfo(CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,&val2);
	if(result == CL_SUCCESS) std::cout << "CL_CONTEXT_PLATFORM [" << plList <<"]: CL_DEVICE_GLOBAL_MEM_CACHE_SIZE [" << deviceNum <<"]:\t" << val2 <<std::endl;
	result = devices[deviceNum].getInfo(CL_DEVICE_LOCAL_MEM_SIZE,&val2);
	if(result == CL_SUCCESS) std::cout << "CL_CONTEXT_PLATFORM " << plList <<": CL_DEVICE_LOCAL_MEM_SIZE ["<< deviceNum <<"]:\t" << val2 << std::endl;
	
	queue = cl::CommandQueue( context, devices[ deviceNum ], 0, &err );
	if( err != CL_SUCCESS ){
		throw std::runtime_error( "failed to create command queue" );
	}

	std::string sourceFileName( OPENCL_PROGRAM_PATH );
	std::ifstream file( sourceFileName.c_str() );
	if( !file.is_open() ){
		throw std::runtime_error( "could not open file " + sourceFileName );
	}
	std::string programSource( std::istreambuf_iterator<char>( file ), ( std::istreambuf_iterator<char>() ));
	cl::Program::Sources source( 1, std::make_pair( programSource.c_str(), programSource.length()+1 ));
	program = cl::Program( context, source );
#if defined(__APPLE__)
	err = program.build( devices, "-g -cl-opt-disable" );
#else
	#if INTEL_OPENCL_DEBUG
		err = program.build( devices, OPENCL_DEBUG_PROGRAM_PATH +  "-g -cl-opt-disable");
	#else
		err = program.build( devices, "-g -cl-opt-disable");
	#endif
#endif
	if( err != CL_SUCCESS ){
		std::string compilationErrors;
		compilationErrors = program.getBuildInfo< CL_PROGRAM_BUILD_LOG >( devices[ 0 ] );
		std::cerr << "Compilation failed: " << std::endl << compilationErrors << std::endl;
		throw std::runtime_error( "failed to build program" );
	}
	return;
}
//Kernels functions definition
unsigned int owOpenCLSolver::_runClearBuffers(owConfigProrerty * config)
{
	// Stage ClearBuffers
	clearBuffers.setArg( 0, neighborMap );
	clearBuffers.setArg( 1, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(clearBuffers, cl::NullRange, cl::NDRange( (int) ( config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_runHashParticles(owConfigProrerty * config)
{
	// Stage HashParticles
	hashParticles.setArg( 0, position );
	hashParticles.setArg( 1, config->gridCellsX );
	hashParticles.setArg( 2, config->gridCellsY );
	hashParticles.setArg( 3, config->gridCellsZ );
	hashParticles.setArg( 4, hashGridCellSizeInv );
	hashParticles.setArg( 5, config->xmin );
	hashParticles.setArg( 6, config->ymin );
	hashParticles.setArg( 7, config->zmin );
	hashParticles.setArg( 8, particleIndex );
	hashParticles.setArg( 9, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(
		hashParticles, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}

unsigned int owOpenCLSolver::_runSort(owConfigProrerty * config)
{
	copy_buffer_from_device( _particleIndex, particleIndex, config->getParticleCount() * 2 * sizeof( int ) );
	qsort( _particleIndex, config->getParticleCount(), 2 * sizeof( int ), myCompare );
	copy_buffer_to_device( _particleIndex, particleIndex, config->getParticleCount() * 2 * sizeof( int ) );
	return 0;
}
unsigned int owOpenCLSolver::_runSortPostPass(owConfigProrerty * config)
{
	// Stage SortPostPass
	sortPostPass.setArg( 0, particleIndex );
	sortPostPass.setArg( 1, particleIndexBack );
	sortPostPass.setArg( 2, position );
	sortPostPass.setArg( 3, velocity );
	sortPostPass.setArg( 4, sortedPosition );
	sortPostPass.setArg( 5, sortedVelocity );
	sortPostPass.setArg( 6, config->getParticleCount()  );
	int err = queue.enqueueNDRangeKernel(
		sortPostPass, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_runIndexx(owConfigProrerty * config)
{
	// Stage Indexx
	indexx.setArg( 0, particleIndex );
	indexx.setArg( 1, config->gridCellCount );
	indexx.setArg( 2, gridCellIndex );
	indexx.setArg( 3, config->getParticleCount() );
	int gridCellCountRoundedUp = ((( config->gridCellCount - 1 ) / local_NDRange_size ) + 1 ) * local_NDRange_size;
	int err = queue.enqueueNDRangeKernel(
		indexx, cl::NullRange, cl::NDRange( (int) ( /**/gridCellCountRoundedUp/**/ ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_runIndexPostPass(owConfigProrerty * config)
{
	// Stage IndexPostPass
	//28aug_Palyanov_start_block
	copy_buffer_from_device( gridNextNonEmptyCellBuffer, gridCellIndex,(config->gridCellCount+1) * sizeof( unsigned int ) * 1 );
	int recentNonEmptyCell = config->gridCellCount;
	for(int i=config->gridCellCount;i>=0;i--)
	{
		if(gridNextNonEmptyCellBuffer[i] == NO_CELL_ID)
			gridNextNonEmptyCellBuffer[i] = recentNonEmptyCell; 
		else recentNonEmptyCell = gridNextNonEmptyCellBuffer[i];
	}
	int err = copy_buffer_to_device( gridNextNonEmptyCellBuffer,gridCellIndexFixedUp,(config->gridCellCount+1) * sizeof( unsigned int ) * 1 );
	return err;
}
unsigned int owOpenCLSolver::_runFindNeighbors(owConfigProrerty * config)
{
	// Stage FindNeighbors
	findNeighbors.setArg( 0, gridCellIndexFixedUp );
	findNeighbors.setArg( 1, sortedPosition );
	findNeighbors.setArg( 2, config->gridCellCount );
	findNeighbors.setArg( 3, config->gridCellsX );
	findNeighbors.setArg( 4, config->gridCellsY );
	findNeighbors.setArg( 5, config->gridCellsZ );
	findNeighbors.setArg( 6, h );
	findNeighbors.setArg( 7, hashGridCellSize );
	findNeighbors.setArg( 8, hashGridCellSizeInv );
	findNeighbors.setArg( 9, simulationScale );
	findNeighbors.setArg( 10, config->xmin );
	findNeighbors.setArg( 11, config->ymin );
	findNeighbors.setArg( 12, config->zmin );
	findNeighbors.setArg( 13, neighborMap );
	findNeighbors.setArg( 14, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(
		findNeighbors, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );/* 
		local_work_size can also be a NULL
		value in which case the OpenCL implementation will
		determine how to be break the global work-items 
		into appropriate work-group instances.
		http://www.khronos.org/registry/cl/specs/opencl-1.0.43.pdf, page 109 
		*/
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_computeDensity(owConfigProrerty * config)
{
	// Stage ComputeDensityPressure
	pcisph_computeDensity.setArg( 0, neighborMap );
	pcisph_computeDensity.setArg( 1, Wpoly6Coefficient );
	//pcisph_computeDensity.setArg( 2, gradWspikyCoefficient );
	pcisph_computeDensity.setArg( 2, h );
	pcisph_computeDensity.setArg( 3, mass );
	pcisph_computeDensity.setArg( 4, rho0 );
	pcisph_computeDensity.setArg( 5, simulationScale );
	pcisph_computeDensity.setArg( 6, stiffness );
	pcisph_computeDensity.setArg( 7, sortedPosition );
	pcisph_computeDensity.setArg( 8, pressure );
	pcisph_computeDensity.setArg( 9, rho );
	pcisph_computeDensity.setArg(10, particleIndexBack );
	pcisph_computeDensity.setArg(11, delta );
	pcisph_computeDensity.setArg(12, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(
		pcisph_computeDensity, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_computeForcesAndInitPressure(owConfigProrerty * config)
{
	pcisph_computeForcesAndInitPressure.setArg( 0, neighborMap );
	pcisph_computeForcesAndInitPressure.setArg( 1, rho );
	pcisph_computeForcesAndInitPressure.setArg( 2, pressure );
	pcisph_computeForcesAndInitPressure.setArg( 3, sortedPosition );
	pcisph_computeForcesAndInitPressure.setArg( 4, sortedVelocity );
	pcisph_computeForcesAndInitPressure.setArg( 5, acceleration );
	pcisph_computeForcesAndInitPressure.setArg( 6, particleIndexBack );
	pcisph_computeForcesAndInitPressure.setArg( 7, surfTensCoeff );
//	pcisph_computeForcesAndInitPressure.setArg( 7, Wpoly6Coefficient );
	pcisph_computeForcesAndInitPressure.setArg( 8, del2WviscosityCoefficient );
	pcisph_computeForcesAndInitPressure.setArg( 9, h );
	pcisph_computeForcesAndInitPressure.setArg(10, mass );
	pcisph_computeForcesAndInitPressure.setArg(11, viscosity );
	pcisph_computeForcesAndInitPressure.setArg(12, simulationScale );
	pcisph_computeForcesAndInitPressure.setArg(13, gravity_x );
	pcisph_computeForcesAndInitPressure.setArg(14, gravity_y );
	pcisph_computeForcesAndInitPressure.setArg(15, gravity_z );
	pcisph_computeForcesAndInitPressure.setArg(16, position );
	pcisph_computeForcesAndInitPressure.setArg(17, particleIndex );
	pcisph_computeForcesAndInitPressure.setArg(18, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(
		pcisph_computeForcesAndInitPressure, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}

unsigned int owOpenCLSolver::_run_pcisph_computeElasticForces(owConfigProrerty * config)
{
	if(numOfElasticP == 0 )
		return 0;
	pcisph_computeElasticForces.setArg( 0, neighborMap );
	pcisph_computeElasticForces.setArg( 1, sortedPosition );
	pcisph_computeElasticForces.setArg( 2, sortedVelocity );
	pcisph_computeElasticForces.setArg( 3, acceleration );
	pcisph_computeElasticForces.setArg( 4, particleIndexBack );
	pcisph_computeElasticForces.setArg( 5, velocity );
	pcisph_computeElasticForces.setArg( 6, h );
	pcisph_computeElasticForces.setArg( 7, mass );
	pcisph_computeElasticForces.setArg( 8, simulationScale );
	pcisph_computeElasticForces.setArg( 9, numOfElasticP );
	pcisph_computeElasticForces.setArg( 10, elasticConnectionsData );
	pcisph_computeElasticForces.setArg( 11, config->getParticleCount() );
	pcisph_computeElasticForces.setArg( 12, MUSCLE_COUNT );
	pcisph_computeElasticForces.setArg( 13, muscle_activation_signal);
	pcisph_computeElasticForces.setArg( 14, position);
	pcisph_computeElasticForces.setArg( 15, elasticityCoefficient);
	int numOfElasticPCountRoundedUp = ((( numOfElasticP - 1 ) / local_NDRange_size ) + 1 ) * local_NDRange_size;
	int err = queue.enqueueNDRangeKernel(
		pcisph_computeElasticForces, cl::NullRange, cl::NDRange( numOfElasticPCountRoundedUp ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}

unsigned int owOpenCLSolver::_run_pcisph_predictPositions(owConfigProrerty * config)
{
	pcisph_predictPositions.setArg( 0, acceleration );
	pcisph_predictPositions.setArg( 1, sortedPosition );
	pcisph_predictPositions.setArg( 2, sortedVelocity );
	pcisph_predictPositions.setArg( 3, particleIndex );
	pcisph_predictPositions.setArg( 4, particleIndexBack );
	pcisph_predictPositions.setArg( 5, gravity_x );
	pcisph_predictPositions.setArg( 6, gravity_y );
	pcisph_predictPositions.setArg( 7, gravity_z );
	pcisph_predictPositions.setArg( 8, simulationScaleInv );
	pcisph_predictPositions.setArg( 9, timeStep );
	pcisph_predictPositions.setArg( 10, config->xmin );
	pcisph_predictPositions.setArg( 11, config->xmax );
	pcisph_predictPositions.setArg( 12, config->ymin );
	pcisph_predictPositions.setArg( 13, config->ymax );
	pcisph_predictPositions.setArg( 14, config->zmin );
	pcisph_predictPositions.setArg( 15, config->zmax );
	pcisph_predictPositions.setArg( 16, damping );
	pcisph_predictPositions.setArg( 17, position );
	pcisph_predictPositions.setArg( 18, velocity );
	pcisph_predictPositions.setArg( 19, r0 );
	pcisph_predictPositions.setArg( 20, neighborMap );
	pcisph_predictPositions.setArg( 21, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(
		pcisph_predictPositions, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_predictDensity(owConfigProrerty * config)
{
	// Stage ComputeDensityPressure
	pcisph_predictDensity.setArg( 0, neighborMap );
	pcisph_predictDensity.setArg( 1, particleIndexBack );
	pcisph_predictDensity.setArg( 2, Wpoly6Coefficient );
	//pcisph_predictDensity.setArg( 3, gradWspikyCoefficient );
	pcisph_predictDensity.setArg( 3, h );
	pcisph_predictDensity.setArg( 4, mass );
	pcisph_predictDensity.setArg( 5, rho0 );
	pcisph_predictDensity.setArg( 6, simulationScale );
	pcisph_predictDensity.setArg( 7, stiffness );
	pcisph_predictDensity.setArg( 8, sortedPosition );
	pcisph_predictDensity.setArg( 9, pressure );
	pcisph_predictDensity.setArg(10, rho );
	pcisph_predictDensity.setArg(11, delta );
	pcisph_predictDensity.setArg(12, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(
		pcisph_predictDensity, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_correctPressure(owConfigProrerty * config)
{
	// Stage ComputeDensityPressure
	pcisph_correctPressure.setArg( 0, neighborMap );
	pcisph_correctPressure.setArg( 1, particleIndexBack );
	//pcisph_correctPressure.setArg( 2, gradWspikyCoefficient );
	pcisph_correctPressure.setArg( 2, h );
	pcisph_correctPressure.setArg( 3, mass );
	pcisph_correctPressure.setArg( 4, rho0 );
	pcisph_correctPressure.setArg( 5, simulationScale );
	pcisph_correctPressure.setArg( 6, stiffness );
	pcisph_correctPressure.setArg( 7, sortedPosition );
	pcisph_correctPressure.setArg( 8, pressure );
	pcisph_correctPressure.setArg( 9, rho );
	pcisph_correctPressure.setArg(10, delta );
	pcisph_correctPressure.setArg(11, position );
	pcisph_correctPressure.setArg(12, particleIndex );
	pcisph_correctPressure.setArg(13, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(
		pcisph_correctPressure, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_computePressureForceAcceleration(owConfigProrerty * config)
{
	// Stage ComputeAcceleration
	pcisph_computePressureForceAcceleration.setArg( 0, neighborMap );
	pcisph_computePressureForceAcceleration.setArg( 1, pressure );
	pcisph_computePressureForceAcceleration.setArg( 2, rho );
	pcisph_computePressureForceAcceleration.setArg( 3, sortedPosition );
	pcisph_computePressureForceAcceleration.setArg( 4, sortedVelocity );
	pcisph_computePressureForceAcceleration.setArg( 5, particleIndexBack );
	pcisph_computePressureForceAcceleration.setArg( 6, delta );
	//pcisph_computePressureForceAcceleration.setArg( 7, del2WviscosityCoefficient );
	pcisph_computePressureForceAcceleration.setArg( 7, gradWspikyCoefficient );
	pcisph_computePressureForceAcceleration.setArg( 8, h );
	pcisph_computePressureForceAcceleration.setArg(  9, mass );
	pcisph_computePressureForceAcceleration.setArg( 10, viscosity );
	pcisph_computePressureForceAcceleration.setArg( 11, simulationScale );
	pcisph_computePressureForceAcceleration.setArg( 12, acceleration );
	pcisph_computePressureForceAcceleration.setArg( 13, rho0 );
	pcisph_computePressureForceAcceleration.setArg( 14, position );
	pcisph_computePressureForceAcceleration.setArg( 15, particleIndex );
	pcisph_computePressureForceAcceleration.setArg( 16, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(
		pcisph_computePressureForceAcceleration, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}

unsigned int owOpenCLSolver::_run_clearMembraneBuffers(owConfigProrerty * config)
{
	clearMembraneBuffers.setArg( 0, position );
	clearMembraneBuffers.setArg( 1, velocity );
	clearMembraneBuffers.setArg( 2, sortedPosition );
	clearMembraneBuffers.setArg( 3, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(
		clearMembraneBuffers, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}

unsigned int owOpenCLSolver::_run_computeInteractionWithMembranes(owConfigProrerty * config)
{
	computeInteractionWithMembranes.setArg( 0, position );
	computeInteractionWithMembranes.setArg( 1, velocity );
	computeInteractionWithMembranes.setArg( 2, sortedPosition );
	computeInteractionWithMembranes.setArg( 3, particleIndex );
	computeInteractionWithMembranes.setArg( 4, particleIndexBack );
	computeInteractionWithMembranes.setArg( 5, neighborMap );
	computeInteractionWithMembranes.setArg( 6, particleMembranesList );
	computeInteractionWithMembranes.setArg( 7, membraneData );
	computeInteractionWithMembranes.setArg( 8, config->getParticleCount() );
	computeInteractionWithMembranes.setArg( 9, numOfElasticP );
	computeInteractionWithMembranes.setArg(10, r0 );
	int err = queue.enqueueNDRangeKernel(
		computeInteractionWithMembranes, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}

unsigned int owOpenCLSolver::_run_computeInteractionWithMembranes_finalize(owConfigProrerty * config)
{
	computeInteractionWithMembranes_finalize.setArg( 0, position );
	computeInteractionWithMembranes_finalize.setArg( 1, velocity );
	computeInteractionWithMembranes_finalize.setArg( 2, particleIndex );
	computeInteractionWithMembranes_finalize.setArg( 3, particleIndexBack );
	computeInteractionWithMembranes_finalize.setArg( 4, config->getParticleCount() );
	int err = queue.enqueueNDRangeKernel(
		computeInteractionWithMembranes_finalize, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}


unsigned int owOpenCLSolver::_run_pcisph_integrate(int iterationCount, owConfigProrerty * config)
{
	// Stage Integrate
	pcisph_integrate.setArg( 0, acceleration );
	pcisph_integrate.setArg( 1, sortedPosition );
	pcisph_integrate.setArg( 2, sortedVelocity );
	pcisph_integrate.setArg( 3, particleIndex );
	pcisph_integrate.setArg( 4, particleIndexBack );
	pcisph_integrate.setArg( 5, gravity_x );
	pcisph_integrate.setArg( 6, gravity_y );
	pcisph_integrate.setArg( 7, gravity_z );
	pcisph_integrate.setArg( 8, simulationScaleInv );
	pcisph_integrate.setArg( 9, timeStep );
	pcisph_integrate.setArg( 10, config->xmin );
	pcisph_integrate.setArg( 11, config->xmax );
	pcisph_integrate.setArg( 12, config->ymin );
	pcisph_integrate.setArg( 13, config->ymax );
	pcisph_integrate.setArg( 14, config->zmin );
	pcisph_integrate.setArg( 15, config->zmax );
	pcisph_integrate.setArg( 16, damping );
	pcisph_integrate.setArg( 17, position );
	pcisph_integrate.setArg( 18, velocity );
	pcisph_integrate.setArg( 19, rho );
	pcisph_integrate.setArg( 20, r0 );
	pcisph_integrate.setArg( 21, neighborMap );
	pcisph_integrate.setArg( 22, config->getParticleCount() );
	pcisph_integrate.setArg( 23, iterationCount );
	int err = queue.enqueueNDRangeKernel(
		pcisph_integrate, cl::NullRange, cl::NDRange( (int) (  config->getParticleCount_RoundUp() ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( local_NDRange_size ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
//end Kernels definition
//Auxiliary methods
int myCompare( const void * v1, const void * v2 ){
	int * f1 = (int *)v1;
	int * f2 = (int *)v2;
	if( f1[ 0 ] < f2[ 0 ] ) return -1;
	if( f1[ 0 ] > f2[ 0 ] ) return +1;
	return 0;
}
void owOpenCLSolver::create_ocl_kernel(const char *name, cl::Kernel &k )
{
	int err;
	k = cl::Kernel(program, name, &err);
	if( err != CL_SUCCESS ){
		std::string error_m = "Kernel creation failed: ";
		error_m.append(name);
		throw std::runtime_error( error_m );
	}
}
void owOpenCLSolver::create_ocl_buffer(const char *name, cl::Buffer &b, const cl_mem_flags flags,const int size)
{
	int err;
	b = cl::Buffer(context, flags, size, NULL, &err);
	if( err != CL_SUCCESS ){
		std::string error_m = "Buffer creation failed: ";
		error_m.append(name);
		throw std::runtime_error( error_m );
	}
}
int owOpenCLSolver::copy_buffer_to_device(const void *host_b, cl::Buffer &ocl_b, const int size )
{
	//Actualy we should check  size and type 
	int err = queue.enqueueWriteBuffer( ocl_b, CL_TRUE, 0, size, host_b );
	if( err != CL_SUCCESS ){
		throw std::runtime_error( "Could not enqueue write" );
	}
	queue.finish();
	return err;
}
int owOpenCLSolver::copy_buffer_from_device(void *host_b, const cl::Buffer &ocl_b, const int size )
{
	//Actualy we should check  size and type 
	int err = queue.enqueueReadBuffer( ocl_b, CL_TRUE, 0, size, host_b );
	if( err != CL_SUCCESS ){
		throw std::runtime_error( "Could not enqueue read" );
	}
	queue.finish();
	return err;
}

unsigned int owOpenCLSolver::updateMuscleActivityData(float *_muscle_activation_signal_cpp)
{
	copy_buffer_to_device( _muscle_activation_signal_cpp, muscle_activation_signal, MUSCLE_COUNT * sizeof( float ) );
	return 0;
}

owOpenCLSolver::~owOpenCLSolver(void)
{
	delete [] gridNextNonEmptyCellBuffer;
	delete [] _particleIndex;
}
