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
#include <sstream>

#include "owOpenCLSolver.h"

int myCompare( const void * v1, const void * v2 );

/** Constructor of class owOpenCLSolver
 *
 *  @param position_cpp
 *  initial position buffer
 *  @param velocity_cpp
 *  initial velocity buffer
 *  @param config
 *  Contain information about simulating configuration
 *  @param elasticConnectionData_cpp
 *  buffer with info about elastic connections
 *  @param membraneData_cpp
 *  buffer with info about membranes
 *  @param particleMembranesList_cpp
 *  buffer with info about sets of membranes in which particular particle is including
 */
owOpenCLSolver::owOpenCLSolver(const float * position_cpp, const float * velocity_cpp, owConfigProperty * config, const float * elasticConnectionsData_cpp, const int * membraneData_cpp, const int * particleMembranesList_cpp)
{
	try{
		initializeOpenCL(config);
		// Create OpenCL buffers
		initializeBuffers(position_cpp, velocity_cpp, config, elasticConnectionsData_cpp, membraneData_cpp, particleMembranesList_cpp);
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
	}catch(std::runtime_error & ex){
		destroy();
		throw ex;
	}
}

/** Reset simulation method
 *
 *  This Method reset all simulation. It's reiniting all buffers with
 *  new start configuration and restart simulation.
 *  NOTE: this method reiniting only buffers not a kernels an so one.
 *
 *  @param position_cpp
 *  initial position buffer
 *  @param velocity_cpp
 *  initial velocity buffer
 *  @param config
 *  Contain information about simulating configuration
 *  @param elasticConnectionData_cpp
 *  buffer with info about elastic connections
 *  @param membraneData_cpp
 *  buffer with info about membranes
 *  @param particleMembranesList_cpp
 *  buffer with info about sets of membranes in which particular particle is including
 */
void owOpenCLSolver::reset(const float * position_cpp, const float * velocity_cpp, owConfigProperty * config, const float * elasticConnectionsData_cpp, const int * membraneData_cpp, const int * particleMembranesList_cpp){
	// Reinitializing all data buffer
	// First freed data for buffer _particleIndex and gridNextNonEmptyCellBuffer
	destroy(); // Clear buffers before relocate it again this is needed because simulation can run different configuration
	initializeBuffers(position_cpp,velocity_cpp, config, elasticConnectionsData_cpp, membraneData_cpp, particleMembranesList_cpp);

}
/** Initialization of data buffers on OpenCL device
*  @param position_cpp
*  initial position buffer
*  @param velocity_cpp
*  initial velocity buffer
*  @param config
*  Contain information about simulating configuration
*  @param elasticConnectionData_cpp
*  buffer with info about elastic connections
*  @param membraneData_cpp
*  buffer with info about membranes
*  @param particleMembranesList_cpp
*  buffer with info about sets of membranes in which particular particle is including
*/
void owOpenCLSolver::initializeBuffers(const float * position_cpp, const float * velocity_cpp, owConfigProperty * config, const float * elasticConnectionsData_cpp, const int * membraneData_cpp, const int * particleMembranesList_cpp){
	// Needed for sortin stuff
	_particleIndex = new   int[ 2 * config->getParticleCount() ];
	gridNextNonEmptyCellBuffer = new int[config->gridCellCount+1];

	create_ocl_buffer( "acceleration", acceleration, CL_MEM_READ_WRITE, ( config->getParticleCount() * sizeof( float ) * 4 * 3 ) );// 4*2-->4*3; third part is to store acceleration[t], while first to are for acceleration[t+delta_t]
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
	create_ocl_buffer( "muscle_activation_signal", muscle_activation_signal, CL_MEM_READ_WRITE, ( config->MUSCLE_COUNT * sizeof( float ) ) );

	if(membraneData_cpp != NULL && particleMembranesList_cpp != NULL)
	{
		create_ocl_buffer( "membraneData", membraneData, CL_MEM_READ_WRITE, ( config->numOfMembranes * sizeof( int ) * 3 ) );
		create_ocl_buffer("particleMembranesList", particleMembranesList,CL_MEM_READ_WRITE, config->numOfElasticP * MAX_MEMBRANES_INCLUDING_SAME_PARTICLE * sizeof(int) );
	}
	if(elasticConnectionsData_cpp != NULL){
		create_ocl_buffer("elasticConnectionsData", elasticConnectionsData,CL_MEM_READ_WRITE, config->numOfElasticP * MAX_NEIGHBOR_COUNT * sizeof(float) * 4);
	}

	// Copy initial position_cpp and velocity_cpp to the OpenCL Device
	copy_buffer_to_device( position_cpp, position, config->getParticleCount() * sizeof( float ) * 4 );
	copy_buffer_to_device( velocity_cpp, velocity, config->getParticleCount() * sizeof( float ) * 4 );
	// Copy membrane data to device memory elastic connections
	if(membraneData_cpp != NULL && particleMembranesList_cpp != NULL)
	{
		copy_buffer_to_device( membraneData_cpp, membraneData, config->numOfMembranes * sizeof( int ) * 3 );
		copy_buffer_to_device( particleMembranesList_cpp, particleMembranesList, config->numOfElasticP * MAX_MEMBRANES_INCLUDING_SAME_PARTICLE * sizeof( int ) );
	}
	// Copy elastic connectiondate to device memory elastic connections
	if(elasticConnectionsData_cpp != NULL){
		copy_buffer_to_device(elasticConnectionsData_cpp, elasticConnectionsData, config->numOfElasticP * MAX_NEIGHBOR_COUNT * sizeof(float) * 4);
	}
}

/** Initialization OpenCL entities
 *
 * Method inits all required entities for work with OpenCL code.
 *
 *  @param config
 *  Contain information about simulating configuration
 */
void owOpenCLSolver::initializeOpenCL(owConfigProperty * config)
{
	cl_int err;
	std::vector< cl::Platform > platformList;
	err = cl::Platform::get( &platformList ); //TODO make check that returned value isn't error
	if( platformList.size() < 1 || err != CL_SUCCESS ){
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
	unsigned int device_type [] = {CL_DEVICE_TYPE_CPU,CL_DEVICE_TYPE_GPU, CL_DEVICE_TYPE_ALL};

	int plList = -1;//selected platform index in platformList array [choose CPU by default]
							//added autodetection of device number corresonding to preferrable device type (CPU|GPU) | otherwise the choice will be made from list of existing devices
	cl_uint ciDeviceCount = 0;
	cl_device_id * devices_t;
	bool bPassed = true, findDevice = false;
	cl_int result;
	cl_uint device_coumpute_unit_num;
	cl_uint device_coumpute_unit_num_current = 0;
	unsigned int deviceNum = 0;
	//Selection of more appropriate device
	while(!findDevice){
		for(int clSelectedPlatformID = 0;clSelectedPlatformID < (int)n_pl;clSelectedPlatformID++){
			//if(findDevice)
			//	break;
			clGetDeviceIDs (cl_pl_id[clSelectedPlatformID], device_type[config->getDeviceType()], 0, NULL, &ciDeviceCount);
			if((devices_t = static_cast<cl_device_id *>(malloc(sizeof(cl_device_id) * ciDeviceCount))) == NULL)
				bPassed = false;
			if(bPassed){
				result= clGetDeviceIDs (cl_pl_id[clSelectedPlatformID], device_type[config->getDeviceType()], ciDeviceCount, devices_t, &ciDeviceCount);
				if( result == CL_SUCCESS){
					for( cl_uint i =0; i < ciDeviceCount; ++i ){
						clGetDeviceInfo(devices_t[i], CL_DEVICE_TYPE, sizeof(type), &type, NULL);
						if( type & device_type[config->getDeviceType()]){
							clGetDeviceInfo(devices_t[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(device_coumpute_unit_num), &device_coumpute_unit_num, NULL);
							if(device_coumpute_unit_num_current <= device_coumpute_unit_num){
								plList = clSelectedPlatformID;
								device_coumpute_unit_num_current = device_coumpute_unit_num;
								findDevice = true;
								deviceNum = i;
							}
							//break;
						}
					}
				}
				free(devices_t);
			}
		}
		if(!findDevice){
			//plList = 0;
			deviceNum = 0;
			std::string deviceTypeName = (config->getDeviceType() == ALL)? "ALL": (config->getDeviceType() == CPU)? "CPU":"GPU";
			std::cout << "Unfortunately OpenCL couldn't find device " << deviceTypeName << std::endl;
			std::cout << "OpenCL try to init existing device " << std::endl;
			if(config->getDeviceType() != ALL)
				config->setDeviceType(ALL);
			else
				throw std::runtime_error("Sibernetic can't find any OpenCL devices. Please check you're environment configuration.");
		}
	}
	cl_context_properties cprops[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties) (platformList[plList])(), 0 };
	context = cl::Context( device_type[config->getDeviceType()], cprops, NULL, NULL, &err );
	devices = context.getInfo< CL_CONTEXT_DEVICES >();
	if( devices.size() < 1 ){
		throw std::runtime_error( "No OpenCL devices found" );
	}
	//Print some information about chosen platform
	int value;
	unsigned long val2;
	size_t val3;
	result = devices[deviceNum].getInfo(CL_DEVICE_NAME,&cBuffer);// CL_INVALID_VALUE = -30;
	if(result == CL_SUCCESS) std::cout << "CL_CONTEXT_PLATFORM ["<< plList << "]: CL_DEVICE_NAME [" << deviceNum << "]:\t" << cBuffer << "\n" << std::endl;
	if(strlen(cBuffer)<1000) config->setDeviceName(cBuffer);
	result = devices[deviceNum].getInfo(CL_DEVICE_TYPE,&cBuffer);
	if(result == CL_SUCCESS) std::cout << "CL_CONTEXT_PLATFORM ["<< plList << "]: CL_DEVICE_TYPE [" << deviceNum << "]:\t" << (((int)cBuffer[0] == CL_DEVICE_TYPE_CPU)? "CPU" : "GPU") << std::endl;
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
		throw std::runtime_error( "Failed to create command queue" );
	}
	std::ifstream file( config->getSourceFileName().c_str() );
	if( !file.is_open() ){
		throw std::runtime_error( "Could not open file with OpenCL program check input arguments oclsourcepath: " + config->getSourceFileName() );
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
		err = program.build( devices, "");
	#endif
#endif
	if( err != CL_SUCCESS ){
		std::string compilationErrors;
		compilationErrors = program.getBuildInfo< CL_PROGRAM_BUILD_LOG >( devices[ 0 ] );
		std::cerr << "Compilation failed: " << std::endl << compilationErrors << std::endl;
		throw std::runtime_error( "failed to build program" );
	}
	std::cout<<"OPENCL program was successfully build. Program file oclsourcepath: " << config->getSourceFileName() << std::endl;
	return;
}
/** Create error message
 */
std::string errorMessage(const char * s, int & error){
	std::stringstream ss;
	ss << s;
	ss << " error code is ";
	ss << error;
	return ss.str();
}
//Kernels functions definition
/** Run clearing neighbor map
 *
 *  This function is depreciated and will be removed in next releases.
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_runClearBuffers(owConfigProperty * config)
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _runClearBuffers ", err ));
	}
	return err;
}
/** Run hash particle kernel
 *
 *  This kernel runs algorithm for calculating cell number
 *  for all particles
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_runHashParticles(owConfigProperty * config)
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _runHashParticles", err ));
	}
	return err;
}
/** Run indexx kernel
 *
 *  Fill up gridCellIndex array with start correct
 *  value in partcileIndex array.
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_runIndexx(owConfigProperty * config)
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _runIndexx", err ));
	}
	return err;
}
/** Run index post pass function
 *
 *  Fill up all empty cell in gridCellIndex. If value in particular cell == -1
 *  it fills with value from last non empty cell. It need for optimization
 *  of neighbor search.
 *  EXAMPLE: particleIndex after sorting [[1,1],[1,2],[2,3],[3,0],[3,7],[4,5],[6,8]...]
 *                                          ^           ^     ^           ^     ^
 *           gridCellIndex               [  0,          2,    3,          5,-1, 6,....]
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to write buffet to a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueWriteBuffer.html)
 */
void owOpenCLSolver::_runIndexPostPass(owConfigProperty * config)
{
	//Stage IndexPostPass
	//28aug_Palyanov_start_block
	copy_buffer_from_device( gridNextNonEmptyCellBuffer, gridCellIndex,(config->gridCellCount+1) * sizeof( int ) * 1 );
	int recentNonEmptyCell = config->gridCellCount;
	for(int i=config->gridCellCount;i>=0;i--)
	{
		if(gridNextNonEmptyCellBuffer[i] == NO_CELL_ID)
			gridNextNonEmptyCellBuffer[i] = recentNonEmptyCell;
		else recentNonEmptyCell = gridNextNonEmptyCellBuffer[i];
	}
	copy_buffer_to_device( gridNextNonEmptyCellBuffer,gridCellIndexFixedUp,(config->gridCellCount+1) * sizeof( int ) * 1 );
}
/** Sorting of particleIndex list
 *
 *  particleIndex list contains for every particle information about current
 *  cell in which it contains in format {cell_id, partcile_id}. Before sorting
 *  particleIndex list arranged by particle id. This method sort particleIndex
 *  in accordance with order of cell.
 *  EXAMPLE: before sorting [[3,0],[1,1],[1,2],[2,3],..]
 *           after sorting  [[1,1],[1,2],[2,3],[3,0],..]
 *  For sorting method's using standard library quick sort algorithm.
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
void owOpenCLSolver::_runSort(owConfigProperty * config)
{
	copy_buffer_from_device( _particleIndex, particleIndex, config->getParticleCount() * 2 * sizeof( int ) );
	qsort( _particleIndex, config->getParticleCount(), 2 * sizeof( int ), myCompare );
	copy_buffer_to_device( _particleIndex, particleIndex, config->getParticleCount() * 2 * sizeof( int ) );
}
/** Run sorting of position and velocity arrays
 *
 *  After sorting particleIndex list's order isn't correspond
 *  to order of particles in position and velocity arrays.
 *  It also need to be arranged in correct order.
 *  NOTE: this kernel doesn't change an order of position and velocity buffers,
 *  it's using special created buffers sortedPosition and sortedVelocity.
 *  NOTE: for tracking info about particular particle using particleIndexBack list
 *  which == particleIndex before sorting.
 *  EXAMPLE:  particleIndex list after sorting  [[1,1],[1,2],[2,3],[3,0],..]
 *            position list                     [ pos_0, pos_1, pos_2, pos_3, ...]
 *            sortingPosition list after sorting[ pos_1, pos_2, pos_3, pos_0, ...]
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_runSortPostPass(owConfigProperty * config)
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _runSortPostPass", err ));
	}
	return err;
}
/** Run search for neighbors kernel
 *
 *  After preparing all required data: particleIndex list,
 *  sortedPosition and sortedVelocity, neighbor search starting.
 *  Kernel looking for a nearest particles for all particles.
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_runFindNeighbors(owConfigProperty * config)
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _runFindNeighbors", err ));
	}
	return err;
}
/** Run pcisph_computeDensity kernel
 *
 *  The kernel's calculating density for every particle.
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_pcisph_computeDensity(owConfigProperty * config)
{
	// Stage ComputeDensityPressure
	pcisph_computeDensity.setArg( 0, neighborMap );
	pcisph_computeDensity.setArg( 1, mass_mult_Wpoly6Coefficient );
	pcisph_computeDensity.setArg( 2, _hScaled2 );
	pcisph_computeDensity.setArg( 3, rho );
	pcisph_computeDensity.setArg( 4, particleIndexBack );
	pcisph_computeDensity.setArg( 5, config->getParticleCount() );
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_pcisph_computeDensity", err ));
	}
	return err;
}
/** Run pcisph_computeForcesAndInitPressure kernel
 *
 *  The kernel initializes pressure by 0.
 *  Calculating viscosity and surface tension forces
 *  and acceleration of particle
 *  acceleration[id] = (ViscosityForces + SurfaceTensiion +GravityForces)/mass
 *  tempAcceleration[id] = acceleration[id + PARTICLE_COUNT]=0.
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_pcisph_computeForcesAndInitPressure(owConfigProperty * config)
{
	pcisph_computeForcesAndInitPressure.setArg( 0, neighborMap );
	pcisph_computeForcesAndInitPressure.setArg( 1, rho );
	pcisph_computeForcesAndInitPressure.setArg( 2, pressure );
	pcisph_computeForcesAndInitPressure.setArg( 3, sortedPosition );
	pcisph_computeForcesAndInitPressure.setArg( 4, sortedVelocity );
	pcisph_computeForcesAndInitPressure.setArg( 5, acceleration );
	pcisph_computeForcesAndInitPressure.setArg( 6, particleIndexBack );
	pcisph_computeForcesAndInitPressure.setArg( 7, surfTensCoeff );
	pcisph_computeForcesAndInitPressure.setArg( 8, mass_mult_divgradWviscosityCoefficient );
	pcisph_computeForcesAndInitPressure.setArg( 9, _hScaled );
	pcisph_computeForcesAndInitPressure.setArg(10, viscosity );
	pcisph_computeForcesAndInitPressure.setArg(11, gravity_x );
	pcisph_computeForcesAndInitPressure.setArg(12, gravity_y );
	pcisph_computeForcesAndInitPressure.setArg(13, gravity_z );
	pcisph_computeForcesAndInitPressure.setArg(14, position );
	pcisph_computeForcesAndInitPressure.setArg(15, particleIndex );
	pcisph_computeForcesAndInitPressure.setArg(16, config->getParticleCount() );
	pcisph_computeForcesAndInitPressure.setArg(17, mass );
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_pcisph_computeForcesAndInitPressure", err ));
	}
	return err;
}
/** Run pcisph_computeElasticForces kernel
 *
 *  The kernel calculates elastic forces and muscle
 *  contraction forces if particle has muscle connection
 *  acceleration[id] += (ElasticForces + MuscleForce)/mass.
 *  NOTE: this kernel works only with elastic particles
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_pcisph_computeElasticForces(owConfigProperty * config)
{
	if(config->numOfElasticP == 0 )
		return 0;
	pcisph_computeElasticForces.setArg( 0, neighborMap );
	pcisph_computeElasticForces.setArg( 1, sortedPosition );
	pcisph_computeElasticForces.setArg( 2, sortedVelocity );
	pcisph_computeElasticForces.setArg( 3, acceleration );
	pcisph_computeElasticForces.setArg( 4, particleIndexBack );
	pcisph_computeElasticForces.setArg( 5, particleIndex );
	pcisph_computeElasticForces.setArg( 6, h );
	pcisph_computeElasticForces.setArg( 7, mass );
	pcisph_computeElasticForces.setArg( 8, simulationScale );
	pcisph_computeElasticForces.setArg( 9, config->numOfElasticP );
	pcisph_computeElasticForces.setArg( 10, elasticConnectionsData );
	pcisph_computeElasticForces.setArg( 11, config->getParticleCount() );
	pcisph_computeElasticForces.setArg( 12, config->MUSCLE_COUNT );
	pcisph_computeElasticForces.setArg( 13, muscle_activation_signal);
	pcisph_computeElasticForces.setArg( 14, position);
	pcisph_computeElasticForces.setArg( 15, elasticityCoefficient);
	int numOfElasticPCountRoundedUp = ((( config->numOfElasticP - 1 ) / local_NDRange_size ) + 1 ) * local_NDRange_size;
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_pcisph_computeElasticForces", err ));
	}
	return err;
}
/** Run pcisph_predictPositions kernel
 *
 *  The kernel predicts possible position value of particles
 *  what leads to incompressibility. Temp value of position
 *  is calculating from temp value of velocity which's tacking from predicted value of
 *  tempacceleration[id].
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_pcisph_predictPositions(owConfigProperty * config)
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
	pcisph_predictPositions.setArg( 9, config->getTimeStep() );
	pcisph_predictPositions.setArg(10, position );
	pcisph_predictPositions.setArg(11, velocity );
	pcisph_predictPositions.setArg(12, r0 );
	pcisph_predictPositions.setArg(13, neighborMap );
	pcisph_predictPositions.setArg(14, config->getParticleCount() );
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_pcisph_predictPositions", err ));
	}
	return err;
}
/** Run pcisph_predictDensity kernel
 *
 *  The kernel predicts possible value of density
 *  taking into account predicted value of particle's position
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_pcisph_predictDensity(owConfigProperty * config)
{
	// Stage ComputeDensityPressure
	pcisph_predictDensity.setArg( 0, neighborMap );
	pcisph_predictDensity.setArg( 1, particleIndexBack );
	pcisph_predictDensity.setArg( 2, mass_mult_Wpoly6Coefficient );
	pcisph_predictDensity.setArg( 3, h );
	pcisph_predictDensity.setArg( 4, rho0 );
	pcisph_predictDensity.setArg( 5, simulationScale );
	pcisph_predictDensity.setArg( 6, sortedPosition );
	pcisph_predictDensity.setArg( 7, pressure );
	pcisph_predictDensity.setArg( 8, rho );
	pcisph_predictDensity.setArg( 9, config->getParticleCount() );
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_pcisph_predictDensity", err ));
	}
	return err;
}
/** Run pcisph_correctPressure kernel
 *
 *  The kernel corrects the pressure
 *  taking into account predicted values of density.
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_pcisph_correctPressure(owConfigProperty * config)
{
	// Stage ComputeDensityPressure
	pcisph_correctPressure.setArg( 0, particleIndexBack );
	pcisph_correctPressure.setArg( 1, rho0 );
	pcisph_correctPressure.setArg( 2, pressure );
	pcisph_correctPressure.setArg( 3, rho );
	pcisph_correctPressure.setArg( 4, config->getDelta() );
	pcisph_correctPressure.setArg( 5, config->getParticleCount() );
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_pcisph_correctPressure", err ));
	}
	return err;
}
/** Run pcisph_computePressureForceAcceleration kernel
 *
 *  The kernel calculating pressure forces
 *  and calculating new value of tempacceleration[id].
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_pcisph_computePressureForceAcceleration(owConfigProperty * config)
{
	// Stage ComputeAcceleration
	pcisph_computePressureForceAcceleration.setArg( 0, neighborMap );
	pcisph_computePressureForceAcceleration.setArg( 1, pressure );
	pcisph_computePressureForceAcceleration.setArg( 2, rho );
	pcisph_computePressureForceAcceleration.setArg( 3, sortedPosition );
	pcisph_computePressureForceAcceleration.setArg( 4, sortedVelocity );
	pcisph_computePressureForceAcceleration.setArg( 5, particleIndexBack );
	pcisph_computePressureForceAcceleration.setArg( 6, config->getDelta() );
	pcisph_computePressureForceAcceleration.setArg( 7, mass_mult_gradWspikyCoefficient );
	pcisph_computePressureForceAcceleration.setArg( 8, h );
	pcisph_computePressureForceAcceleration.setArg( 9, simulationScale );
	pcisph_computePressureForceAcceleration.setArg(10, viscosity );
	pcisph_computePressureForceAcceleration.setArg(11, acceleration );
	pcisph_computePressureForceAcceleration.setArg(12, rho0 );
	pcisph_computePressureForceAcceleration.setArg(13, position );
	pcisph_computePressureForceAcceleration.setArg(14, particleIndex );
	pcisph_computePressureForceAcceleration.setArg(15, config->getParticleCount());
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_pcisph_computePressureForceAcceleration", err ));
	}
	return err;
}
/** Run clearMembraneBuffers kernel
 *
 *  The kernel clears membranes data buffers.
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_clearMembraneBuffers(owConfigProperty * config)
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_clearMembraneBuffers", err ));
	}
	return err;
}
/** Run computeInteractionWithMembranes kernel
 *
 *  The kernel handles interaction particles and membranes.
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_computeInteractionWithMembranes(owConfigProperty * config)
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
	computeInteractionWithMembranes.setArg( 9, config->numOfElasticP );
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_computeInteractionWithMembranes", err ));
	}
	return err;
}
/** Run computeInteractionWithMembranes_finalize kernel
 *
 *  The kernel corrects position and velocity
 *  of particles after interaction with membranes.
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_computeInteractionWithMembranes_finalize(owConfigProperty * config)
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_computeInteractionWithMembranes_finalize", err ));
	}
	return err;
}
/** Run pcisph_integrate kernel
 *
 *  The kernel run numerical integration method.
 *  Calculating value of position and velocity on step (t+1)
 *  NOTE: for now simulation using Semi-implicit Euler method
 *        for integration 1th order
 *  NOTE: soon we plan to add Leap-frog 2th order
 *
 *  @param config
 *  Contain information about simulating configuration
 *  @return value taking after enqueue a command to execute a kernel on a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html)
 */
unsigned int owOpenCLSolver::_run_pcisph_integrate(int iterationCount, int pcisph_integrate_mode, owConfigProperty * config)
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
	pcisph_integrate.setArg( 9, config->getTimeStep() );
	pcisph_integrate.setArg( 10, config->xmin );
	pcisph_integrate.setArg( 11, config->xmax );
	pcisph_integrate.setArg( 12, config->ymin );
	pcisph_integrate.setArg( 13, config->ymax );
	pcisph_integrate.setArg( 14, config->zmin );
	pcisph_integrate.setArg( 15, config->zmax );
	pcisph_integrate.setArg( 16, position );
	pcisph_integrate.setArg( 17, velocity );
	pcisph_integrate.setArg( 18, rho );
	pcisph_integrate.setArg( 19, r0 );
	pcisph_integrate.setArg( 20, neighborMap );
	pcisph_integrate.setArg( 21, config->getParticleCount() );
	pcisph_integrate.setArg( 22, iterationCount );
	pcisph_integrate.setArg( 23, pcisph_integrate_mode );
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
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("An ERROR is appearing during work of kernel _run_pcisph_integrate", err ));
	}
	return err;
}
//end Kernels definition
//Auxiliary methods
/** Comparator method
 *
 *  This needs for quick sort standard method.
 *  More info (http://www.cplusplus.com/reference/cstdlib/qsort/).
 *
 *  @param v1
 *  @param v2
 *  @return -1 if value v1[0] > v2[0], +1 v1[0] < v2[0] else 0.
 */
int myCompare( const void * v1, const void * v2 ){
	const int * f1 = static_cast<const int *>(v1);
	const int * f2 = static_cast<const int *>(v2);
	if( f1[ 0 ] < f2[ 0 ] ) return -1;
	if( f1[ 0 ] > f2[ 0 ] ) return +1;
	return 0;
}
/** OpenCL kernel creator method
 *
 *  @param name
 *  Name of kernel should be the same that in OpenCL program file
 *  @param k
 *  reference to cl::Kernel object
 */
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
/** OpenCL buffer creator method
 *
 *  @param name
 *  name of buffer
 *  @param b
 *  reference to cl::Buffer object
 *  @param flags
 *  OpenCL memory flags more info (https://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/enums.html)
 */
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
/** Copy openCL buffer from host program into a memory of device
 *
 *  @param host_b
 *  host buffer
 *  @param ocl_b
 *  reference to cl::Buffer object in which
 *  @param flags
 *  OpenCL memory flags more info (https://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/enums.html)
 *  @param size
 *  buffer's size
 *  @return value taking after enqueue a command to write buffet to a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueWriteBuffer.html)
 */
void owOpenCLSolver::copy_buffer_to_device(const void *host_b, cl::Buffer &ocl_b, const int size )
{
	//Actualy we should check  size and type
	int err = queue.enqueueWriteBuffer( ocl_b, CL_TRUE, 0, size, host_b );
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("Could not enqueue read data from buffer  error code is", err));
	}
	queue.finish();
}
/** Copy openCL buffer from device to host program memory
 *
 *  @param host_b
 *  host buffer
 *  @param ocl_b
 *  reference to cl::Buffer object in which
 *  OpenCL memory flags more info (https://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/enums.html)
 *  @param size
 *  buffer's size
 *  @return value taking after enqueue a command to read buffet from a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueReadBuffer.html)
 */
void owOpenCLSolver::copy_buffer_from_device(void *host_b, const cl::Buffer &ocl_b, const int size )
{
	//Actualy we should check  size and type
	int err = queue.enqueueReadBuffer( ocl_b, CL_TRUE, 0, size, host_b );
	if( err != CL_SUCCESS ){
		throw std::runtime_error( errorMessage("Could not enqueue read data from buffer  error code is", err) );
	}
	queue.finish();
}
/** Copy openCL buffer from device to host program memory
 *
 *  @param host_b
 *  host buffer
 *  @param ocl_b
 *  reference to cl::Buffer object in which
 *  OpenCL memory flags more info (https://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/enums.html)
 *  @param size
 *  buffer's size
 *  @return value taking after enqueue a command to read buffet from a device.
 *  More info here (http://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueReadBuffer.html)
 */
void owOpenCLSolver::updateMuscleActivityData(float *_muscle_activation_signal_cpp, owConfigProperty * config)
{
	copy_buffer_to_device( _muscle_activation_signal_cpp, muscle_activation_signal, config->MUSCLE_COUNT * sizeof( float ) );
}


/** Destructor
 *
 */
owOpenCLSolver::~owOpenCLSolver(void)
{
	destroy();
}
