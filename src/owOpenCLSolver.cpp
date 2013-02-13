#include <stdexcept>
#include <iostream>
#include <fstream>

#include "owOpenCLSolver.h"

const float xmin = XMIN;
const float xmax = XMAX;
const float ymin = YMIN;
const float ymax = YMAX;
const float zmin = ZMIN;
const float zmax = ZMAX;

int gridCellsX = (int)( ( XMAX - XMIN ) / h ) + 1;
int gridCellsY = (int)( ( YMAX - YMIN ) / h ) + 1;
int gridCellsZ = (int)( ( ZMAX - ZMIN ) / h ) + 1;
int gridCellCount = gridCellsX * gridCellsY * gridCellsZ;
int * _particleIndex = new int[ PARTICLE_COUNT * 2 ];
unsigned int * gridNextNonEmptyCellBuffer = new unsigned int[gridCellCount+1];
extern int numOfElasticConnections;
int myCompare( const void * v1, const void * v2 ); 

owOpenCLSolver::owOpenCLSolver(const float * positionBuffer, const float * velocityBuffer, const float * elasticConnections)
{
	try{
		initializeOpenCL();
		create_ocl_buffer( "acceleration", acceleration, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 * 2 ) );
		create_ocl_buffer( "gridCellIndex", gridCellIndex, CL_MEM_READ_WRITE, ( ( gridCellCount + 1 ) * sizeof( unsigned int ) * 1 ) );
		create_ocl_buffer( "gridCellIndexFixedUp", gridCellIndexFixedUp, CL_MEM_READ_WRITE, ( ( gridCellCount + 1 ) * sizeof( unsigned int ) * 1 ) );
		create_ocl_buffer( "neighborMap", neighborMap, CL_MEM_READ_WRITE, ( NK * sizeof( float ) * 2 ) );
		create_ocl_buffer( "particleIndex", particleIndex, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( unsigned int ) * 2 ) );
		create_ocl_buffer( "particleIndexBack", particleIndexBack, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( unsigned int ) ) );
		create_ocl_buffer( "position", position, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 ) );
		create_ocl_buffer( "pressure", pressure, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 ) );
		create_ocl_buffer( "rho", rho, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 2 ) );
		create_ocl_buffer( "rhoInv", rhoInv, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) ) );
		create_ocl_buffer( "sortedPosition", sortedPosition, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 * 2 ) );
		create_ocl_buffer( "sortedVelocity", sortedVelocity, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 ) );
		create_ocl_buffer( "velocity", velocity, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 ) );
		// Create kernels
		create_ocl_kernel("clearBuffers", clearBuffers);
		create_ocl_kernel("findNeighbors", findNeighbors);
		create_ocl_kernel("hashParticles", hashParticles);
		create_ocl_kernel("indexx", indexx);
		create_ocl_kernel("sortPostPass", sortPostPass);
		// Additional kernels PCISPH
		create_ocl_kernel("pcisph_computeForcesAndInitPressure", pcisph_computeForcesAndInitPressure);
		create_ocl_kernel("pcisph_integrate", pcisph_integrate);
		create_ocl_kernel("pcisph_predictPositions", pcisph_predictPositions);
		create_ocl_kernel("pcisph_predictDensity", pcisph_predictDensity);
		create_ocl_kernel("pcisph_correctPressure", pcisph_correctPressure);
		create_ocl_kernel("pcisph_computePressureForceAcceleration", pcisph_computePressureForceAcceleration);
		create_ocl_kernel("pcisph_computeDensity", pcisph_computeDensity);
		create_ocl_kernel("pcisph_computeElasticForces", pcisph_computeElasticForces);
		//Copy positionBuffer and velocityBuffer to the OpenCL Device
		copy_buffer_to_device( positionBuffer, position, PARTICLE_COUNT * sizeof( float ) * 4 );
		copy_buffer_to_device( velocityBuffer, velocity, PARTICLE_COUNT * sizeof( float ) * 4 );
		//elastic connections
		if(elasticConnections != NULL){
			create_ocl_buffer("elasticConnectionsData", elasticConnectionsData,CL_MEM_READ_WRITE,numOfElasticConnections * sizeof(float) * 4);
			copy_buffer_to_device(elasticConnections,elasticConnectionsData,numOfElasticConnections * sizeof(float) * 4);
		}
	}catch( std::exception &e ){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}
void owOpenCLSolver::initializeOpenCL()
{
	cl_int err;
	std::vector< cl::Platform > platformList;
	err = cl::Platform::get( &platformList );
	if( platformList.size() < 1 ){
		throw std::runtime_error( "No OpenCL platforms found" );
	}
	char cBuffer[1024];
	cl_platform_id clSelectedPlatformID = NULL;
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
	//0-CPU, 1-GPU// depends on order appropriet drivers was instaled
	int plList=1;//selected platform index in platformList array
	cl_context_properties cprops[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties) (platformList[plList])(), 0 };
	context = cl::Context( CL_DEVICE_TYPE_ALL, cprops, NULL, NULL, &err );
	devices = context.getInfo< CL_CONTEXT_DEVICES >();
	if( devices.size() < 1 ){
		throw std::runtime_error( "No OpenCL devices found" );
	}
	//Print some infrmation about chouse platform
	int value;
	cl_int result = devices[0].getInfo(CL_DEVICE_NAME,&cBuffer);// CL_INVALID_VALUE = -30;		
	if(result == CL_SUCCESS) printf("CL_PLATFORM_VERSION [%d]: CL_DEVICE_NAME [%d]: \t%s\n",plList, 0, cBuffer);
	result = devices[0].getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE,&value);
	if(result == CL_SUCCESS) printf("CL_PLATFORM_VERSION [%d]: CL_DEVICE_MAX_WORK_GROUP_SIZE [%d]: \t%d\n",plList, 0, value);
	result = devices[0].getInfo(CL_DEVICE_MAX_COMPUTE_UNITS,&value);
	if(result == CL_SUCCESS) printf("CL_PLATFORM_VERSION [%d]: CL_DEVICE_MAX_COMPUTE_UNITS [%d]: \t%d\n",plList, 0, value);
	result = devices[0].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE,&value);
	if(result == CL_SUCCESS) printf("CL_PLATFORM_VERSION [%d]: CL_DEVICE_GLOBAL_MEM_SIZE [%d]: \t%d\n",plList, 0, value);
	result = devices[0].getInfo(CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,&value);
	if(result == CL_SUCCESS) printf("CL_PLATFORM_VERSION [%d]: CL_DEVICE_GLOBAL_MEM_CACHE_SIZE [%d]: \t%d\n",plList, 0, value);

	queue = cl::CommandQueue( context, devices[ 0 ], 0, &err );
	if( err != CL_SUCCESS ){
		throw std::runtime_error( "failed to create command queue" );
	}

	std::string sourceFileName( OPENCL_PORGRAMM_PATH );
	std::ifstream file( sourceFileName.c_str() );
	if( !file.is_open() ){
		throw std::runtime_error( "could not open file " + sourceFileName );
	}
	std::string programSource( std::istreambuf_iterator<char>( file ), ( std::istreambuf_iterator<char>() ));
	cl::Program::Sources source( 1, std::make_pair( programSource.c_str(), programSource.length()+1 ));
	program = cl::Program( context, source );
#if INTEL_OPENCL_DEBUG
	err = program.build( devices, OPENCL_DEBUG_PORGRAMM_PATH );
#else
	err = program.build( devices, "" );
#endif
	if( err != CL_SUCCESS ){
		std::string compilationErrors;
		compilationErrors = program.getBuildInfo< CL_PROGRAM_BUILD_LOG >( devices[ 0 ] );
		std::cerr << "Compilation failed: " << std::endl << compilationErrors << std::endl;
		throw std::runtime_error( "failed to build program" );
	}
	return;
}
//Kernes function definition
unsigned int owOpenCLSolver::_runClearBuffers()
{
	// Stage ClearBuffers
	clearBuffers.setArg( 0, neighborMap );
	int err = queue.enqueueNDRangeKernel(clearBuffers, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_runHashParticles()
{
	// Stage HashParticles
	hashParticles.setArg( 0, position );
	hashParticles.setArg( 1, gridCellsX );
	hashParticles.setArg( 2, gridCellsY );
	hashParticles.setArg( 3, gridCellsZ );
	hashParticles.setArg( 4, hashGridCellSizeInv );
	hashParticles.setArg( 5, xmin );
	hashParticles.setArg( 6, ymin );
	hashParticles.setArg( 7, zmin );
	hashParticles.setArg( 8, particleIndex );
	int err = queue.enqueueNDRangeKernel(
		hashParticles, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}

unsigned int owOpenCLSolver::_runSort()
{
	copy_buffer_from_device( _particleIndex, particleIndex,PARTICLE_COUNT * 2 * sizeof( int ) );
	qsort( _particleIndex, PARTICLE_COUNT, 2 * sizeof( int ), myCompare );
	copy_buffer_to_device( _particleIndex, particleIndex, PARTICLE_COUNT * 2 * sizeof( int ) );
	return 0;
}
unsigned int owOpenCLSolver::_runSortPostPass()
{
	// Stage SortPostPass
	sortPostPass.setArg( 0, particleIndex );
	sortPostPass.setArg( 1, particleIndexBack );
	sortPostPass.setArg( 2, position );
	sortPostPass.setArg( 3, velocity );
	sortPostPass.setArg( 4, sortedPosition );
	sortPostPass.setArg( 5, sortedVelocity );
	int err = queue.enqueueNDRangeKernel(
		sortPostPass, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_runIndexx()
{
	// Stage Indexx
	indexx.setArg( 0, particleIndex );
	gridCellCount = ((gridCellsX) * (gridCellsY)) * (gridCellsZ);
	indexx.setArg( 1, gridCellCount );
	indexx.setArg( 2, gridCellIndex );
	int gridCellCountRoundedUp = ((( gridCellCount - 1 ) / 256 ) + 1 ) * 256;
	int err = queue.enqueueNDRangeKernel(
		indexx, cl::NullRange, cl::NDRange( (int) ( gridCellCountRoundedUp ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_runIndexPostPass()
{
	// Stage IndexPostPass
	//28aug_Palyanov_start_block
	copy_buffer_from_device( gridNextNonEmptyCellBuffer, gridCellIndex,(gridCellCount+1) * sizeof( unsigned int ) * 1 );
	int recentNonEmptyCell = gridCellCount;
	for(int i=gridCellCount;i>=0;i--)
	{
		if(gridNextNonEmptyCellBuffer[i]==NO_CELL_ID)
			gridNextNonEmptyCellBuffer[i] = recentNonEmptyCell; 
		else recentNonEmptyCell = gridNextNonEmptyCellBuffer[i];
	}
	int err = copy_buffer_to_device( gridNextNonEmptyCellBuffer,gridCellIndexFixedUp,(gridCellCount+1) * sizeof( unsigned int ) * 1 );
	return err;
}
unsigned int owOpenCLSolver::_runFindNeighbors()
{
	// Stage FindNeighbors
	findNeighbors.setArg( 0, gridCellIndexFixedUp );
	findNeighbors.setArg( 1, sortedPosition );
	gridCellCount = ((gridCellsX) * (gridCellsY)) * (gridCellsZ);
	findNeighbors.setArg( 2, gridCellCount );
	findNeighbors.setArg( 3, gridCellsX );
	findNeighbors.setArg( 4, gridCellsY );
	findNeighbors.setArg( 5, gridCellsZ );
	findNeighbors.setArg( 6, h );
	findNeighbors.setArg( 7, hashGridCellSize );
	findNeighbors.setArg( 8, hashGridCellSizeInv );
	findNeighbors.setArg( 9, simulationScale );
	findNeighbors.setArg( 10, xmin );
	findNeighbors.setArg( 11, ymin );
	findNeighbors.setArg( 12, zmin );
	findNeighbors.setArg( 13, neighborMap );
	int err = queue.enqueueNDRangeKernel(
		findNeighbors, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_computeDensity()
{
	// Stage ComputeDensityPressure
	pcisph_computeDensity.setArg( 0, neighborMap );
	pcisph_computeDensity.setArg( 1, Wpoly6Coefficient );
	pcisph_computeDensity.setArg( 2, gradWspikyCoefficient );
	pcisph_computeDensity.setArg( 3, h );
	pcisph_computeDensity.setArg( 4, mass );
	pcisph_computeDensity.setArg( 5, rho0 );
	pcisph_computeDensity.setArg( 6, simulationScale );
	pcisph_computeDensity.setArg( 7, stiffness );
	pcisph_computeDensity.setArg( 8, sortedPosition );
	pcisph_computeDensity.setArg( 9, pressure );
	pcisph_computeDensity.setArg(10, rho );
	pcisph_computeDensity.setArg(11, rhoInv );
	pcisph_computeDensity.setArg(12, particleIndexBack );
	pcisph_computeDensity.setArg(13, delta );
	int err = queue.enqueueNDRangeKernel(
		pcisph_computeDensity, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_computeForcesAndInitPressure()
{
	pcisph_computeForcesAndInitPressure.setArg( 0, neighborMap );
	pcisph_computeForcesAndInitPressure.setArg( 1, rho );
	pcisph_computeForcesAndInitPressure.setArg( 2, pressure );
	pcisph_computeForcesAndInitPressure.setArg( 3, sortedPosition );
	pcisph_computeForcesAndInitPressure.setArg( 4, sortedVelocity );
	pcisph_computeForcesAndInitPressure.setArg( 5, acceleration );
	pcisph_computeForcesAndInitPressure.setArg( 6, particleIndexBack );
	pcisph_computeForcesAndInitPressure.setArg( 7, Wpoly6Coefficient );
	pcisph_computeForcesAndInitPressure.setArg( 8, del2WviscosityCoefficient );
	pcisph_computeForcesAndInitPressure.setArg( 9, h );
	pcisph_computeForcesAndInitPressure.setArg(10, mass );
	pcisph_computeForcesAndInitPressure.setArg(11, mu );
	pcisph_computeForcesAndInitPressure.setArg(12, simulationScale );
	pcisph_computeForcesAndInitPressure.setArg(13, gravity_x );
	pcisph_computeForcesAndInitPressure.setArg(14, gravity_y );
	pcisph_computeForcesAndInitPressure.setArg(15, gravity_z );
	int err = queue.enqueueNDRangeKernel(
		pcisph_computeForcesAndInitPressure, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_computeElasticForces()
{
	if(numOfElasticConnections == 0 )
		return 0;
	pcisph_computeElasticForces.setArg( 0, neighborMap );
	pcisph_computeElasticForces.setArg( 1, sortedPosition );
	pcisph_computeElasticForces.setArg( 2, sortedVelocity );
	pcisph_computeElasticForces.setArg( 3, acceleration );
	pcisph_computeElasticForces.setArg( 4, particleIndexBack );
	pcisph_computeElasticForces.setArg( 5, h );
	pcisph_computeElasticForces.setArg( 6, mass );
	pcisph_computeElasticForces.setArg( 7, simulationScale );
	pcisph_computeElasticForces.setArg( 8, numOfElasticConnections );
	pcisph_computeElasticForces.setArg( 9, elasticConnectionsData );
	int elasticConnectionsCountRoundedUp = ((( numOfElasticConnections - 1 ) / 256 ) + 1 ) * 256;
	int err = queue.enqueueNDRangeKernel(
		pcisph_computeElasticForces, cl::NullRange, cl::NDRange( (int) ( elasticConnectionsCountRoundedUp ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_predictPositions()
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
	pcisph_predictPositions.setArg( 10, xmin );
	pcisph_predictPositions.setArg( 11, xmax );
	pcisph_predictPositions.setArg( 12, ymin );
	pcisph_predictPositions.setArg( 13, ymax );
	pcisph_predictPositions.setArg( 14, zmin );
	pcisph_predictPositions.setArg( 15, zmax );
	pcisph_predictPositions.setArg( 16, damping );
	pcisph_predictPositions.setArg( 17, position );
	pcisph_predictPositions.setArg( 18, velocity );
	pcisph_predictPositions.setArg( 19, r0 );
	pcisph_predictPositions.setArg( 20, neighborMap );
	int err = queue.enqueueNDRangeKernel(
		pcisph_predictPositions, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_predictDensity()
{
	// Stage ComputeDensityPressure
	pcisph_predictDensity.setArg( 0, neighborMap );
	pcisph_predictDensity.setArg( 1, particleIndexBack );
	pcisph_predictDensity.setArg( 2, Wpoly6Coefficient );
	pcisph_predictDensity.setArg( 3, gradWspikyCoefficient );
	pcisph_predictDensity.setArg( 4, h );
	pcisph_predictDensity.setArg( 5, mass );
	pcisph_predictDensity.setArg( 6, rho0 );
	pcisph_predictDensity.setArg( 7, simulationScale );
	pcisph_predictDensity.setArg( 8, stiffness );
	pcisph_predictDensity.setArg( 9, sortedPosition );
	pcisph_predictDensity.setArg(10, pressure );
	pcisph_predictDensity.setArg(11, rho );
	pcisph_predictDensity.setArg(12, rhoInv );
	pcisph_predictDensity.setArg(13, delta );
	int err = queue.enqueueNDRangeKernel(
		pcisph_predictDensity, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_correctPressure()
{
	// Stage ComputeDensityPressure
	pcisph_correctPressure.setArg( 0, neighborMap );
	pcisph_correctPressure.setArg( 1, particleIndexBack );
	pcisph_correctPressure.setArg( 2, Wpoly6Coefficient );
	pcisph_correctPressure.setArg( 3, gradWspikyCoefficient );
	pcisph_correctPressure.setArg( 4, h );
	pcisph_correctPressure.setArg( 5, mass );
	pcisph_correctPressure.setArg( 6, rho0 );
	pcisph_correctPressure.setArg( 7, simulationScale );
	pcisph_correctPressure.setArg( 8, stiffness );
	pcisph_correctPressure.setArg( 9, sortedPosition );
	pcisph_correctPressure.setArg(10, pressure );
	pcisph_correctPressure.setArg(11, rho );
	pcisph_correctPressure.setArg(12, rhoInv );
	pcisph_correctPressure.setArg(13, delta );
	pcisph_correctPressure.setArg(14, position );
	pcisph_correctPressure.setArg(15, particleIndex );
	int err = queue.enqueueNDRangeKernel(
		pcisph_correctPressure, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_computePressureForceAcceleration()
{
	// Stage ComputeAcceleration
	pcisph_computePressureForceAcceleration.setArg( 0, neighborMap );
	pcisph_computePressureForceAcceleration.setArg( 1, pressure );
	pcisph_computePressureForceAcceleration.setArg( 2, rho );
	pcisph_computePressureForceAcceleration.setArg( 3, rhoInv );
	pcisph_computePressureForceAcceleration.setArg( 4, sortedPosition );
	pcisph_computePressureForceAcceleration.setArg( 5, sortedVelocity );
	pcisph_computePressureForceAcceleration.setArg( 6, particleIndexBack );
	pcisph_computePressureForceAcceleration.setArg( 7, CFLLimit );
	pcisph_computePressureForceAcceleration.setArg( 8, del2WviscosityCoefficient );
	pcisph_computePressureForceAcceleration.setArg( 9, gradWspikyCoefficient );
	pcisph_computePressureForceAcceleration.setArg( 10, h );
	pcisph_computePressureForceAcceleration.setArg( 11, mass );
	pcisph_computePressureForceAcceleration.setArg( 12, mu );
	pcisph_computePressureForceAcceleration.setArg( 13, simulationScale );
	pcisph_computePressureForceAcceleration.setArg( 14, acceleration );
	pcisph_computePressureForceAcceleration.setArg( 15, rho0 );
	pcisph_computePressureForceAcceleration.setArg( 16, position );
	pcisph_computePressureForceAcceleration.setArg( 17, particleIndex );
	int err = queue.enqueueNDRangeKernel(
		pcisph_computePressureForceAcceleration, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#if QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return err;
}
unsigned int owOpenCLSolver::_run_pcisph_integrate()
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
	pcisph_integrate.setArg( 10, xmin );
	pcisph_integrate.setArg( 11, xmax );
	pcisph_integrate.setArg( 12, ymin );
	pcisph_integrate.setArg( 13, ymax );
	pcisph_integrate.setArg( 14, zmin );
	pcisph_integrate.setArg( 15, zmax );
	pcisph_integrate.setArg( 16, damping );
	pcisph_integrate.setArg( 17, position );
	pcisph_integrate.setArg( 18, velocity );
	pcisph_integrate.setArg( 19, rho );
	pcisph_integrate.setArg( 20, r0 );
	pcisph_integrate.setArg( 21, neighborMap );
	int err = queue.enqueueNDRangeKernel(
		pcisph_integrate, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
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
owOpenCLSolver::~owOpenCLSolver(void)
{
	delete [] gridNextNonEmptyCellBuffer;
	delete [] _particleIndex;
}
