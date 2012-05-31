
#pragma comment( lib, "opencl.lib" )							// Подключается библиотека opencl.lib

//#define QUEUE_EACH_KERNEL//debugging feature

#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <math.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <windows.h>
#include <time.h>

#include "sph.h"

//AP2012//#include "dx10_render.h"
/*
#ifdef NDEBUG
#define ENABLE_OPENCL_RADIXSORT
#include "radixsort.hpp"
#endif
*/

//see opencl.hpp to find this:
/*! \class BufferGL
 * \brief Memory buffer interface for GL interop.
 */

//AP2012//#define USE_DX_INTEROP
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#include <OpenCL/cl_d3d10.h>
#else
#include <CL/cl.hpp>
//AP2012//#include <CL/cl_d3d10.h>
#endif

const float xmin = XMIN;
const float xmax = XMAX;
const float ymin = YMIN;
const float ymax = YMAX;
const float zmin = ZMIN;
const float zmax = ZMAX;

const float rho0 = 1000.0f;
const float stiffness = 0.75f;
const float h = 3.34f;
const float hashGridCellSize = 2.0f * h;
const float hashGridCellSizeInv = 1.0f / hashGridCellSize;
const float mass = 0.001f;
const float simulationScale = 0.004f;
const float simulationScaleInv = 1.0f / simulationScale;
const float mu = 10.0f;
const float timeStep = 0.0042f;
const float CFLLimit = 100.0f;
const int NK = NEIGHBOR_COUNT * PARTICLE_COUNT;
const float damping = 0.75f;

const float Wpoly6Coefficient = 315.0f / ( 64.0f * M_PI * pow( h * simulationScale, 9.0f ) );
const float gradWspikyCoefficient= -45.0f / ( M_PI * pow( h * simulationScale, 6.0f ) );
const float del2WviscosityCoefficient = -gradWspikyCoefficient;

const float gravity_x = 0.0f;
const float gravity_y = -9.8f;
const float gravity_z = 0.0f;

//int particleCount = PARTICLE_COUNT;//AP2012

int gridCellsX;
int gridCellsY;
int gridCellsZ;
int gridCellCount;


float * positionBuffer;
float * velocityBuffer;

cl::Context context;
std::vector< cl::Device > devices;
cl::CommandQueue queue;
cl::Program program;

// Buffers
cl::Buffer acceleration;
cl::Buffer gridCellIndex;
cl::Buffer gridCellIndexFixedUp;
cl::Buffer neighborMap;
cl::Buffer particleIndex;
cl::Buffer position;
cl::Buffer pressure;
cl::Buffer rho;
cl::Buffer rhoInv;
cl::Buffer sortedPosition;
cl::Buffer sortedVelocity;
cl::Buffer velocity;

// Kernels
cl::Kernel clearBuffers;
cl::Kernel computeAcceleration;
cl::Kernel computeDensityPressure;
cl::Kernel findNeighbors;
cl::Kernel hashParticles;
cl::Kernel indexPostPass;
cl::Kernel indexx;
cl::Kernel integrate;
cl::Kernel sortPostPass;

/*//AP2012
#ifdef NDEBUG
amd::RadixSortCL radixSort;
#endif
*///AP2012

/*
#if !defined( NDEBUG )
bool testBeginClearBuffers( cl::CommandQueue queue );
bool testEndClearBuffers( cl::CommandQueue queue, int _NK, cl::Buffer neighborMap );
bool testBeginComputeAcceleration( cl::CommandQueue queue );
bool testEndComputeAcceleration( cl::CommandQueue queue, int _NK,
								float simulationScale, float gradWspikyCoefficient,
								float del2WviscosityCoefficient, float CFLLimit,
								cl::Buffer neighborMap, cl::Buffer pressure, cl::Buffer rho,
								cl::Buffer rhoInv, cl::Buffer sortedPosition, cl::Buffer sortedVelocity,
								cl::Buffer acceleration,
								float gravity_x, float gravity_y, float gravity_z,
								float mass, float h, float mu );
bool testBeginComputeDensityPressure( cl::CommandQueue queue );
bool testEndComputeDensityPressure( cl::CommandQueue queue, int _NK, float mass, float h,
								   float Wpoly6Coefficient, float rho0, float stiffness,
								   cl::Buffer neighborMap, cl::Buffer pressure, cl::Buffer rho,
								   cl::Buffer rhoInv, float simulationScale );
bool testBeginFindNeighbors( cl::CommandQueue queue );
bool testEndFindNeighbors( cl::CommandQueue queue, int gridCellCount, float hashGridCellSize, float hashGridCellSizeInv,
						  cl::Buffer gridCellIndexFixedUp, cl::Buffer sortedPosition, cl::Buffer neighborMap, int _NK,
						  float simulationScale, float h );
bool testBeginHashParticles( cl::CommandQueue queue );
bool testEndHashParticles( cl::CommandQueue queue, float hashGridCellSize, int gridCellCount,
						  int gridCellsX, int gridCellsY, int gridCellsZ,
						  float xmin, float ymin, float zmin,
						  cl::Buffer position, cl::Buffer particleIndex );
bool testBeginIndexPostPass( cl::CommandQueue queue );
bool testEndIndexPostPass( cl::CommandQueue queue, int gridCellCount, cl::Buffer gridCellIndex,
						  cl::Buffer gridCellIndexFixedUp, cl::Buffer particleIndex );
bool testBeginIndexx( cl::CommandQueue queue );
bool testEndIndexx( cl::CommandQueue queue, int gridCellCount, cl::Buffer particleIndex,
				  cl::Buffer gridCellIndex );
bool testBeginIntegrate( cl::CommandQueue queue );
bool testEndIntegrate( cl::CommandQueue queue, cl::Buffer acceleration, cl::Buffer sortedPosition,
					  cl::Buffer sortedVelocity, cl::Buffer position, cl::Buffer velocity,
					  float xmin, float ymin, float zmin, float xmax, float ymax, float zmax );
bool testBeginSort( cl::CommandQueue queue, cl::Buffer particleIndex, int gridCellCount );
bool testEndSort( cl::CommandQueue queue, cl::Buffer particleIndex );
bool testBeginSortPostPass( cl::CommandQueue queue );
bool testEndSortPostPass( cl::CommandQueue queue, cl::Buffer particleIndex, cl::Buffer position,
						 cl::Buffer velocity, cl::Buffer sortedPosition, cl::Buffer sortedVelocity,
						 int gridCellCount );
#endif//NDEBUG
*/

/*
Code:
for (int i=0; i < n; i++)
{
     //your code
}


OpenCL code:
Code:
kernel ( /your arguments/ )
{
    i = get_global_id(0);
   //your code
}
EnqueueNDRange(yourkernel, work_items, etc.)

work_dim = {n, p}
Regular Code:

Code:
for (int i=0; i < n; i++)
{
   for (int j=0; j < p; j++)
   {
     //your code
   }
}


OpenCL code:
Code:
kernel ( /your arguments/ )
{
    i = get_global_id(0);
    j = get_global_id(1);
   //your code
}
EnqueueNDRange(yourkernel, work_items, etc.)
*/


//sphFluidDemo.cl: __kernel void clearBuffers( __global float2 * neighborMap )


unsigned int
_runClearBuffers( cl::CommandQueue queue ){
	// Stage ClearBuffers
	clearBuffers.setArg( 0, neighborMap );
	queue.enqueueNDRangeKernel(
		clearBuffers, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#ifdef QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return 0;
}

/*
#if !defined( NDEBUG )

unsigned int
_testClearBuffers( cl::CommandQueue queue ){
	queue.finish();
	// Invoke stage ClearBuffers
	bool thisTestPassed = true;
	std::cerr << "Testing stage ClearBuffers:" << std::endl;
	thisTestPassed &= testBeginClearBuffers( queue );
	unsigned int result = _runClearBuffers( queue );
	queue.finish();
	thisTestPassed &= testEndClearBuffers( queue, NK, neighborMap );

	if( !thisTestPassed ){
		std::cerr << "Test FAILED in stage ClearBuffers" << std::endl;
	}else{
		std::cerr << "Test SUCCEEDED in stage ClearBuffers" << std::endl;
	}
	return result;
}

#endif
*/

unsigned int
_runComputeAcceleration( cl::CommandQueue queue ){
	// Stage ComputeAcceleration
	computeAcceleration.setArg( 0, neighborMap );
	computeAcceleration.setArg( 1, pressure );
	computeAcceleration.setArg( 2, rho );
	computeAcceleration.setArg( 3, rhoInv );
	computeAcceleration.setArg( 4, sortedPosition );
	computeAcceleration.setArg( 5, sortedVelocity );
	computeAcceleration.setArg( 6, CFLLimit );
	computeAcceleration.setArg( 7, del2WviscosityCoefficient );
	computeAcceleration.setArg( 8, gradWspikyCoefficient );
	computeAcceleration.setArg( 9, h );
	computeAcceleration.setArg( 10, mass );
	computeAcceleration.setArg( 11, mu );
	computeAcceleration.setArg( 12, simulationScale );
	computeAcceleration.setArg( 13, acceleration );
	queue.enqueueNDRangeKernel(
		computeAcceleration, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#ifdef QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return 0;
}

/*
#if !defined( NDEBUG )

unsigned int
_testComputeAcceleration( cl::CommandQueue queue ){
	queue.finish();
	// Invoke stage ComputeAcceleration
	bool thisTestPassed = true;
	std::cerr << "Testing stage ComputeAcceleration:" << std::endl;
	thisTestPassed &= testBeginComputeAcceleration( queue );
	unsigned int result = _runComputeAcceleration( queue );
	queue.finish();
	thisTestPassed &= testEndComputeAcceleration( queue, NK, 
		simulationScale, gradWspikyCoefficient, del2WviscosityCoefficient, CFLLimit,
		neighborMap, pressure, rho, rhoInv, sortedPosition, sortedVelocity, acceleration,
		gravity_x, gravity_y, gravity_z, mass, h, mu );

	if( !thisTestPassed ){
		std::cerr << "Test FAILED in stage ComputeAcceleration" << std::endl;
	}else{
		std::cerr << "Test SUCCEEDED in stage ComputeAcceleration" << std::endl;
	}
	return result;
}


#endif
*/

unsigned int
_runComputeDensityPressure( cl::CommandQueue queue ){
	// Stage ComputeDensityPressure
	computeDensityPressure.setArg( 0, neighborMap );
	computeDensityPressure.setArg( 1, Wpoly6Coefficient );
	computeDensityPressure.setArg( 2, h );
	computeDensityPressure.setArg( 3, mass );
	computeDensityPressure.setArg( 4, rho0 );
	computeDensityPressure.setArg( 5, simulationScale );
	computeDensityPressure.setArg( 6, stiffness );
	computeDensityPressure.setArg( 7, pressure );
	computeDensityPressure.setArg( 8, rho );
	computeDensityPressure.setArg( 9, rhoInv );
	queue.enqueueNDRangeKernel(
		computeDensityPressure, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#ifdef QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return 0;
}
/*
#if !defined( NDEBUG )

unsigned int
_testComputeDensityPressure( cl::CommandQueue queue ){
	queue.finish();
	// Invoke stage ComputeDensityPressure
	bool thisTestPassed = true;
	std::cerr << "Testing stage ComputeDensityPressure:" << std::endl;
	thisTestPassed &= testBeginComputeDensityPressure( queue );
	unsigned int result = _runComputeDensityPressure( queue );
	queue.finish();
	thisTestPassed &= testEndComputeDensityPressure( queue, NK, mass, h, Wpoly6Coefficient,
		rho0, stiffness, neighborMap, pressure, rho, rhoInv, simulationScale );

	if( !thisTestPassed ){
		std::cerr << "Test FAILED in stage ComputeDensityPressure" << std::endl;
	}else{
		std::cerr << "Test SUCCEEDED in stage ComputeDensityPressure" << std::endl;
	}
	return result;
}

#endif
*/

unsigned int
_runFindNeighbors( cl::CommandQueue queue ){
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
	queue.enqueueNDRangeKernel(
		findNeighbors, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#ifdef QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return 0;
}

/*
#if !defined( NDEBUG )

unsigned int
_testFindNeighbors( cl::CommandQueue queue ){
	queue.finish();
	// Invoke stage FindNeighbors
	bool thisTestPassed = true;
	std::cerr << "Testing stage FindNeighbors:" << std::endl;
	thisTestPassed &= testBeginFindNeighbors( queue );
	unsigned int result = _runFindNeighbors( queue );
	queue.finish();
	thisTestPassed &= testEndFindNeighbors( queue, gridCellCount, hashGridCellSize,
		hashGridCellSizeInv, gridCellIndexFixedUp, sortedPosition, neighborMap, NK,
		simulationScale, h );

	if( !thisTestPassed ){
		std::cerr << "Test FAILED in stage FindNeighbors" << std::endl;
	}else{
		std::cerr << "Test SUCCEEDED in stage FindNeighbors" << std::endl;
	}
	return result;
}

#endif
*/



unsigned int
_runHashParticles( cl::CommandQueue queue ){
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
	queue.enqueueNDRangeKernel(
		hashParticles, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#ifdef QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return 0;
}

/*
#if !defined( NDEBUG )

unsigned int
_testHashParticles( cl::CommandQueue queue ){
	queue.finish();
	// Invoke stage HashParticles
	bool thisTestPassed = true;
	std::cerr << "Testing stage HashParticles:" << std::endl;
	thisTestPassed &= testBeginHashParticles( queue );
	unsigned int result = _runHashParticles( queue );
	queue.finish();
	thisTestPassed &= testEndHashParticles( queue, hashGridCellSize, gridCellCount,
		gridCellsX, gridCellsY, gridCellsZ, xmin, ymin, zmin, position, particleIndex );

	if( !thisTestPassed ){
		std::cerr << "Test FAILED in stage HashParticles" << std::endl;
	}else{
		std::cerr << "Test SUCCEEDED in stage HashParticles" << std::endl;
	}
	return result;
}

#endif
*/


unsigned int
_runIndexPostPass( cl::CommandQueue queue ){
	// Stage IndexPostPass
	indexPostPass.setArg( 0, gridCellIndex );
	gridCellCount = ((gridCellsX) * (gridCellsY)) * (gridCellsZ);
	indexPostPass.setArg( 1, gridCellCount );
	indexPostPass.setArg( 2, gridCellIndexFixedUp );
	int gridCellCountRoundedUp = ((( gridCellCount - 1 ) / 256 ) + 1 ) * 256;
	queue.enqueueNDRangeKernel(
		indexPostPass, cl::NullRange, cl::NDRange( (int) ( gridCellCountRoundedUp ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#ifdef QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return 0;
}

/*
#if !defined( NDEBUG )

unsigned int
_testIndexPostPass( cl::CommandQueue queue ){
	queue.finish();
	// Invoke stage IndexPostPass
	bool thisTestPassed = true;
	std::cerr << "Testing stage IndexPostPass:" << std::endl;
	thisTestPassed &= testBeginIndexPostPass( queue );
	unsigned int result = _runIndexPostPass( queue );
	queue.finish();
	thisTestPassed &= testEndIndexPostPass( queue, gridCellCount, gridCellIndex,
		gridCellIndexFixedUp, particleIndex );

	if( !thisTestPassed ){
		std::cerr << "Test FAILED in stage IndexPostPass" << std::endl;
	}else{
		std::cerr << "Test SUCCEEDED in stage IndexPostPass" << std::endl;
	}
	return result;
}

#endif
*/

unsigned int
_runIndexx( cl::CommandQueue queue ){
	// Stage Indexx
	indexx.setArg( 0, particleIndex );
	gridCellCount = ((gridCellsX) * (gridCellsY)) * (gridCellsZ);
	indexx.setArg( 1, gridCellCount );
	indexx.setArg( 2, gridCellIndex );
	int gridCellCountRoundedUp = ((( gridCellCount - 1 ) / 256 ) + 1 ) * 256;
	queue.enqueueNDRangeKernel(
		indexx, cl::NullRange, cl::NDRange( (int) ( gridCellCountRoundedUp ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#ifdef QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return 0;
}

/*
#if !defined( NDEBUG )

unsigned int
_testIndexx( cl::CommandQueue queue ){
	queue.finish();
	// Invoke stage Indexx
	bool thisTestPassed = true;
	std::cerr << "Testing stage Indexx:" << std::endl;
	thisTestPassed &= testBeginIndexx( queue );
	unsigned int result = _runIndexx( queue );
	queue.finish();
	thisTestPassed &= testEndIndexx( queue, gridCellCount, particleIndex, gridCellIndex );

	if( !thisTestPassed ){
		std::cerr << "Test FAILED in stage Indexx" << std::endl;
	}else{
		std::cerr << "Test SUCCEEDED in stage Indexx" << std::endl;
	}
	return result;
}

#endif
*/


unsigned int
_runIntegrate( cl::CommandQueue queue ){
	// Stage Integrate
	integrate.setArg( 0, acceleration );
	integrate.setArg( 1, sortedPosition );
	integrate.setArg( 2, sortedVelocity );
	integrate.setArg( 3, gravity_x );
	integrate.setArg( 4, gravity_y );
	integrate.setArg( 5, gravity_z );
	integrate.setArg( 6, simulationScaleInv );
	integrate.setArg( 7, timeStep );
	integrate.setArg( 8, xmin );
	integrate.setArg( 9, xmax );
	integrate.setArg( 10, ymin );
	integrate.setArg( 11, ymax );
	integrate.setArg( 12, zmin );
	integrate.setArg( 13, zmax );
	integrate.setArg( 14, damping );
	integrate.setArg( 15, position );
	integrate.setArg( 16, velocity );
	queue.enqueueNDRangeKernel(
		integrate, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#ifdef QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return 0;
}

/*
#if !defined( NDEBUG )

unsigned int
_testIntegrate( cl::CommandQueue queue ){
	queue.finish();
	// Invoke stage Integrate
	bool thisTestPassed = true;
	std::cerr << "Testing stage Integrate:" << std::endl;
	thisTestPassed &= testBeginIntegrate( queue );
	unsigned int result = _runIntegrate( queue );
	queue.finish();
	thisTestPassed &= testEndIntegrate( queue, acceleration, sortedPosition, sortedVelocity,
		position, velocity, xmin, ymin, zmin, xmax, ymax, zmax );

	if( !thisTestPassed ){
		std::cerr << "Test FAILED in stage Integrate" << std::endl;
	}else{
		std::cerr << "Test SUCCEEDED in stage Integrate" << std::endl;
	}
	return result;
}

#endif
*/

//#ifndef NDEBUG
int myCompare( const void * v1, const void * v2 ){
	int * f1 = (int *)v1;
	int * f2 = (int *)v2;
	if( f1[ 0 ] < f2[ 0 ] ) return -1;
	if( f1[ 0 ] > f2[ 0 ] ) return +1;
	return 0;
}
//#endif

int * _particleIndex = new int[ PARTICLE_COUNT * 2 ];

unsigned int
_runSort( cl::CommandQueue queue ){
/*#ifdef NDEBUG
	radixSort.sort( particleIndex );
	radixSort.wait();
#else*/
	
	queue.enqueueReadBuffer( particleIndex, CL_TRUE, 0, PARTICLE_COUNT * 2 * sizeof( int ), _particleIndex );
	queue.finish();
	qsort( _particleIndex, PARTICLE_COUNT, 2 * sizeof( int ), myCompare );
	queue.enqueueWriteBuffer( particleIndex, CL_TRUE, 0, PARTICLE_COUNT * 2 * sizeof( int ), _particleIndex );
	queue.finish();	
	//delete [] _particleIndex;
//#endif
	return 0;
}

/*
#if !defined( NDEBUG )

unsigned int
_testSort( cl::CommandQueue queue ){
	// Invoke stage Sort
	bool thisTestPassed = true;
	std::cerr << "Testing stage Sort:" << std::endl;
	thisTestPassed &= testBeginSort( queue, particleIndex, gridCellCount );
	unsigned int result = _runSort( queue );
	thisTestPassed &= testEndSort( queue, particleIndex );
	if( !thisTestPassed ){
		std::cerr << "Test FAILED in stage Sort" << std::endl;
	}else{
		std::cerr << "Test SUCCEEDED in stage Sort" << std::endl;
	}
	return result;
}

#endif
*/

unsigned int
_runSortPostPass( cl::CommandQueue queue ){
	// Stage SortPostPass
	sortPostPass.setArg( 0, particleIndex );
	sortPostPass.setArg( 1, position );
	sortPostPass.setArg( 2, velocity );
	sortPostPass.setArg( 3, sortedPosition );
	sortPostPass.setArg( 4, sortedVelocity );
	queue.enqueueNDRangeKernel(
		sortPostPass, cl::NullRange, cl::NDRange( (int) ( PARTICLE_COUNT ) ),
#if defined( __APPLE__ )
		cl::NullRange, NULL, NULL );
#else
		cl::NDRange( (int)( 256 ) ), NULL, NULL );
#endif
#ifdef QUEUE_EACH_KERNEL
	queue.finish();
#endif
	return 0;
}

/*
#if !defined( NDEBUG )

unsigned int
_testSortPostPass( cl::CommandQueue queue ){
	queue.finish();
	// Invoke stage SortPostPass
	bool thisTestPassed = true;
	std::cerr << "Testing stage SortPostPass:" << std::endl;
	thisTestPassed &= testBeginSortPostPass( queue );
	unsigned int result = _runSortPostPass( queue );
	queue.finish();
	thisTestPassed &= testEndSortPostPass( queue, particleIndex, position, velocity,
		sortedPosition, sortedVelocity, gridCellCount );

	if( !thisTestPassed ){
		std::cerr << "Test FAILED in stage SortPostPass" << std::endl;
	}else{
		std::cerr << "Test SUCCEEDED in stage SortPostPass" << std::endl;
	}
	return result;
}

#endif
*/



void step()
{
	LARGE_INTEGER frequency;				// ticks per second
    LARGE_INTEGER t0, t1, t2;				// ticks
    double elapsedTime;
	int err;

    QueryPerformanceFrequency(&frequency);	// get ticks per second
    QueryPerformanceCounter(&t1);			// start timer
	t0 = t1;

//#ifdef NDEBUG
	printf("\n");

	_runClearBuffers( queue ); queue.finish();
	QueryPerformanceCounter(&t2);
	printf("_runClearBuffers: \t%9.3f ms\n",	(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart); t1 = t2;

	_runHashParticles( queue ); queue.finish();
	QueryPerformanceCounter(&t2);
	printf("_runHashParticles: \t%9.3f ms\n",	(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart); t1 = t2;

	_runSort( queue ); queue.finish();// <--here 
	QueryPerformanceCounter(&t2);
	printf("_runSort: \t%9.3f ms\n",	(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart); t1 = t2;

	_runSortPostPass( queue ); queue.finish();
	QueryPerformanceCounter(&t2);
	printf("_runSortPostPass: \t%9.3f ms\n",	(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart); t1 = t2;

	_runIndexx( queue ); queue.finish();
	QueryPerformanceCounter(&t2);
	printf("_runIndexx: \t%9.3f ms\n",	(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart); t1 = t2;

	_runIndexPostPass( queue ); queue.finish();
	QueryPerformanceCounter(&t2);
	printf("_runIndexPostPass: \t%9.3f ms\n",	(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart); t1 = t2;

	_runFindNeighbors( queue ); queue.finish();
	QueryPerformanceCounter(&t2);
	printf("_runFindNeighbors: \t%9.3f ms\n",	(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart); t1 = t2;

    _runComputeDensityPressure( queue ); queue.finish();
	QueryPerformanceCounter(&t2);
	printf("_runComputeDensityPressure: \t%9.3f ms\n",	(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart); t1 = t2;

    _runComputeAcceleration( queue ); queue.finish();
	QueryPerformanceCounter(&t2);
	printf("_runComputeAcceleration: \t%9.3f ms\n",	(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart); t1 = t2;

	_runIntegrate( queue ); queue.finish();
	QueryPerformanceCounter(&t2);
	printf("_runIntegrate: \t%9.3f ms\n",	(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart); t1 = t2;

/*
#else

	_testClearBuffers( queue );
	_testHashParticles( queue );
	_testSort( queue );
	_testSortPostPass( queue );
	_testIndexx( queue );
	_testIndexPostPass( queue );
	_testFindNeighbors( queue );
	_testComputeDensityPressure( queue );
	_testComputeAcceleration( queue );
	_testIntegrate( queue );

#endif
*/

	//printf("enter <queue.enqueueReadBuffer>, line 700 at main.cpp\n");
	err = queue.enqueueReadBuffer( position, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, positionBuffer );
	if( err != CL_SUCCESS ){
		throw std::runtime_error( "could not enqueue position read" );
	}
	queue.finish();
	QueryPerformanceCounter(&t2);// stop timer
	elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
	printf("_readBuffer: \t\t\t%9.3f ms\n",elapsedTime);


	elapsedTime = (t2.QuadPart - t0.QuadPart) * 1000.0 / frequency.QuadPart;
	printf("===============================\n");
	printf("_Total_step_time:\t\t%9.3f ms\n",elapsedTime);
	printf("");
}


void initializeOpenCL(
					  cl::Context & context,
					  std::vector< cl::Device > & devices,
					  cl::CommandQueue & queue,
					  cl::Program & program
					  )
{
	cl_int err;
	std::vector< cl::Platform > platformList;
	err = cl::Platform::get( &platformList );
	if( platformList.size() < 1 ){
		throw std::runtime_error( "no OpenCL platforms found" );
	}

	
///////////////////AP2012///////////////

	char cBuffer[1024];
	cl_platform_id clSelectedPlatformID = NULL;
	cl_platform_id cl_pl_id[10];
	cl_uint n_pl;
	clGetPlatformIDs(10,cl_pl_id,&n_pl);
		
	cl_int ciErrNum;// = oclGetPlatformID (&clSelectedPlatformID);
	//oclCheckError(ciErrNum, CL_SUCCESS);
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

///////////////////AP2012////////////////
/*
	cl_context_properties *cprops;
	cprops = new cl_context_properties[ 6 ];
	cprops[ 0 ] = CL_CONTEXT_D3D10_DEVICE_KHR; 
	cprops[ 1 ] = (intptr_t) DXUTGetD3D10Device();
	cprops[ 2 ] = CL_CONTEXT_PLATFORM;
	cprops[ 3 ] = (cl_context_properties)(platformList[0])();
	cprops[ 4 ] = cprops[ 5 ] = 0;

#ifdef NDEBUG
	context = cl::Context( CL_DEVICE_TYPE_ALL, cprops, NULL, NULL, &err );
#else
	context = cl::Context( CL_DEVICE_TYPE_ALL, cprops, NULL, NULL, &err );
#endif
*/

	//0-CPU, 1-GPU// depends on order appropriet drivers was instaled
	int plList=1;//selected platform index in platformList array

	cl_context_properties cprops[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties) (platformList[plList])(), 0 };

	//cl_context context1 = clCreateContext(cprops, CL_DEVICE_TYPE_CPU, NULL, NULL, NULL,&err);//( cl::Context( CL_DEVICE_TYPE_ALL, cprops, NULL, NULL, &err );
	context = cl::Context( CL_DEVICE_TYPE_ALL, cprops, NULL, NULL, &err );
	
	devices = context.getInfo< CL_CONTEXT_DEVICES >();
	if( devices.size() < 1 ){
		throw std::runtime_error( "no OpenCL devices found" );
	}

	///////////////////AP2012////////////////
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
	///////////////////AP2012////////////////

	queue = cl::CommandQueue( context, devices[ 0 ], 0, &err );
	if( err != CL_SUCCESS ){
		throw std::runtime_error( "failed to create command queue" );
	}

	std::string sourceFileName( "sphFluidDemo.cl" );
	std::ifstream file( sourceFileName.c_str() );
	if( !file.is_open() ){
		throw std::runtime_error( "could not open file " + sourceFileName );
	}

	std::string programSource( std::istreambuf_iterator<char>( file ), ( std::istreambuf_iterator<char>() ));
	cl::Program::Sources source( 1, std::make_pair( programSource.c_str(), programSource.length()+1 ));

	program = cl::Program( context, source );   
//#ifdef NDEBUG
/*work*///	err = program.build( devices, "-g -s \"D:\\MyProject\\OpenWorm\\SPH_OPENCL\\SPH_test4\\sphFluidDemo.cl\"" );
/*home*///err = program.build( devices,"-g -s \"C:\\Users\\Sergey\\Desktop\\SPH_test4\\sphFluidDemo.cl\"" );
/*#else*/
	//err = program.build( devices, "-g" );
	err = program.build( devices, "" );
/*#endif*/
	if( err != CL_SUCCESS ){
		std::string compilationErrors;
		compilationErrors = program.getBuildInfo< CL_PROGRAM_BUILD_LOG >( devices[ 0 ] );
		std::cerr << "Compilation failed: " << std::endl << compilationErrors << std::endl;
		throw std::runtime_error( "failed to build program" );
	}

	return;
}


int sph_fluid_main_start ( /*int argc, char **argv*/ )
{
	int err;
	positionBuffer = new float[ 4 * PARTICLE_COUNT ];
	velocityBuffer = new float[ 4 * PARTICLE_COUNT ];

	//preInitDX10( PARTICLE_COUNT, positionBuffer );

	try{

		initializeOpenCL( context, devices, queue, program );
//AP2012
/*
#ifdef NDEBUG
		radixSort.initializeSort( context, queue, PARTICLE_COUNT, 16, true );		
#endif
*/
//AP2012

		// initialize buffers
		gridCellsX = (int)( ( XMAX - XMIN ) / h ) + 1;
		gridCellsY = (int)( ( YMAX - YMIN ) / h ) + 1;
		gridCellsZ = (int)( ( ZMAX - ZMIN ) / h ) + 1;
		gridCellCount = gridCellsX * gridCellsY * gridCellsZ;

		acceleration = cl::Buffer( context, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer acceleration creation failed" );
		}
		gridCellIndex = cl::Buffer( context, CL_MEM_READ_WRITE, ( ( gridCellCount + 1 ) * sizeof( unsigned int ) * 1 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer gridCellIndex creation failed" );
		}
		gridCellIndexFixedUp = cl::Buffer( context, CL_MEM_READ_WRITE, ( ( gridCellCount + 1 ) * sizeof( unsigned int ) * 1 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer gridCellIndexFixedUp creation failed" );
		}
		neighborMap = cl::Buffer( context, CL_MEM_READ_WRITE, ( NK * sizeof( float ) * 2 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer neighborMap creation failed" );
		}
		particleIndex = cl::Buffer( context, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( unsigned int ) * 2 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer particleIndex creation failed" );
		}
		position = cl::Buffer( context, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer position creation failed" );
		}
		pressure = cl::Buffer( context, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 1 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer pressure creation failed" );
		}
		rho = cl::Buffer( context, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 1 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer rho creation failed" );
		}
		rhoInv = cl::Buffer( context, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 1 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer rhoInv creation failed" );
		}
		sortedPosition = cl::Buffer( context, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer sortedPosition creation failed" );
		}
		sortedVelocity = cl::Buffer( context, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer sortedVelocity creation failed" );
		}
		velocity = cl::Buffer( context, CL_MEM_READ_WRITE, ( PARTICLE_COUNT * sizeof( float ) * 4 ), NULL, &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "buffer velocity creation failed" );
		}

		// create kernels
		clearBuffers = cl::Kernel( program, "clearBuffers", &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "kernel clearBuffers creation failed" );
		}
		computeAcceleration = cl::Kernel( program, "computeAcceleration", &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "kernel computeAcceleration creation failed" );
		}
		computeDensityPressure = cl::Kernel( program, "computeDensityPressure", &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "kernel computeDensityPressure creation failed" );
		}
		findNeighbors = cl::Kernel( program, "findNeighbors", &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "kernel findNeighbors creation failed" );
		}
		hashParticles = cl::Kernel( program, "hashParticles", &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "kernel hashParticles creation failed" );
		}
		indexPostPass = cl::Kernel( program, "indexPostPass", &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "kernel indexPostPass creation failed" );
		}
		indexx = cl::Kernel( program, "indexx", &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "kernel indexx creation failed" );
		}
		integrate = cl::Kernel( program, "integrate", &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "kernel integrate creation failed" );
		}
		sortPostPass = cl::Kernel( program, "sortPostPass", &err );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "kernel sortPostPass creation failed" );
		}


		for( int i = 0; i < PARTICLE_COUNT; ++i ){
			float x, y, z;
			float r;
			r = ( (float)rand() / (float)RAND_MAX );
#define SCALE( MIN, MAX, X ) ( (MIN) + (X) * ( (MAX) - (MIN) ) )
			x = SCALE( XMIN, ( XMAX / 10 ), r );
			r = ( (float)rand() / (float)RAND_MAX );
			y = SCALE( YMIN, YMAX, r );
			r = ( (float)rand() / (float)RAND_MAX );
			z = SCALE( ZMIN, ZMAX, r );
			float * positionVector = positionBuffer + 4 * i;
			positionVector[ 0 ] = x;
			positionVector[ 1 ] = y;
			positionVector[ 2 ] = z;
			positionVector[ 3 ] = 0;
			float * velocityVector = velocityBuffer + 4 * i;
			r = ( (float)rand() / (float)RAND_MAX );
			velocityVector[ 0 ] = SCALE( -1.0f, 1.0f, r );
			r = ( (float)rand() / (float)RAND_MAX );
			velocityVector[ 1 ] = SCALE( -1.0f, 1.0f, r );
			r = ( (float)rand() / (float)RAND_MAX );
			velocityVector[ 2 ] = SCALE( -1.0f, 1.0f, r );
			velocityVector[ 3 ] = 0;
		}//for

		err = queue.enqueueWriteBuffer( position, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, positionBuffer );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "could not enqueue position write" );
		}
		err = queue.enqueueWriteBuffer( velocity, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, velocityBuffer );
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "could not enqueue velocity write" );
		}
		err = queue.finish();
		if( err != CL_SUCCESS ){
			throw std::runtime_error( "failed queue.finish" );
		}

	}catch( std::exception &e ){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}

	printf("Entering main loop\n");

	int nIter = 0;
/*
	while(1)
	{
		nIter++;
		printf("\n[[ Step %d ]]\n",nIter);
		step();
	}
*/
	//goDX10();

	//delete [] positionBuffer;
	//delete [] velocityBuffer;

	return err;//AP2012
}

int nIter=0;

extern int frames_counter;

void sph_fluid_main_step ()
{
	int c = clock();
	int work_time;
	nIter++;
	printf("\n[[ Step %d ]]",nIter);
	//printf("\n[[ Step %d ]], OpenGL_frames: %d",nIter,frames_counter);
	step();
	//printf("\nsph_fluid_main_step:%d\n",clock() - c);
}

void sph_fluid_main_stop ()
{
	delete [] positionBuffer;
	delete [] velocityBuffer;
}
