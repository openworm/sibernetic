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

#ifndef OW_OPENCL_SOLVER_H
#define OW_OPENCL_SOLVER_H

#if defined(_WIN32) || defined (_WIN64)
#pragma comment( lib, "opencl.lib" )							// opencl.lib
#endif

#if defined(__APPLE__) || defined(__MACOSX)
	#include "../inc/OpenCL/cl.hpp"
//	#include <OpenCL/cl_d3d10.h>
#else
	#include <CL/cl.hpp>
#endif

#include "owOpenCLConstant.h"
#include "owPhysicsConstant.h"
#include "owConfigProperty.h"

extern int MUSCLE_COUNT;

#if generateWormBodyConfiguration
	#define PY_NETWORK_SIMULATION
#endif

#if INTEL_OPENCL_DEBUG
//FOR DEBUGGING OPENCL CODE ON YOUR INTEL DEVICE PUT FULL PATH TO OPENCL FILE
//#define  OPENCL_DEBUG_PROGRAM_PATH "-g -s \"C:\\Users\\Serg\\Documents\\GitHub\\Smoothed-Particle-Hydrodynamics\\src\\sphFluid.cl\"" // if you debuging with intel opencl debuger you need past here full path to you opencl program
//#define  OPENCL_DEBUG_PROGRAM_PATH "-g -s \"C:\\Users\\������\\Documents\\GitHub\\Smoothed-Particle-Hydrodynamics\\src\\sphFluid.cl\"" // if you debuging with intel opencl debuger you need past here full path to you opencl program
#define  OPENCL_DEBUG_PROGRAM_PATH "-g -s \"C:\\GitHub\\Smoothed-Particle-Hydrodynamics\\src\\sphFluid.cl\"" // if you debuging with intel opencl debuger you need past here full path to you opencl program
#endif
#define OPENCL_PROGRAM_PATH "src/sphFluid.cl"

class owOpenCLSolver
{
public:
	owOpenCLSolver(const float * position_cpp, const float * velocity_cpp, owConfigProrerty * config, const float * elasticConnectionsData_cpp = NULL, const int * membraneData_cpp = NULL, const int * particleMembranesList_cpp = NULL);
	owOpenCLSolver(void);
	~owOpenCLSolver(void);
	// Initialize OPENCL device, context, queue, program...
	void initializeOpenCL(owConfigProrerty * config);
	//PCISPH kernels for data structures support and management
	//Kernels functions definition
	unsigned int _runClearBuffers(owConfigProrerty * config);
	unsigned int _runHashParticles(owConfigProrerty * config);
	unsigned int _runSort(owConfigProrerty * config);
	unsigned int _runSortPostPass(owConfigProrerty * config);
	unsigned int _runIndexx(owConfigProrerty * config);
	unsigned int _runIndexPostPass(owConfigProrerty * config);
	unsigned int _runFindNeighbors(owConfigProrerty * config);
	//PCISPH kernels for physics-related calculations
	unsigned int _run_pcisph_computeDensity(owConfigProrerty * config);
	unsigned int _run_pcisph_computeForcesAndInitPressure(owConfigProrerty * config);
	unsigned int _run_pcisph_computeElasticForces(owConfigProrerty * config);
	unsigned int _run_pcisph_predictPositions(owConfigProrerty * config);
	unsigned int _run_pcisph_predictDensity(owConfigProrerty * config);
	unsigned int _run_pcisph_correctPressure(owConfigProrerty * config);
	unsigned int _run_pcisph_computePressureForceAcceleration(owConfigProrerty * config);
	unsigned int _run_pcisph_integrate(int iterationCount, owConfigProrerty * config);
	//
	unsigned int _run_clearMembraneBuffers(owConfigProrerty * config);
	unsigned int _run_computeInteractionWithMembranes(owConfigProrerty * config);
	unsigned int _run_computeInteractionWithMembranes_finalize(owConfigProrerty * config);
	//
	unsigned int updateMuscleActivityData(float *_muscle_activation_signal_cpp);
	
	void read_position_buffer( float * position_cpp, owConfigProrerty * config) { copy_buffer_from_device( position_cpp, position, config->getParticleCount() * sizeof( float ) * 4 ); };
	void read_velocity_buffer( float * velocity_cpp, owConfigProrerty * config) { copy_buffer_from_device( velocity_cpp, velocity, config->getParticleCount() * sizeof( float ) * 4 ); };
	void read_density_buffer( float * density_cpp, owConfigProrerty * config ) { copy_buffer_from_device( density_cpp, rho, config->getParticleCount() * sizeof( float ) * 1 ); }; // This need only for visualization current density of particle (graphic effect)
	void read_particleIndex_buffer( unsigned int * particleIndexBuffer, owConfigProrerty * config ) { copy_buffer_from_device( particleIndexBuffer, particleIndex, config->getParticleCount() * sizeof( unsigned int ) * 2 ); }; // This need only for visualization current density of particle (graphic effect)
	void refresh(const float * position_cpp, const float * velocity_cpp, owConfigProrerty * config, const float * elasticConnectionsData_cpp = NULL, const int * membraneData_cpp = NULL, const int * particleMembranesList_cpp = NULL);
private:
	void create_ocl_kernel( const char *name, cl::Kernel &k );
	void create_ocl_buffer(const char *name, cl::Buffer &b, const cl_mem_flags flags,const int size);
	int copy_buffer_to_device(const void *host_b, cl::Buffer &ocl_b,const int size);
	int copy_buffer_from_device(void *host_b, const cl::Buffer &ocl_b, const int size );
	cl::Context context;
	std::vector< cl::Device > devices;
	cl::CommandQueue		  queue;
	cl::Program				  program;
	// Buffers
	cl::Buffer  muscle_activation_signal;   // array storing data (activation signals) for an array of muscles. 
											// now each can be activated by user independently

	cl::Buffer acceleration;				// forceAcceleration and pressureForceAcceleration
	cl::Buffer gridCellIndex;
	cl::Buffer gridCellIndexFixedUp;
	cl::Buffer neighborMap;
	cl::Buffer particleIndex;				// list of pairs [CellIndex, particleIndex]
	cl::Buffer particleIndexBack;			// list of indexes of particles before sort 
	cl::Buffer position;
	cl::Buffer pressure;					// size * (1+1extra[for membrane handling])
	cl::Buffer rho;							// size * 2
	cl::Buffer sortedPosition;				// size * 2
	cl::Buffer sortedVelocity;				
	cl::Buffer velocity;					// size * (1+1extra[for membrane handling])
	cl::Buffer elasticConnectionsData;		// list of particle pairs connected with springs and rest distance between them

	cl::Buffer membraneData;				// elementary membrane is built on 3 adjacent particles (i,j,k) and should have a form of triangle
											// highly recommended that i-j, j-k and k-i are already connected with springs to keep them close 
											// to each other during whole lifetime of the simulation (user should control this by him(her)self)

	cl::Buffer particleMembranesList;		// potentially any particle can be connected with others via membrane(s)
											// this buffer contains MAX_MEMBRANES_INCLUDING_SAME_PARTICLE integer data cells per particle
											// each cell can contain -1 in case when no or no more membranes are associated with this particle,
											// or the index of corresponding membrane in membraneData list othewize

	// Kernels
	cl::Kernel clearBuffers;
	cl::Kernel findNeighbors;
	cl::Kernel hashParticles;
	cl::Kernel indexx;
	cl::Kernel sortPostPass;

	cl::Kernel pcisph_computeDensity;
	cl::Kernel pcisph_computeForcesAndInitPressure;
	cl::Kernel pcisph_integrate;
	cl::Kernel pcisph_predictPositions;
	cl::Kernel pcisph_predictDensity;
	cl::Kernel pcisph_correctPressure;
	cl::Kernel pcisph_computePressureForceAcceleration;
	cl::Kernel pcisph_computeElasticForces;
	//
	//cl::Kernel prepareMembranesList;
	cl::Kernel clearMembraneBuffers;
	cl::Kernel computeInteractionWithMembranes;
	cl::Kernel computeInteractionWithMembranes_finalize;

	//Needed for sorting stuff
	int * _particleIndex;
	unsigned int * gridNextNonEmptyCellBuffer;
};

#endif //OW_OPENCL_SOLVER_H
