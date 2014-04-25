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

#pragma comment( lib, "opencl.lib" )							// opencl.lib

#if defined(__APPLE__) || defined(__MACOSX)
	#include "../inc/OpenCL/cl.hpp"
//	#include <OpenCL/cl_d3d10.h>
#else
	#include <CL/cl.hpp>
#endif
#include "owPhysicsConstant.h"

extern int PARTICLE_COUNT;
extern int PARTICLE_COUNT_RoundedUp;
extern int local_NDRange_size;
extern int MUSCLE_COUNT;

#define generateWormBodyConfiguration 0 //or load from file otherwise [0/1]
#if generateWormBodyConfiguration
	#define PY_NETWORK_SIMULATION
#endif

#if INTEL_OPENCL_DEBUG
//#define  OPENCL_DEBUG_PROGRAM_PATH "-g -s \"C:\\Users\\Serg\\Documents\\GitHub\\Smoothed-Particle-Hydrodynamics\\src\\sphFluid.cl\"" // if you debuging with intel opencl debuger you need past here full path to you opencl program
//#define  OPENCL_DEBUG_PROGRAM_PATH "-g -s \"C:\\Users\\������\\Documents\\GitHub\\Smoothed-Particle-Hydrodynamics\\src\\sphFluid.cl\"" // if you debuging with intel opencl debuger you need past here full path to you opencl program
#define  OPENCL_DEBUG_PROGRAM_PATH "-g -s \"C:\\GitHub\\Smoothed-Particle-Hydrodynamics\\src\\sphFluid.cl\"" // if you debuging with intel opencl debuger you need past here full path to you opencl program
#endif
#define OPENCL_PROGRAM_PATH "src/sphFluid.cl"

class owOpenCLSolver
{
public:
	owOpenCLSolver(const float * position_cpp, const float * velocity_cpp, const float * elasticConnectionsData_cpp = NULL, const int * membraneData_cpp = NULL, const int * particleMembranesList_cpp = NULL);
	owOpenCLSolver(void);
	~owOpenCLSolver(void);
	// Initialize OPENCL device, context, queue, program...
	void initializeOpenCL();
	//PCISPH kernels for data structures support and management
	unsigned int _runClearBuffers();
	unsigned int _runHashParticles();
	unsigned int _runSort();
	unsigned int _runSortPostPass();
	unsigned int _runIndexx();
	unsigned int _runIndexPostPass();
	unsigned int _runFindNeighbors();
	//PCISPH kernels for physics-related calculations
	unsigned int _run_pcisph_computeDensity();
	unsigned int _run_pcisph_computeForcesAndInitPressure();
	unsigned int _run_pcisph_computeElasticForces();
	unsigned int _run_pcisph_predictPositions();
	unsigned int _run_pcisph_predictDensity();
	unsigned int _run_pcisph_correctPressure();
	unsigned int _run_pcisph_computePressureForceAcceleration();
	unsigned int _run_pcisph_integrate(int iterationCount);
	//
	unsigned int _run_clearMembraneBuffers();
	unsigned int _run_computeInteractionWithMembranes();
	unsigned int _run_computeInteractionWithMembranes_finalize();
	//
	unsigned int updateMuscleActivityData(float *_muscle_activation_signal_cpp);
	
	void read_position_buffer( float * position_cpp ) { copy_buffer_from_device( position_cpp, position, PARTICLE_COUNT * sizeof( float ) * 4 ); };
	void read_density_buffer( float * density_cpp ) { copy_buffer_from_device( density_cpp, rho, PARTICLE_COUNT * sizeof( float ) * 1 ); }; // This need only for visualization current density of particle (graphic effect)
	void read_particleIndex_buffer( unsigned int * particleIndexBuffer ) { copy_buffer_from_device( particleIndexBuffer, particleIndex, PARTICLE_COUNT * sizeof( unsigned int ) * 2 ); }; // This need only for visualization current density of particle (graphic effect)
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
	cl::Kernel computeAcceleration;
	cl::Kernel computeDensityPressure;
	cl::Kernel findNeighbors;
	cl::Kernel hashParticles;
	cl::Kernel indexx;
	cl::Kernel integrate;
	cl::Kernel sortPostPass;
	// Additional kernels for PCISPH and for calculation elastic forces
	//cl::Kernel preElasticMatterPass;
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

};

#endif //OW_OPENCL_SOLVER_H
