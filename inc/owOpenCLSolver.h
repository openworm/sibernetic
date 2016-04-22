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

//OpenCL solver class
class owOpenCLSolver
{
public:
	owOpenCLSolver(const float * position_cpp, const float * velocity_cpp, owConfigProperty * config, const float * elasticConnectionsData_cpp = NULL, const int * membraneData_cpp = NULL, const int * particleMembranesList_cpp = NULL);
	owOpenCLSolver(void);
	~owOpenCLSolver(void);
	//Kernels functions definition for neighbor search algorithm
	unsigned int _runClearBuffers(owConfigProperty * config);
	unsigned int _runHashParticles(owConfigProperty * config);
	void _runSort(owConfigProperty * config);
	unsigned int _runSortPostPass(owConfigProperty * config);
	unsigned int _runIndexx(owConfigProperty * config);
	void _runIndexPostPass(owConfigProperty * config);
	unsigned int _runFindNeighbors(owConfigProperty * config);
	//PCISPH kernels for physics-related calculations
	unsigned int _run_pcisph_computeDensity(owConfigProperty * config);
	unsigned int _run_pcisph_computeForcesAndInitPressure(owConfigProperty * config);
	unsigned int _run_pcisph_computeElasticForces(owConfigProperty * config);
	unsigned int _run_pcisph_predictPositions(owConfigProperty * config);
	unsigned int _run_pcisph_predictDensity(owConfigProperty * config);
	unsigned int _run_pcisph_correctPressure(owConfigProperty * config);
	unsigned int _run_pcisph_computePressureForceAcceleration(owConfigProperty * config);
	unsigned int _run_pcisph_integrate(int iterationCount, int pcisph_integrate_mode, owConfigProperty * config);
	//Kernels for membrane handling interaction
	unsigned int _run_clearMembraneBuffers(owConfigProperty * config);
	unsigned int _run_computeInteractionWithMembranes(owConfigProperty * config);
	unsigned int _run_computeInteractionWithMembranes_finalize(owConfigProperty * config);
	//
	void updateMuscleActivityData(float *_muscle_activation_signal_cpp, owConfigProperty * config);

	void read_position_buffer( float * position_cpp, owConfigProperty * config) { copy_buffer_from_device( position_cpp, position, config->getParticleCount() * sizeof( float ) * 4 ); };
	void read_velocity_buffer( float * velocity_cpp, owConfigProperty * config) { copy_buffer_from_device( velocity_cpp, velocity, config->getParticleCount() * sizeof( float ) * 4 ); };
	void read_density_buffer( float * density_cpp, owConfigProperty * config ) { copy_buffer_from_device( density_cpp, rho, config->getParticleCount() * sizeof( float ) * 1 ); }; // This need only for visualization current density of particle (graphic effect)
	void read_particleIndex_buffer( unsigned int * particleIndexBuffer, owConfigProperty * config ) { copy_buffer_from_device( particleIndexBuffer, particleIndex, config->getParticleCount() * sizeof( unsigned int ) * 2 ); }; // This need only for visualization current density of particle (graphic effect)
	void reset(const float * position_cpp, const float * velocity_cpp, owConfigProperty * config, const float * elasticConnectionsData_cpp = NULL, const int * membraneData_cpp = NULL, const int * particleMembranesList_cpp = NULL);
private:
	void create_ocl_kernel( const char *name, cl::Kernel &k );
	void create_ocl_buffer(const char *name, cl::Buffer &b, const cl_mem_flags flags,const int size);
	void copy_buffer_to_device(const void *host_b, cl::Buffer &ocl_b,const int size);
	void copy_buffer_from_device(void *host_b, const cl::Buffer &ocl_b, const int size );
	void destroy(){
		delete [] gridNextNonEmptyCellBuffer;
		delete [] _particleIndex;
	}
	// Initialization of openCl data buffers
	void initializeBuffers(const float * , const float * , owConfigProperty * , const float * , const int * , const int * );
	// Initialize OPENCL device, context, queue, program...
	void initializeOpenCL(owConfigProperty * config);
	cl::Context context;
	std::vector< cl::Device > devices;
	cl::CommandQueue		  queue;
	cl::Program				  program;
	// Buffers
	cl::Buffer  muscle_activation_signal;   // array storing data (activation signals) for an array of muscles.
                                            // now each can be activated by user independently

	cl::Buffer acceleration;                // Acceleration buffer
	cl::Buffer gridCellIndex;               // buffer with position of in particleIndex from which  located in the cell right now gridCellIndex[i] = someNumber, if cell has no particles it's equal -1
	cl::Buffer gridCellIndexFixedUp;        // the same that gridCellIndex but without empty cells
	cl::Buffer neighborMap;                 // Contains information about neighbors for all particles size = PARTICLE_COUNT * MAX_NEIGHBOR_COUNT
	cl::Buffer particleIndex;               // list of pairs [CellIndex, particleIndex]
	cl::Buffer particleIndexBack;           // list of indexes of particles before sort
	cl::Buffer position;                    // Buffer with position
	cl::Buffer pressure;                    // Pressure buffer size * (1+1extra[for membrane handling])
	cl::Buffer rho;                         // density buffer size * 2
	cl::Buffer sortedPosition;              // buffer with sorted position size * 2
	cl::Buffer sortedVelocity;              // buffer with sorted velocity size * 2
	cl::Buffer velocity;                    // buffer with velocity size * (1+1extra[for membrane handling])
	cl::Buffer elasticConnectionsData;      // list of particle pairs connected with springs and rest distance between them

	cl::Buffer membraneData;                // elementary membrane is built on 3 adjacent particles (i,j,k) and should have a form of triangle
                                            // highly recommended that i-j, j-k and k-i are already connected with springs to keep them close
                                            // to each other during whole lifetime of the simulation (user should control this by him(her)self)

	cl::Buffer particleMembranesList;       // potentially any particle can be connected with others via membrane(s)
                                            // this buffer contains MAX_MEMBRANES_INCLUDING_SAME_PARTICLE integer data cells per particle
                                            // each cell can contain -1 in case when no or no more membranes are associated with this particle,
                                            // or the index of corresponding membrane in membraneData list otherwise

	// Kernels
	cl::Kernel clearBuffers;
	cl::Kernel findNeighbors;
	cl::Kernel hashParticles;
	cl::Kernel indexx;
	cl::Kernel sortPostPass;
	//PCISPH kernels
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
	int * gridNextNonEmptyCellBuffer;
};

#endif //OW_OPENCL_SOLVER_H
