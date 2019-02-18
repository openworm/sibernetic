
/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2017 OpenWorm.
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
 
#ifdef cl_amd_printf
	#pragma OPENCL EXTENSION cl_amd_printf : enable
	#define PRINTF_ON // this comment because using printf leads to very slow work on Radeon r9 290x on my machine
    	              // don't know why
#elif defined(cl_intel_printf)
	#pragma OPENCL EXTENSION cl_intel_printf : enable
	#define PRINTF_ON
#endif

#ifdef _DOUBLE_PRECISION
#ifdef cl_khr_fp64
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable
	
#elif defined(cl_amd_fp64)
	#pragma OPENCL EXTENSION cl_amd_fp64 : enable
	#define _DOUBLE_PRECISION
#else
	#error "Double precision floating point not supported by OpenCL implementation."
#endif
#endif
#if defined(_WIN32) || defined(_WIN64)
#include "inc\\ocl_struct.h"
#else
#include "inc/ocl_struct.h"
#endif

typedef struct particle{
#ifdef _DOUBLE_PRECISION
	double4 pos;
	double4 vel;
#else
	float4 pos;
	float4 vel;
#endif
	size_t type_;
	size_t cell_id;
#ifdef _DOUBLE_PRECISION
	double density;
	double pressure;
#else
	float density;
	float pressure;
#endif
} particle;


/** Just for test
*/
__kernel void k_check_copy(__global struct extend_particle * ext_particles,
							   __global struct particle	* particles){
	int id = get_global_id(0);
#ifdef PRINTF_ON
	if(id == 0){
		printf("sizeof() of particles_f is %d\n", sizeof(particle) );
	}
	if(id == 0 && particles[0].pos.x == 1.67 && particles[0].pos.y == 1.67 && particles[0].pos.z == 1.67 ){
		printf("\nTEST PASSED.\n");
	}
#endif
}

/**Initialization of neighbour list by -1 
* what means that it's no neighbours. 
*/
__kernel void k_init_ext_particles(__global struct extend_particle * ext_particles){
	int id = get_global_id(0);
	ext_particles[id].p_id = id;
	for(int i=0;i<NEIGHBOUR_COUNT;++i){
		ext_particles[id].neighbour_list[i] = -1;
	}
}

/** Calc current cell id for each particles
*/
__kernel void k_calc_cell_id(__global struct
								particle * particles){

}

/** Searchin for neigbours foe each particles
*/
__kernel void k_neighbour_search(__global struct extend_particle * ext_particles,
							   __global struct 
							   	particle * particles){
}

int cell_id(
		   int4 cell_factors_,
		   uint grid_cells_X,
		   uint grid_cells_Y,
		   uint grid_cells_Z//doesn't use
		   )
{
	int cell_id = cell_factors_.x + cell_factors_.y * grid_cells_X
		+ cell_factors_.z * grid_cells_X * grid_cells_Y;
	return cell_id;
}
/** Caculation spatial hash cellId for every particle
 *  Kernel fill up particleIndex buffer.
 */
 int4 cell_factors(
				 __global struct 
				 particle * particle,
				 float x_min,
				 float y_min,
				 float z_min,
				 float hash_grid_cell_size_inv
				 )
{
	//xmin, ymin, zmin
	int4 result;
	result.x = (int)( particle->pos.x *  hash_grid_cell_size_inv );
	result.y = (int)( particle->pos.y *  hash_grid_cell_size_inv );
	result.z = (int)( particle->pos.z *  hash_grid_cell_size_inv );
	return result;
}
__kernel void k_hash_particles(
							__global struct 
							particle * particles,
							uint gridCellsX,
							uint gridCellsY,
							uint gridCellsZ,
							float hashGridCellSizeInv,
							float xmin,
							float ymin,
							float zmin,
							__global uint2 * particleIndex,
							uint   PARTICLE_COUNT
							)
{
	int id = get_global_id( 0 );
	if( id >= PARTICLE_COUNT ) return;
	int4 cellFactors_ = cell_factors( &particles[ id ], xmin, ymin, zmin, hashGridCellSizeInv );
	int cellId_ = cell_id( cellFactors_, gridCellsX, gridCellsY, gridCellsZ ) & 0xffffff; // truncate to low 16 bits
	particles[ id ].cell_id = cellId_;
}

