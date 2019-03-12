
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
__kernel void k_init_ext_particles(__globalPOSITION_CELL_ID( position_ ) struct extend_particle * ext_particles){
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

int cell_id(
		   int4 cell_factors_,
		   uint grid_cells_X,
		   uint grid_cells_Y,
		   uint grid_cells_Z//doesn't use
		   )
{
	int cell_id = cell_factors_.y + cell_factors_.x * grid_cells_Y
		+ cell_factors_.z * grid_cells_X * grid_cells_Y;
	return cell_id;
}
/** Caculation spatial hash cellId for every particle
 *  Kernel fill up particleIndex buffer.
 */
 int4 cell_factors(
				 __global struct 
				 particle * particle,
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
							__global struct particle * particles,
							uint gridCellsX,
							uint gridCellsY,
							uint gridCellsZ,
							float hashGridCellSizeInv
							)
{
	int id = get_global_id( 0 );
	int4 cellFactors_ = cell_factors( &particles[ id ], hashGridCellSizeInv );
	int cellId_ = cell_id( cellFactors_, gridCellsX, gridCellsY, gridCellsZ ) & 0xffffff; // truncate to low 16 bits
	particles[ id ].cell_id = cellId_;
}

// Neighbour Search algorithm function and kernel block
/** Searching for neighbour in particular spatial cell for particular particle
 *  It takes every particles from particular cell and check if it satisfy
 *  a condition that distance between particles is <= closest_distance
 */
int searchForNeighbors_b(
	int searchCell_,
	__global uint * gridCellIndex,
	float4 position_,
	int myParticleId,
	__global float4 * sortedPosition,
	__global float2 * neighborMap,
	int * closest_indexes,
	float * closest_distances,
	int last_farthest,
	int *found_count){
	int baseParticleId = gridCellIndex[ searchCell_ ];
	int nextParticleId = gridCellIndex[ searchCell_ + 1 ];
	int particleCountThisCell = nextParticleId - baseParticleId; // Calcuating how many particle is containining in particular cell
	int i = 0;
	float _distanceSquared;
	int neighborParticleId;
	int farthest_neighbor = last_farthest;
	while( i < particleCountThisCell ){
		neighborParticleId = baseParticleId + i;
		if(myParticleId != neighborParticleId)
		{
			float4 d = position_ - sortedPosition[ neighborParticleId ];
			_distanceSquared = d.x*d.x + d.y*d.y + d.z*d.z; // inlined openCL dot(d,d)
			if( _distanceSquared <= closest_distances[farthest_neighbor])
			{
				closest_distances[farthest_neighbor] = _distanceSquared;
				closest_indexes[farthest_neighbor] = neighborParticleId;
				if(*found_count < MAX_NEIGHBOR_COUNT-1){
					(*found_count)++;
					farthest_neighbor = *found_count;
				} else{
					farthest_neighbor = getMaxIndex(closest_distances);
				}
			}
		}
		i++;
	}//while
	return farthest_neighbor;
}
/** Return value of cellId from gridCellIndexFixedUp
 *  for particular cell and offset (deltaX,Y,Z)
 */
int searchCell(
		int cellId,
		int deltaX,
		int deltaY,
		int deltaZ,
		uint gridCellsX,
		uint gridCellsY,
		uint gridCellsZ,
		uint gridCellCount
)
{
	int dx = deltaX;
	int dy = deltaY * gridCellsX;
	int dz = deltaZ * gridCellsX * gridCellsY;
	int newCellId = cellId + dx + dy + dz;
	newCellId = newCellId < 0 ? newCellId + gridCellCount : newCellId;
	newCellId = newCellId >= gridCellCount ? newCellId - gridCellCount : newCellId;
	return newCellId;
}


/*Fill Cell particle
 * */

/** Searchin for neigbours foe each particles
*/
__kernel void k_neighbour_search(
		__global struct extend_particle * ext_particles,
		__global struct
				particle * particles,
		uint grid_cells_X,
		uint grid_cells_Y,
		uint grid_cells_Z,
		float h,
		float hashGridCellSize,
		float hashGridCellSizeInv,
		float simulationScale,
		float xmin,
		float ymin,
		float zmin,
){
	int id = get_global_id( 0 );
	__global uint * gridCellIndex = gridCellIndexFixedUp;
	float4 position_ = particles[ id ];
	int myCellId = (int)position_.cell_id & 0xffffff;// truncate to low 16 bits
	int searchCells[8];
	float r_thr2 = h * h;
	float closest_distances[NEIGHBOUR_COUNT];
	int closest_indexes[NEIGHBOUR_COUNT];
	int found_count = 0;
	for(int k=0;k<NEIGHBOUR_COUNT;k++){
		closest_distances[k] = r_thr2;
		closest_indexes[k] = -1;
	}
	searchCells[0] = myCellId;

	// p is the current particle position within the bounds of the hash grid
//	float4 p;
//	float4 p0 = (float4)( xmin, ymin, zmin, 0.0f );
//	p = position_ - p0;
//
//	// cf is the min,min,min corner of the current cell
//	int4 cellFactors_ = cellFactors( position_, xmin, ymin, zmin, hashGridCellSizeInv );
//	float4 cf;
//	cf.x = cellFactors_.x * hashGridCellSize;
//	cf.y = cellFactors_.y * hashGridCellSize;
//	cf.z = cellFactors_.z * hashGridCellSize;
//
//	// lo.A is true if the current position is in the low half of the cell for dimension A
//	int4 lo;
//	lo = (( p - cf ) < h );
//
//	int4 delta;
//	int4 one = (int4)( 1, 1, 1, 1 );
//	delta = one + 2 * lo;
	//searchCells[1] = searchCell( myCellId, delta.x, 0, 0, gridCellsX, gridCellsY, gridCellsZ, gridCellCount );
	int last_farthest = 0;
	// Search neighbour particles in every cells from searchCells list
	last_farthest = searchForNeighbors_b( searchCells[0], gridCellIndex, position_,
										  id, sortedPosition, neighborMap,
										  closest_indexes, closest_distances, last_farthest, &found_count );

	// storing all found neighbors into neighborMap buffer
	for(int j=0; j<NEIGHBOUR_COUNT; j++){
		float2 neighbor_data;
		neighbor_data.x = closest_indexes[j];
		if(closest_indexes[j] >= 0){
			neighbor_data.y = SQRT( closest_distances[j] ) * simulationScale; // scaled, OK
		}else{
			neighbor_data.y = -1.f;
		}
		ext_particles[id].neighbour_list[j] = neighbor_data;
	}
}