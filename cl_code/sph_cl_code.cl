
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

#define SQRT( x ) native_sqrt( x )
#define DIVIDE( a, b ) native_divide( a, b )
#define DOT( a, b ) dot( a, b )

#define NO_PARTICLE_ID -1

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
__kernel void k_init_ext_particles(
		__global struct extend_particle * ext_particles,
		int PARTICLE_COUNT
){
	int id = get_global_id(0);
	if(id > PARTICLE_COUNT){
		return;
	}
	ext_particles[id].p_id = id;
	for(int i=0;i<NEIGHBOUR_COUNT;++i){
		ext_particles[id].neighbour_list[i][0] = NO_PARTICLE_ID;
		ext_particles[id].neighbour_list[i][1] = NO_PARTICLE_ID;
	}
}

int cell_id(
		   int4 cell_factors_,
		   uint grid_cells_X,
		   uint grid_cells_Y,
		   uint grid_cells_Z//doesn't use
		   )
{
	int cell_id = cell_factors_.y + cell_factors_.z * grid_cells_Y
		+ cell_factors_.x * grid_cells_Z * grid_cells_Y;
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
							float hashGridCellSizeInv,
							int PARTICLE_COUNT,
							int OFFSET,
							int LIMIT
							)
{
	int id = get_global_id( 0 );
	if(id > PARTICLE_COUNT && id + OFFSET > LIMIT){
		return;
	}
	id += OFFSET;
	int4 cellFactors_ = cell_factors( &particles[ id ], hashGridCellSizeInv );
	int cellId_ = cell_id( cellFactors_, gridCellsX, gridCellsY, gridCellsZ ) & 0xffffff; // truncate to low 16 bits
	particles[ id ].cell_id = cellId_;
}

int getMaxIndex(
		float *d_array
)
{
	int result;
	float max_d = -1.f;
	for(int i=0; i<NEIGHBOUR_COUNT; i++){
		if (d_array[i] > max_d){
			max_d = d_array[i];
			result = i;
		}
	}
	return result;
}

// Neighbour Search algorithm function and kernel block
/** Searching for neighbour in particular spatial cell for particular particle
 *  It takes every particles from particular cell and check if it satisfy
 *  a condition that distance between particles is <= closest_distance
 */
int searchForNeighbors_b(
	int searchCell_,
	__global struct particle * particles,
	__global int * b_grid_cell_id_list,
	float4 position_,
	int myParticleId,
	int * closest_indexes,
	float * closest_distances,
	int last_farthest,
	int *found_count)
{
	int baseParticleId = b_grid_cell_id_list[ searchCell_ ];
	int i = 0;
	float _distanceSquared;
	int neighborParticleId;
	int farthest_neighbor = last_farthest;
	while( particles[i].cell_id == searchCell_){
		neighborParticleId = i;
		if(myParticleId != neighborParticleId)
		{
			float4 d = position_ - particles[neighborParticleId].pos;
			_distanceSquared = d.x*d.x + d.y*d.y + d.z*d.z; // inlined openCL dot(d,d)
			if( _distanceSquared <= closest_distances[farthest_neighbor])
			{
				closest_distances[farthest_neighbor] = _distanceSquared;
				closest_indexes[farthest_neighbor] = neighborParticleId;
				if(*found_count < NEIGHBOUR_COUNT - 1){
					(*found_count)++;
					farthest_neighbor = *found_count;
				} else{
					farthest_neighbor = getMaxIndex(closest_distances);
				}
			}
		}
		++i;
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
	int dy = deltaY;
	int dz = deltaZ * gridCellsY;
	int dx = deltaX * gridCellsZ * gridCellsY;
	int newCellId = cellId + dx + dy + dz;
	newCellId = newCellId < 0 ? newCellId + gridCellCount : newCellId;
	newCellId = newCellId >= gridCellCount ? newCellId - gridCellCount : newCellId;
	return newCellId;
}

/** Caculation spatial hash cellId for every particle
 *  Kernel fill up particleIndex buffer.
 */
int4 cellFactors(
        float4 position,
        float xmin,
        float ymin,
        float zmin,
        float hashGridCellSizeInv
)
{
    //xmin, ymin, zmin
    int4 result;
    result.x = (int)( position.x *  hashGridCellSizeInv );
    result.y = (int)( position.y *  hashGridCellSizeInv );
    result.z = (int)( position.z *  hashGridCellSizeInv );
    return result;
}

/*Fill Cell particle
 * */

__kernel void k_clear_grid_hash(
		__global int * b_grid_cell_id_list,
		uint GRID_CELL_COUNT
		){
	int id = get_global_id(0);
	if(id >= GRID_CELL_COUNT){
		return;
	}
    b_grid_cell_id_list[id] = -1;
}

__kernel void k_fill_particle_cell_hash(
		__global int * b_grid_cell_id_list,
		__global struct particle * particles,
        unsigned int DEVICE_CELL_OFFSET,
		uint PARTICLE_COUNT
){
	int id = get_global_id(0);
	if(id >= PARTICLE_COUNT){
		return;
	}
	unsigned int particle_cell_id = particles[id].cell_id - DEVICE_CELL_OFFSET;
	if(id == 0){
        b_grid_cell_id_list[particle_cell_id] = 0;
		return;
	}
	if(particles[id].cell_id != particles[id - 1].cell_id){
        b_grid_cell_id_list[particle_cell_id] = id;
	}
}

/** Searchin for neigbours foe each particles
*/
__kernel void k_neighbour_search(
		__global struct extend_particle * ext_particles,
		__global struct
				particle * particles,
        __global int * b_grid_cell_id_list,
		uint grid_cells_X,
		uint grid_cells_Y,
		uint grid_cells_Z,
		uint grid_cell_count,
        uint grid_offset,
		float h,
		float hashGridCellSize,
		float hashGridCellSizeInv,
		float simulationScale,
		float xmin,
		float ymin,
		float zmin,
		int PARTICLE_COUNT,
		int OFFSET,
		int LIMIT
){
	int id = get_global_id( 0 );
	if(id >= PARTICLE_COUNT && id + OFFSET > LIMIT){
		return;
	}
	id += OFFSET;
	float4 position_ = particles[ id ].pos;
	int myCellId = particles[id].cell_id;//& 0xffffff;// truncate to low 16 bits
	int searchCells[8];
	float r_thr2 = h * h;
	float closest_distances[NEIGHBOUR_COUNT];
	int closest_indexes[NEIGHBOUR_COUNT];
	int found_count = 0;
	for(int k=0;k<NEIGHBOUR_COUNT;k++){
		closest_distances[k] = r_thr2;
		closest_indexes[k] = -1;
	}
	searchCells[0] = myCellId - grid_offset;

	// p is the current particle position within the bounds of the hash grid
	float4 p;
	float4 p0 = (float4)( xmin, ymin, zmin, 0.0f );
	p = position_ - p0;
//
//	// cf is the min,min,min corner of the current cell
	int4 cellFactors_ = cellFactors( position_, xmin, ymin, zmin, hashGridCellSizeInv );
	float4 cf;
	cf.x = cellFactors_.x * hashGridCellSize;
	cf.y = cellFactors_.y * hashGridCellSize;
	cf.z = cellFactors_.z * hashGridCellSize;
//
	// lo.A is true if the current position is in the low half of the cell for dimension A
	int4 lo;
	lo = (( p - cf ) < h );

	int4 delta;
	int4 one = (int4)( 1, 1, 1, 1 );
	delta = one + 2 * lo;
	searchCells[1] = searchCell( myCellId, delta.x, 0, 0, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[2] = searchCell( myCellId, 0, delta.y, 0, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[3] = searchCell( myCellId, 0, 0, delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[4] = searchCell( myCellId, delta.x, delta.y, 0, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[5] = searchCell( myCellId, delta.x, 0, delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[6] = searchCell( myCellId, 0, delta.y, delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;
	searchCells[7] = searchCell( myCellId, delta.x, delta.y, delta.z, grid_cells_X, grid_cells_Y, grid_cells_Z, grid_cell_count ) - grid_offset;

 	int last_farthest = 0;
	// Search neighbour particles in every cells from searchCells list
	last_farthest = searchForNeighbors_b( searchCells[0], particles, b_grid_cell_id_list, position_,
	                                      id, closest_indexes, closest_distances,
	                                      last_farthest, &found_count );

	last_farthest = searchForNeighbors_b( searchCells[1], particles, b_grid_cell_id_list, position_,
	                                      id, closest_indexes, closest_distances,
	                                      last_farthest, &found_count );

	last_farthest = searchForNeighbors_b( searchCells[2], particles, b_grid_cell_id_list, position_,
	                                      id, closest_indexes, closest_distances,
	                                      last_farthest, &found_count );

	last_farthest = searchForNeighbors_b( searchCells[3], particles, b_grid_cell_id_list, position_,
	                                      id, closest_indexes, closest_distances,
	                                      last_farthest, &found_count );

	last_farthest = searchForNeighbors_b( searchCells[4], particles, b_grid_cell_id_list, position_,
	                                      id, closest_indexes, closest_distances,
	                                      last_farthest, &found_count );
	last_farthest = searchForNeighbors_b( searchCells[5], particles, b_grid_cell_id_list, position_,
	                                      id, closest_indexes, closest_distances,
	                                      last_farthest, &found_count );

	last_farthest = searchForNeighbors_b( searchCells[6], particles, b_grid_cell_id_list, position_,
	                                      id, closest_indexes, closest_distances,
	                                      last_farthest, &found_count );

	last_farthest = searchForNeighbors_b( searchCells[7], particles, b_grid_cell_id_list, position_,
	                                      id, closest_indexes, closest_distances,
	                                      last_farthest, &found_count );

	// storing all found neighbors into neighborMap buffer
	for(int j=0; j<NEIGHBOUR_COUNT; j++){
		float2 neighbor_data;
		neighbor_data.x = closest_indexes[j];
		if(closest_indexes[j] >= 0){
			neighbor_data.y = SQRT( closest_distances[j] ) * simulationScale; // scaled, OK
		}else{
			neighbor_data.y = -1.f;
		}
		ext_particles[id].neighbour_list[j][0] = neighbor_data.x;
		ext_particles[id].neighbour_list[j][1] = neighbor_data.y;
	}
}

//=================================
// PCI SPH KERNELS BELOW
//=================================
/** Run pcisph_computeDensity kernel
 *  The kernel's calculating density for every particle.
 */
__kernel void k_compute_density(
	__global struct extend_particle * ext_particles,
	__global struct
			particle * particles,
	float mass_mult_Wpoly6Coefficient,
	float hScaled2,
	int PARTICLE_COUNT,
	int OFFSET,
	int LIMIT
)
{
	int id = get_global_id( 0 );
	if(id >= PARTICLE_COUNT){
		return;
	}
	int nc=0;							//neighbor counter
	float density = 0.0f;
	float r_ij2;						//squared r_ij
	float hScaled6 = hScaled2*hScaled2*hScaled2;
	int real_nc = 0;

	for(int i=0;i<NEIGHBOUR_COUNT;++i){
		if( ext_particles[id].neighbour_list[i][1] != NO_PARTICLE_ID ) {
			r_ij2 = ext_particles[id].neighbour_list[i][1];
			r_ij2 *= r_ij2;
			if (r_ij2 < hScaled2) {
				density += (hScaled2 - r_ij2) * (hScaled2 - r_ij2) * (hScaled2 - r_ij2);
				//if(r_ij2>hScaled2) printf("=Error: r_ij/h = %f\n", NEIGHBOR_MAP_DISTANCE( neighborMap[ idx + nc ] ) / hScaled);
				real_nc++;
			}
		} else {
			break;
		}
	}

	if(density<hScaled6)
		density = hScaled6;
	density *= mass_mult_Wpoly6Coefficient; // since all particles are same fluid type, factor this out to here
	particles[id + OFFSET].density = density;
}
