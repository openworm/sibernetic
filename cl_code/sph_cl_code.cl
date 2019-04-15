
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

#define LIQUID_PARTICLE   1
#define ELASTIC_PARTICLE  2
#define BOUNDARY_PARTICLE 3

typedef struct particle{
#ifdef _DOUBLE_PRECISION
	double4 pos;
	double4 pos_n_1;
	double4 vel;
	double4 acceleration;
	double4 acceleration_n_1;
#else
	float4 pos;
	float4 pos_n_1;
	float4 vel;
	float4 acceleration;
	float4 acceleration_n_1;
	float4 acceleration_n_0_5;
#endif
	size_t type_;
	size_t cell_id;
#ifdef _DOUBLE_PRECISION
	double density;
	double pressure;
	double viscosity;
	double mass;
#else
	float density;
	float pressure;
	float viscosity;
	float mass;
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
		ext_particles[id - OFFSET].neighbour_list[j][0] = neighbor_data.x;
		ext_particles[id - OFFSET].neighbour_list[j][1] = neighbor_data.y;
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


/** Run pcisph_computeForcesAndInitPressure kernel
 *  The kernel initializes pressure by 0.
 *  Calculating viscosity and surface tension forces
 *  and acceleration of particle
 *  acceleration[id] = (ViscosityForces + SurfaceTensiion +GravityForces)/mass
 *  tempAcceleration[id] = acceleration[id + PARTICLE_COUNT]=0.
 */
__kernel void k_compute_forces_init_pressure(
		__global struct extend_particle * ext_particles,
		__global struct	particle * particles,
		float surf_tens_coeff,
		float mass_mult_divgradWviscosityCoefficient,
		float hScaled,
		float gravity_x,
		float gravity_y,
		float gravity_z,
		uint PARTICLE_COUNT,
		int OFFSET,
		int LIMIT
)
{
	int id = get_global_id( 0 );
	if(id >= PARTICLE_COUNT){
		return;
	}
	if(particles[ id + OFFSET].type_ == BOUNDARY_PARTICLE){
		particles[ id + OFFSET].acceleration = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
		return;
	}
	float4 acceleration_i;
	float4 result = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
	float r_ij, r_ij2;
	float4 vr_ij;
	int jd;
	float value;
	float4 accel_viscosityForce = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
	float4 vi, vj;
	float rho_i,rho_j;
	float4 accel_surfTensForce = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
	float not_bp;
	float hScaled2 = hScaled * hScaled;
	for(int i=0; i< NEIGHBOUR_COUNT; ++i)
	{
		if( (jd = (int)(ext_particles[id].neighbour_list[i][0])) != NO_PARTICLE_ID)
		{
			r_ij = ext_particles[id].neighbour_list[i][1];
			r_ij2 = r_ij * r_ij;
			if(r_ij<hScaled)
			{
				rho_i = particles[id + OFFSET].density;
				rho_j = particles[jd].density;
				vi = particles[id + OFFSET].vel;
				vj = particles[jd].vel;
				not_bp = (float)(particles[jd].type_ != BOUNDARY_PARTICLE);
				accel_viscosityForce += 1.0e-4f * (vj*not_bp-vi)*(hScaled-r_ij)/1000;
				float surffKern = (hScaled2 - r_ij2) * (hScaled2 - r_ij2) * (hScaled2 - r_ij2);
				//high resolution?
				//accel_surfTensForce += -1.7e-09f * surfTensCoeff * surffKern * (sortedPosition[id]-sortedPosition[jd]);

				//low (1/2) resolution?
				accel_surfTensForce += - 1.7e-09f * surf_tens_coeff * surffKern * (particles[id + OFFSET].pos-particles[jd].pos);
			}
		}
	}

	accel_surfTensForce.w = 0.f;
	accel_surfTensForce /= particles[id + OFFSET].mass;
	accel_viscosityForce *= 1.5f * mass_mult_divgradWviscosityCoefficient / particles[id + OFFSET].density;
	// apply external forces
	acceleration_i = accel_viscosityForce;
	acceleration_i += (float4)( gravity_x, gravity_y, gravity_z, 0.0f );
	acceleration_i +=  accel_surfTensForce; //29aug_A.Palyanov

	particles[id + OFFSET].acceleration = acceleration_i;
	// 1st half of 'acceleration' array is used to store acceleration corresponding to gravity, visc. force etc.
	particles[id + OFFSET].acceleration_n_1 = (float4)(0.0f, 0.0f, 0.0f, 0.0f );
	particles[id + OFFSET].pressure = 0.f;
}

// Boundary handling, according to the following article:
// M. Ihmsen, N. Akinci, M. Gissler, M. Teschner, Boundary Handling and Adaptive Time-stepping for PCISPH Proc. VRIPHYS, Copenhagen, Denmark, pp. 79-88, Nov 11-12, 2010.
// short citation: Ihmsen et. al., 2010
// The article chapter 3.2 describes new boundary method that combines the idea of direct-forcing [BTT09]
// with the pressure-based frozen-particles method. The proposed boundary method enforces non-penetration
// of rigid objects even for large time steps. By incorporating density estimates at the boundary into the
// pressure force, unnatural accelerations resulting from high pressure ratios are avoided.
// Boundary handling, according to the following article:
// M. Ihmsen, N. Akinci, M. Gissler, M. Teschner, Boundary Handling and Adaptive Time-stepping for PCISPH Proc. VRIPHYS, Copenhagen, Denmark, pp. 79-88, Nov 11-12, 2010.
// short citation: Ihmsen et. al., 2010
// The article chapter 3.2 describes new boundary method that combines the idea of direct-forcing [BTT09]
// with the pressure-based frozen-particles method. The proposed boundary method enforces non-penetration
// of rigid objects even for large time steps. By incorporating density estimates at the boundary into the
// pressure force, unnatural accelerations resulting from high pressure ratios are avoided.
void computeInteractionWithBoundaryParticles(
									   int id,
									   float r0,
									   __global struct extend_particle * ext_particles,
									   __global struct	particle * particles,
									   float4 * pos_,
									   bool tangVel,
									   float4 * vel
									   )
{
	//track selected particle (indices are not shuffled anymore)
	int id_b;//index of id's particle neighbour which is a boundary particle
	int id_b_source_particle, nc = 0;
	float4 n_c_i = (float4)(0.f,0.f,0.f,0.f);
	float4 n_b;
	float w_c_ib, w_c_ib_sum = 0.f, w_c_ib_second_sum = 0.f;
	float4 delta_pos;
	float n_c_i_length,x_ib_dist;

	for(int i=0; i< NEIGHBOUR_COUNT; ++i)
	{
		if( (id_b = ext_particles[id].neighbour_list[i][0]) == NO_PARTICLE_ID )
		{
			if(particles[id_b].type_ == BOUNDARY_PARTICLE){
				x_ib_dist  = ((*pos_).x - particles[id_b].pos.x) * ((*pos_).x - particles[id_b].pos.x);
				x_ib_dist += ((*pos_).y - particles[id_b].pos.y) * ((*pos_).y - particles[id_b].pos.y);
				x_ib_dist += ((*pos_).z - particles[id_b].pos.z) * ((*pos_).z - particles[id_b].pos.z);
				x_ib_dist = SQRT(x_ib_dist);
				w_c_ib = max(0.f,(r0-x_ib_dist)/r0);			//Ihmsen et. al., 2010, page 4, formula (10)
				n_b = particles[id_b].vel;			//ATTENTION! for boundary, non-moving particles velocity has no sense, but instead we need to store normal vector. We keep it in velocity data structure for memory economy.
				n_c_i += n_b * w_c_ib;							//Ihmsen et. al., 2010, page 4, formula (9)
				w_c_ib_sum += w_c_ib;							//Ihmsen et. al., 2010, page 4, formula (11), sum #1
				w_c_ib_second_sum += w_c_ib * (r0 - x_ib_dist); //Ihmsen et. al., 2010, page 4, formula (11), sum #2
			}
		}
	}
	n_c_i_length = DOT(n_c_i,n_c_i);
	if(n_c_i_length != 0.f && w_c_ib_sum != 0.f){
		n_c_i_length = sqrt(n_c_i_length);
		delta_pos = ((n_c_i/n_c_i_length) * w_c_ib_second_sum)/w_c_ib_sum;	//
		(*pos_).x += delta_pos.x;								//
		(*pos_).y += delta_pos.y;								// Ihmsen et. al., 2010, page 4, formula (11)
		(*pos_).z += delta_pos.z;								//
		if(tangVel){// tangential component of velocity
			float eps = 0.99f; //eps should be <= 1.0			// controls the friction of the collision
			float vel_n_len = n_c_i.x * (*vel).x + n_c_i.y * (*vel).y + n_c_i.z * (*vel).z;
			if(vel_n_len < 0){
				(*vel).x -= n_c_i.x * vel_n_len;
				(*vel).y -= n_c_i.y * vel_n_len;
				(*vel).z -= n_c_i.z * vel_n_len;
				(*vel) = (*vel) * eps;							// Ihmsen et. al., 2010, page 4, formula (12)
			}
		}
	}
}

/** The kernel predicts possible position value of particles
*  what leads to incompressibility. Temp value of position
*  is calculating from temp value of velocity which's tacking from predicted value of
*  tempacceleration[id].
*/
__kernel void ker_predict_positions(
		__global struct extend_particle * ext_particles,
        __global struct	particle * particles,
        float simulationScaleInv,
        float timeStep,
        float r0,
		uint PARTICLE_COUNT,
        int OFFSET,
        int LIMIT
)
{
    int id = get_global_id( 0 );
    if( id >= PARTICLE_COUNT ) return;
    float4 position_t = particles[ id  + OFFSET].pos;
    if(particles[id  + OFFSET].type_ == BOUNDARY_PARTICLE){  //stationary (boundary) particles, right?
        return;
    }
    //  pressure force (dominant)            + all other forces
    float4 acceleration_t_dt = particles[id + OFFSET].acceleration + particles[id + OFFSET].acceleration_n_1;
    float4 velocity_t = particles[id + OFFSET].vel;
//    // Semi-implicit Euler integration
    float4 velocity_t_dt = velocity_t + timeStep * acceleration_t_dt; //newVelocity_.w = 0.f;
    float posTimeStep = timeStep * simulationScaleInv;
    float4 position_t_dt = position_t + posTimeStep * velocity_t_dt;  //newPosition_.w = 0.f;
//    //sortedVelocity[id] = newVelocity_;// sorted position, as well as velocity,
//    computeInteractionWithBoundaryParticles(id, r0, ext_particles, particles, &position_t_dt, false, &velocity_t_dt);
    particles[OFFSET + id].pos_n_1 = position_t_dt;// in current version sortedPosition array has double size,
//    // PARTICLE_COUNT*2, to store both x(t) and x*(t+1)
}

/** The kernel predicts possible value of density
 *  taking into account predicted value of particle's position
 */
__kernel void ker_predict_density(
		__global struct extend_particle * ext_particles,
		__global struct	particle * particles,
		float mass_mult_Wpoly6Coefficient,
		float h,
		float simulation_scale,
		uint PARTICLE_COUNT,
		int OFFSET,
		int LIMIT
		)
{
	int id = get_global_id( 0 );
	if( id >= PARTICLE_COUNT ) return;
	float density = 0.0f;
	float density_accum = 0.0f;
	float4 r_ij;
	float r_ij2;//squared r_ij
	float h2 = h*h;
	float hScaled = h * simulation_scale;//scaled smoothing radius
	float hScaled2 = hScaled*hScaled;	//squared scaled smoothing radius
	float hScaled6 = hScaled2*hScaled2*hScaled2;
	float simulation_scale6 = simulation_scale*simulation_scale;
		  simulation_scale6 = simulation_scale6*simulation_scale6*simulation_scale6;
	int jd;
	for(int i=0;i < NEIGHBOUR_COUNT; ++i)// gather density contribution from all neighbors (if they exist)
	{
		if( (jd = (int)ext_particles[id].neighbour_list[i][0]) != NO_PARTICLE_ID )
		{
			r_ij = particles[id + OFFSET].pos - particles[jd + OFFSET].pos;
			r_ij2 = (r_ij.x*r_ij.x+r_ij.y*r_ij.y+r_ij.z*r_ij.z);
			if(r_ij2<h2)
			{
				density_accum += (h2-r_ij2)*(h2-r_ij2)*(h2-r_ij2);
			}
		}
	}
	density = density_accum * simulation_scale6;
	if(density<hScaled6)
	{
		//density += hScaled6;
		density = hScaled6;
	}
	density *= mass_mult_Wpoly6Coefficient; // since all particles are same fluid type, factor this out to here

	particles[OFFSET + id].density = density;
}
/** The kernel corrects the pressure
 *  taking into account predicted values of density.
 */
__kernel void ker_correct_pressure(
		__global struct	particle * particles,
		 float rho0,
		 float delta,
		 uint PARTICLE_COUNT,
		 int OFFSET,
		 int LIMIT
		 )
{
	int id = get_global_id( 0 );
	if( id >= PARTICLE_COUNT ) return;
	int nc = 0;// neigbor counter
	float rho_err;
	float p_corr;

	rho_err = particles[OFFSET + id].density - rho0;
	p_corr = rho_err*delta;
	if(p_corr < 0.0f) p_corr = 0.0f;//non-negative pressure
	particles[OFFSET + id].pressure += p_corr;
}

/** The kernel calculating pressure forces
 *  and calculating new value of tempacceleration[id].
 */
__kernel void ker_compute_pressure_force_acceleration(
		__global struct extend_particle * ext_particles,
		__global struct	particle * particles,
		float mass_mult_gradWspikyCoefficient,
		float h_scaled,
		float simulation_scale,
		float delta,
		float rho0,
		uint PARTICLE_COUNT,
		int OFFSET,
		int LIMIT
		)
{
	int id = get_global_id( 0 );
	if( id >= PARTICLE_COUNT ) return;
	int id_source_particle = id + OFFSET;
	if(particles[ id_source_particle ].type_ == BOUNDARY_PARTICLE){
		return;
	}
	float pressure_i  = particles[id_source_particle].pressure;
	float rho_i		  = particles[ id_source_particle ].density;
	float4 result = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
	float4 gradW_ij;
	float r_ij,rho_err;
	float4 vr_ij;
	int jd;
	float value;
	for(int i=0;i<NEIGHBOUR_COUNT; ++i)
	{
		if( (jd = ext_particles[id].neighbour_list[i][0]) != NO_PARTICLE_ID)
		{
			r_ij = ext_particles[id].neighbour_list[i][1];
			if(r_ij<h_scaled)
			{	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// Variant 1 corresponds to http://www.ifi.uzh.ch/vmml/publications/older-puclications/Solenthaler_sca08.pdf, formula (5) at page 3
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// Variant 2 corresponds http://www.ifi.uzh.ch/vmml/publications/older-puclications/Solenthaler_sca08.pdf, formula (6) at page 3
				// in more details here: http://www.ifi.uzh.ch/pax/uploads/pdf/publication/1299/Solenthaler.pdf, formula (3.3), end of page 29
				// (B. Solenthaler's dissertation "Incompressible Fluid Simulation and Advanced Surface Handling with SPH")
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/*1*/value = -(h_scaled-r_ij)*(h_scaled-r_ij)*0.5f*(particles[id + OFFSET].pressure+particles[jd].pressure)/particles[jd].density;
				/*2*///value = -(h_scaled-r_ij)*(h_scaled-r_ij)*( pressure[id]/(rho[PARTICLE_COUNT+id]*rho[PARTICLE_COUNT+id])
				/*2*///										+pressure[jd]/(rho[PARTICLE_COUNT+id]*rho[PARTICLE_COUNT+id]) );
				vr_ij = (particles[id + OFFSET].pos - particles[jd].pos)*simulation_scale; vr_ij.w = 0.0f;


				if(r_ij<0.5f*(h_scaled/2))//h_scaled/2 = r0
				{
					value = -(h_scaled*0.25f-r_ij)*(h_scaled*0.25f-r_ij)*0.5f*(rho0*delta)/particles[jd].density;
					vr_ij = (particles[id_source_particle].pos - particles[jd].pos)*simulation_scale; vr_ij.w = 0.0f;
				}

				result += value*vr_ij/r_ij;
			}
		}

	}

	result *= mass_mult_gradWspikyCoefficient / particles[OFFSET + id].density;

	particles[id + OFFSET].acceleration_n_1 = result;
}



/** The kernel run numerical integration method.
 *  Calculating value of position and velocity on step (t+1)
 *  NOTE: for now simulation using Semi-implicit Euler method
 *  	  for integration 1th order
 *  NOTE: soon we plan to add Leap-frog 2th order
 */
__kernel void ker_integrate(
		__global struct extend_particle * ext_particles,
		__global struct	particle * particles,
		float simulation_scale_inv,
		float time_step,
		float r0,
		int iteration_count,
		int mode,
		uint PARTICLE_COUNT,
		int OFFSET,
		int LIMIT
		)
{
	int id = get_global_id( 0 );
	if(id>=PARTICLE_COUNT) return;
	int id_source_particle = id + OFFSET;
	if(particles[id_source_particle].type_ == BOUNDARY_PARTICLE)
	{
		return;
	}
	if(iteration_count==0)
	{
		particles[id_source_particle].acceleration_n_0_5 = particles[id_source_particle].acceleration_n_1
				+ particles[id_source_particle].acceleration;
		return;
	}
	float4 acceleration_t = particles[id_source_particle].acceleration_n_0_5;    acceleration_t.w    = 0.f;
	float4 velocity_t = particles[id_source_particle].vel;
	int particleType = particles[ id_source_particle ].type_;
    if(mode == 2){
		float4 acceleration_t_dt = particles[id_source_particle].acceleration_n_1
		                           + particles[id_source_particle].acceleration; acceleration_t_dt.w = 0.f;
		float4 position_t = particles[id_source_particle].pos;
		float4 velocity_t_dt = velocity_t + (acceleration_t_dt)*time_step;						//
		float4 position_t_dt = position_t + (velocity_t_dt)*time_step*simulation_scale_inv;		//

//	    computeInteractionWithBoundaryParticles(id,r0, ext_particles, particles, &position_t_dt,false, &velocity_t_dt,PARTICLE_COUNT);

		particles[ id_source_particle ].vel = velocity_t_dt;
	    particles[ id_source_particle ].pos = position_t_dt;
		//velocity[ id_source_particle ] = (float4)((float)velocity_t_dt_x, (float)velocity_t_dt_y, (float)velocity_t_dt_z, 0.f);
		//position[ id_source_particle ] = (float4)((float)position_t_dt_x, (float)position_t_dt_y, (float)position_t_dt_z, particleType);

	    particles[id_source_particle].acceleration_n_0_5 = acceleration_t_dt;
		return;
	}
	/**///	float4 velocity_t_dt = velocity_t + (acceleration_t_dt)*time_step;						//
	/**///	float4 position_t_dt = position_t + (velocity_t_dt)*time_step*simulation_scale_inv;		//
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//printf("\n===[ time_step= %5e ]===",time_step);
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//	LEAPFROG METHOD		2-nd order(!)		symplectic(!)		obviously best choice			//
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// leapfrog is also symplectic
	// in the leapfrog method, the recipe changes a little bit.
    //Find the forces.
    //Find the new momentum based on the force and HALF of the small time step interval (not the whole time step)
    //Find the new position.
    //Find the next new momentum with the other half of the time step.
	//A second way to write the leapfrog looks quite different at first sight. Defining all quantities only at integer times, we can write:
	/**///	float4 position_t_dt = position_t + (velocity_t*time_step + acceleration_t*time_step*time_step/2.f)*simulation_scale_inv;
	/**///	float4 velocity_t_dt = velocity_t + (acceleration_t + acceleration_t_dt)*time_step/2.f;
	// for floats it works with a significant error, which neglects all advantages of this really nice method
	// for example, at first time step we get 2.17819E-03 instead of 2.18000E-03, and such things occur at every step and accumulate.
	// switching to doubles.
  if(mode==0/*positions_mode*/)
	{
		float4 position_t = particles[ id ].pos;
		float4 position_t_dt = position_t + (velocity_t*time_step + acceleration_t*time_step*time_step/2.f)*simulation_scale_inv;
		particles[ id_source_particle ].pos_n_1 = position_t_dt;
	}
	else
	if(mode==1/*velocities_mode*/)
	{
		float4 position_t_dt = particles[id_source_particle].pos;//necessary for computeInteractionsWithBoundaryParticles()
		float4 acceleration_t_dt = particles[id_source_particle].acceleration_n_1
		                           + particles[id_source_particle].acceleration; acceleration_t_dt.w = 0.f;
		float4 velocity_t_dt = velocity_t + (acceleration_t + acceleration_t_dt)*time_step/2.f;

		//computeInteractionWithBoundaryParticles(id,r0,neighborMap,particleIndexBack,particleIndex,position,velocity,&position_t_dt, true, &velocity_t_dt,PARTICLE_COUNT);
		particles[ id_source_particle ].vel = velocity_t_dt;
		particles[id_source_particle].acceleration_n_0_5 = acceleration_t_dt;

		particles[ id_source_particle ].pos = position_t_dt;
	}
  return;
}

