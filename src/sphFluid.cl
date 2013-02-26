// sphFluid.cl  OpenWorm Project   31/2/2013
//
// Equations referenced here are from:
// "Particle-based fluid simulation for interactive applications", Muller, Charypar & Gross,
// Eurographics/SIGGRAPH Symposium on Computer Animation (2003).
//TODO: write all papers and dissertation what we use in work

#include "src//owOpenCLConstant.h"
//#define PARTICLE_COUNT ( 32 * 1024 )//( 32 * 1024 )
//#define NEIGHBOR_COUNT 32



#define POSITION_CELL_ID( i ) i.w

#define PI_CELL_ID( name ) name.x
#define PI_SERIAL_ID( name ) name.y

#define NEIGHBOR_MAP_ID( nm ) nm.x
#define NEIGHBOR_MAP_DISTANCE( nm ) nm.y

#define RHO( i ) i.x
#define RHO_INV( i ) i.y
//#define P( i ) i.z

#define DIVIDE( a, b ) native_divide( a, b )
#define SQRT( x ) native_sqrt( x )
#define DOT( a, b ) dot( a, b )


#if 1
#define SELECT( A, B, C ) select( A, B, (C) * 0xffffffff )
#else
#define SELECT( A, B, C ) C ? B : A
#endif

#pragma OPENCL EXTENSION cl_intel_printf : enable

__kernel void clearBuffers(
						   __global float2 * neighborMap
						   )
{
	int id = get_global_id( 0 );
	__global float4 * nm = (__global float4 *)neighborMap;
	int outIdx = ( id * NEIGHBOR_COUNT ) >> 1;//int4 versus int2 addressing
	float4 fdata = (float4)( -1, -1, -1, -1 );

	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
	nm[ outIdx++ ] = fdata;
}

// Gradient of equation 21.  Vector result.

float4 gradWspiky(
				  float r,
				  float h,
				  float gradWspikyCoefficient,
				  float4 position_i,
				  float4 position_j,
				  float simulationScale
				  )
{
	float4 rVec = position_i - position_j; rVec.w = 0.f;
	float4 scaledVec = rVec * simulationScale;
	if(r>=0.f)
	scaledVec /= r; 
	//rVec.w = 0.0f;  // what for? never used any longer
	float x = h - r; // r & h are scaled, as I see at least in basic SPH mode
	float4 result = 0.f;
	if(x>=0.f) result = (x) * (x) * scaledVec * gradWspikyCoefficient;
	return result;
}

float4 contributeGradP(
					   int id,
					   int neighborParticleId,						
					   float p_i,
					   float p_j,
					   float rho_j_inv,
					   float4 position_i,
					   __global float * pressure,
					   __global float * rho,
					   __global float4 * sortedPosition,
					   float r,
					   float mass,
					   float h,
					   float gradWspikyCoefficient,
					   float simulationScale
					   )
{
	// Following Muller Charypar and Gross ( 2003 ) Particle-Based Fluid Simulation for Interactive Applications
	// -grad p_i = - sum_j m_j ( p_i + p_j ) / ( 2 rho_j ) grad Wspiky
	// Equation 10.
	// AP2012: seem to be en error here: for correct acceleration dimension [m/s^2]  rho in denominator should be in 2nd degree
	// mass*pressure*gradW/rho^2 = (kg*Neuton/m^2)*(1/m^4)/(kg/m^3)^2 = (kg*kg*m/(m*s^2))*(1/m^4)/(kg/m^3)^2 = 1/(m*s^2)*(m^6/m^4) = m/s^2 [OK]
	// but here rho is only in 1st degree, not 2nd one
	// ===
	// Finally everything is ok here. gradP is further multiplied by rho_inv which makes acceleration dimensions correct


	float4 neighborPosition;
	neighborPosition = sortedPosition[ neighborParticleId ];
	float4 smoothingKernel = gradWspiky( r, h, gradWspikyCoefficient, position_i, neighborPosition, simulationScale );
	float4 result = mass * ( p_i + p_j ) * 0.5f * rho_j_inv * smoothingKernel;
	return result;
}


// Laplacian of equation 22.  Scalar result.

float del2Wviscosity(
					 float r,
					 float h,
					 float del2WviscosityCoefficient
					 )
{
	// equation 22
	float result = 0.f;
	if((r>=0)&&(r<h)) result = ( h - r ) * del2WviscosityCoefficient;
	return result;
}



float4 contributeDel2V(
					   int id,
					   float4 v_i,
					   int neighborParticleId,
					   __global float4 * sortedVelocity,
					   float rho_j_inv,
					   float r,
					   float mass,
					   float h,
					   float del2WviscosityCoefficient
					   )
{
	// mu del^2 v = mu sum_j m_j ( v_j - v_i ) / rho_j del^2 Wviscosity
	// Equation 14.
	float4 v_j = sortedVelocity[ neighborParticleId ];
	float4 d = v_j - v_i;
	d.w = 0.f;
	float4 result = mass * d * rho_j_inv * del2Wviscosity( r, h, del2WviscosityCoefficient );
	return result;
}



// Mueller et al equation 3.  Scalar result.

float Wpoly6(
			 float rSquared,
			 float hSquared,
			 float Wpoly6Coefficient
			 )
{
	float x = hSquared - rSquared;
	float result = 0.f;
	if(x>0) result = x * x * x * Wpoly6Coefficient;
	return result;
}

float densityContribution(
						  int idx,
						  int i,
						  __global float2 * neighborMap,
						  float mass,
						  float hSquared,
						  float Wpoly6Coefficient
						  )
{
	float2 nm = neighborMap[ idx + i ];
	int neighborParticleId = NEIGHBOR_MAP_ID( nm );
	float r = NEIGHBOR_MAP_DISTANCE( nm );	
	float smoothingKernel = Wpoly6( r*r, hSquared, Wpoly6Coefficient );
	float result = SELECT( smoothingKernel, 0.0f, ( neighborParticleId == NO_PARTICLE_ID ) );
	return result;
}


int searchCell( 
			   int cellId,
			   int deltaX,
			   int deltaY,
			   int deltaZ,
			   int gridCellsX, 
			   int gridCellsY, 
			   int gridCellsZ,
			   int gridCellCount
			   )
{
	int dx = deltaX;
	int dy = deltaY * gridCellsX;
	int dz = deltaZ * gridCellsX * gridCellsY;
	int newCellId = cellId + dx + dy + dz;
	newCellId = SELECT( newCellId, newCellId + gridCellCount, ( newCellId < 0 ) );
	newCellId = SELECT( newCellId, newCellId - gridCellCount, ( newCellId >= gridCellCount ) );
	return newCellId;
}


#define FOUND_NO_NEIGHBOR 0
#define FOUND_ONE_NEIGHBOR 1

#define radius_segments 30


int searchForNeighbors( 
					   int searchCell_, 
					   __global uint * gridCellIndex, 
					   float4 position_, 
					   int myParticleId, 
					   __global float4 * sortedPosition,
					   __global float2 * neighborMap,
					   int spaceLeft,
					   float h,
					   float simulationScale,
					   int mode,
					   int * radius_distrib,
					   float r_thr
					   )
{
	int baseParticleId = gridCellIndex[ searchCell_ ];
	int nextParticleId = gridCellIndex[ searchCell_ + 1 ];
	int particleCountThisCell = nextParticleId - baseParticleId;
	int potentialNeighbors = particleCountThisCell;
	int foundCount = 0;
	bool loop = true;
	int i = 0,j;
	float _distance,_distanceSquared;
	float r_thr_Squared = r_thr*r_thr;
	float2 neighbor_data;
	int neighborParticleId;
	int myOffset;
	if(spaceLeft>0)
		while( i < particleCountThisCell ){

			neighborParticleId = baseParticleId + i;
			
			//if(myParticleId == neighborParticleId) continue;

			if(myParticleId != neighborParticleId)
			{
				float4 d = position_ - sortedPosition[ neighborParticleId ];
				d.w = 0.0f;
				_distanceSquared = DOT( d, d );
				if( _distanceSquared <= r_thr_Squared )
				{
					_distance = SQRT( _distanceSquared );
					j = (int)(_distance*radius_segments/h);
					if(j<radius_segments) radius_distrib[j]++; 

					// searchForNeighbors runs twice
					// first time with mode = 0, to build distribution
					// and 2nd time with mode = 1, to select 32 nearest neighbors
					if(mode)
					{
						myOffset = NEIGHBOR_COUNT - spaceLeft + foundCount;
						if(myOffset>=NEIGHBOR_COUNT) break;// New line fixing the bug with indeterminism. A. Palyanov 22.02.2013
						neighbor_data.x = neighborParticleId;
						neighbor_data.y = _distance * simulationScale; // scaled, OK
						neighborMap[ myParticleId*NEIGHBOR_COUNT + myOffset ] = neighbor_data;
						foundCount++;
					}
				}
			
			}

			i++;
			
		}//while

	return foundCount;
}


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





__kernel void findNeighbors(
							__global uint * gridCellIndexFixedUp,
							__global float4 * sortedPosition,
							int gridCellCount,
							int gridCellsX,
							int gridCellsY,
							int gridCellsZ,
							float h,
							float hashGridCellSize,
							float hashGridCellSizeInv,
							float simulationScale,
							float xmin,
							float ymin,
							float zmin,
							__global float2 * neighborMap
							//int   nIter
							)
{
	__global uint * gridCellIndex = gridCellIndexFixedUp;
	int id = get_global_id( 0 );
	float4 position_ = sortedPosition[ id ];
	int myCellId = (int)POSITION_CELL_ID( position_ ) & 0xffff;// truncate to low 16 bits
	int searchCell_;
	int foundCount = 0;
	int mode = 0;
	int distrib_sum = 0;
	int radius_distrib[radius_segments];
	int i=0,j;
	float r_thr = h;
	
	while( i<radius_segments )
	{
		radius_distrib[i]=0;
		i++;
	}
	
	while( mode<2 )
	{
		// search surrounding cell 1

		searchCell_ = myCellId;
		foundCount += searchForNeighbors( searchCell_, gridCellIndex, position_, 
			id, sortedPosition, neighborMap, NEIGHBOR_COUNT - foundCount, 
			h, simulationScale, mode, radius_distrib, r_thr );

		// p is the current particle position within the bounds of the hash grid
		float4 p;
		float4 p0 = (float4)( xmin, ymin, zmin, 0.0f );
		p = position_ - p0;

		// cf is the min,min,min corner of the current cell
		int4 cellFactors_ = cellFactors( position_, xmin, ymin, zmin, hashGridCellSizeInv );
		float4 cf;
		cf.x = cellFactors_.x * hashGridCellSize;
		cf.y = cellFactors_.y * hashGridCellSize;
		cf.z = cellFactors_.z * hashGridCellSize;

		// lo.A is true if the current position is in the low half of the cell for dimension A
		int4 lo;
		lo = (( p - cf ) < h );

		int4 delta;
		int4 one = (int4)( 1, 1, 1, 1 );
		delta = one + 2 * lo;

		// search surrounding cells 2..8

		searchCell_ = searchCell( myCellId, delta.x, 0, 0, gridCellsX, gridCellsY, gridCellsZ, gridCellCount );
		foundCount += searchForNeighbors( searchCell_, gridCellIndex, position_, 
			id, sortedPosition, neighborMap, NEIGHBOR_COUNT - foundCount, 
			h, simulationScale, mode, radius_distrib, r_thr  );

		searchCell_ = searchCell( myCellId, 0, delta.y, 0, gridCellsX, gridCellsY, gridCellsZ, gridCellCount );
		foundCount += searchForNeighbors( searchCell_, gridCellIndex, position_, 
			id, sortedPosition, neighborMap, NEIGHBOR_COUNT - foundCount, 
			h, simulationScale, mode, radius_distrib, r_thr  );

		searchCell_ = searchCell( myCellId, 0, 0, delta.z, gridCellsX, gridCellsY, gridCellsZ, gridCellCount );
		foundCount += searchForNeighbors( searchCell_, gridCellIndex, position_, 
			id, sortedPosition, neighborMap, NEIGHBOR_COUNT - foundCount, 
			h, simulationScale, mode, radius_distrib, r_thr  );

		searchCell_ = searchCell( myCellId, delta.x, delta.y, 0, gridCellsX, gridCellsY, gridCellsZ, gridCellCount );
		foundCount += searchForNeighbors( searchCell_, gridCellIndex, position_, 
			id, sortedPosition, neighborMap, NEIGHBOR_COUNT - foundCount, 
			h, simulationScale, mode, radius_distrib, r_thr  );

		searchCell_ = searchCell( myCellId, delta.x, 0, delta.z, gridCellsX, gridCellsY, gridCellsZ, gridCellCount );
		foundCount += searchForNeighbors( searchCell_, gridCellIndex, position_, 
			id, sortedPosition, neighborMap, NEIGHBOR_COUNT - foundCount, 
			h, simulationScale, mode, radius_distrib, r_thr  );

		searchCell_ = searchCell( myCellId, 0, delta.y, delta.z, gridCellsX, gridCellsY, gridCellsZ, gridCellCount );
		foundCount += searchForNeighbors( searchCell_, gridCellIndex, position_, 
			id, sortedPosition, neighborMap, NEIGHBOR_COUNT - foundCount, 
			h, simulationScale, mode, radius_distrib, r_thr  );

		searchCell_ = searchCell( myCellId, delta.x, delta.y, delta.z, gridCellsX, gridCellsY, gridCellsZ, gridCellCount );
		foundCount += searchForNeighbors( searchCell_, gridCellIndex, position_, 
			id, sortedPosition, neighborMap, NEIGHBOR_COUNT - foundCount, 
			h, simulationScale, mode, radius_distrib, r_thr );

		if(mode==0)
		{
			j=0;
			
			while(j<radius_segments)
			{
				distrib_sum += radius_distrib[j];
				if(distrib_sum==NEIGHBOR_COUNT) break;
				if(distrib_sum> NEIGHBOR_COUNT) { j--; break; }
				j++;
			}

			r_thr = (j+1)*h/radius_segments;

		}

		mode++;
	}

}


int cellId( 
		   int4 cellFactors_,
		   int gridCellsX,
		   int gridCellsY,
		   int gridCellsZ//don't use
		   )
{
	int cellId_ = cellFactors_.x + cellFactors_.y * gridCellsX
		+ cellFactors_.z * gridCellsX * gridCellsY;
	return cellId_;
}



__kernel void hashParticles(
							__global float4 * position,
							int gridCellsX,
							int gridCellsY,
							int gridCellsZ,
							float hashGridCellSizeInv,
							float xmin,
							float ymin,
							float zmin,
							__global uint2 * particleIndex
							)
{
	int id = get_global_id( 0 );
	if( id >= PARTICLE_COUNT ){
		uint2 result;
		int gridCellCount = gridCellsX * gridCellsY * gridCellsZ;
		PI_CELL_ID( result ) = gridCellCount + 1;
		PI_SERIAL_ID( result ) = id;
		particleIndex[ id ] = result;
		return;
	}

	//position[id].w = 0.f;
	//position[PARTICLE_COUNT+id].w = 0.f;

	float4 _position = position[ id ];
	int4 cellFactors_ = cellFactors( _position, xmin, ymin, zmin, hashGridCellSizeInv ); 
	//int cellId_ = cellId( cellFactors_, gridCellsX, gridCellsY, gridCellsZ );
	int cellId_ = cellId( cellFactors_, gridCellsX, gridCellsY, gridCellsZ ) & 0xffff; // truncate to low 16 bits
	uint2 result;
	PI_CELL_ID( result ) = cellId_;
	PI_SERIAL_ID( result ) = id;
	particleIndex[ id ] = result;

}


__kernel void indexx(
					 __global uint2 * particleIndex,
					 int gridCellCount,
					 __global uint * gridCellIndex
					 )
{
	//fill up gridCellIndex
	int id = get_global_id( 0 );
	if( id > gridCellCount  ){
		return;
	}

	if( id == gridCellCount ){
		// add the nth+1 index value
		gridCellIndex[ id ] = PARTICLE_COUNT;
		return;
	}		
	if( id == 0 ){
		gridCellIndex[ id ] = 0;
		return;
	}

	// binary search for the starting position in sortedParticleIndex
	int low = 0;
	int high = PARTICLE_COUNT - 1;
	bool converged = false;

	int cellIndex = NO_PARTICLE_ID;
	while( !converged ){
		if( low > high ){
			converged = true;
			cellIndex = NO_PARTICLE_ID;
			continue;
		}

		int idx = ( high - low ) * 0.5f + low;
		uint2 sample = particleIndex[ idx ];
		int sampleCellId = PI_CELL_ID( sample );
		bool isHigh = ( sampleCellId > id );
		high = SELECT( high, idx - 1, isHigh );
		bool isLow = ( sampleCellId < id );
		low = SELECT( low, idx + 1, isLow );
		bool isMiddle = !( isHigh || isLow );

		uint2 zero2 = (uint2)( 0, 0 );
		uint2 sampleMinus1;
		int sampleM1CellId = 0;
		bool zeroCase = ( idx == 0 && isMiddle ); //it means that we in middle or 
		sampleMinus1 = SELECT( (uint2)particleIndex[ idx - 1 ], zero2, (uint2)zeroCase );//if we in middle this return zero2 else (uint2)particleIndex[ idx - 1 ]
		sampleM1CellId = SELECT( PI_CELL_ID( sampleMinus1 ), (uint)(-1), zeroCase );//if we in middle this return (uint)(-1) else sampleMinus1.x (index of cell)
		bool convergedCondition = isMiddle && ( zeroCase || sampleM1CellId < sampleCellId );
		converged = convergedCondition;
		cellIndex = SELECT( cellIndex, idx, convergedCondition );
		high = SELECT( high, idx - 1, ( isMiddle && !convergedCondition ) );
	}//while

	gridCellIndex[ id ] = cellIndex;//
}



void handleBoundaryConditions(
							  float4 position,
							  float4 * newVelocity,
							  float timeStep,
							  float4 * newPosition,
							  float xmin,
							  float xmax,
							  float ymin,
							  float ymax,
							  float zmin,
							  float zmax,
							  float damping
							  )
{
	if( (*newPosition).x < xmin ){
		float intersectionDistance = -position.x / (*newVelocity).x;
		float4 intersection = position + intersectionDistance * *newVelocity;
		float4 normal = (float4)( 1, 0, 0, 0 );
		float4 reflection = *newVelocity - 2.0f * DOT( *newVelocity, normal ) * normal;
		float remaining = timeStep - intersectionDistance;
		position = intersection;
		*newVelocity = reflection;
		*newPosition = intersection + remaining * damping * reflection;
	}
	else if( (*newPosition).x > xmax ){
		float intersectionDistance = ( xmax - position.x ) / (*newVelocity).x;
		float4 intersection = position + intersectionDistance * *newVelocity;
		float4 normal = (float4)( -1, 0, 0, 0 );
		float4 reflection = *newVelocity - 2.0f * DOT( *newVelocity, normal ) * normal;
		float remaining = timeStep - intersectionDistance;
		position = intersection;
		*newVelocity = reflection;
		*newPosition = intersection + remaining * damping * reflection;
	}

	if( (*newPosition).y < ymin ){
		float intersectionDistance = -position.y / (*newVelocity).y;
		float4 intersection = position + intersectionDistance * *newVelocity;
		float4 normal = (float4)( 0, 1, 0, 0 );
		float4 reflection = *newVelocity - 2.0f * DOT( *newVelocity, normal ) * normal;
		float remaining = timeStep - intersectionDistance;
		position = intersection;
		*newVelocity = reflection;
		*newPosition = intersection + remaining * damping * reflection;
	}
	else if( (*newPosition).y > ymax ){
		float intersectionDistance = ( ymax - position.y ) / (*newVelocity).y;
		float4 intersection = position + intersectionDistance * *newVelocity;
		float4 normal = (float4)( 0, -1, 0, 0 );
		float4 reflection = *newVelocity - 2.0f * DOT( *newVelocity, normal ) * normal;
		float remaining = timeStep - intersectionDistance;
		position = intersection;
		*newVelocity = reflection;
		*newPosition = intersection + remaining * damping * reflection;
	}

	if( (*newPosition).z < zmin ){
		float intersectionDistance = -position.z / (*newVelocity).z;
		float4 intersection = position + intersectionDistance * *newVelocity;
		float4 normal = (float4)( 0, 0, 1, 0 );
		float4 reflection = *newVelocity - 2.0f * DOT( *newVelocity, normal ) * normal;
		float remaining = timeStep - intersectionDistance;
		position = intersection;
		*newVelocity = reflection;
		*newPosition = intersection + remaining * damping * reflection;
	}
	else if( (*newPosition).z > zmax ){
		float intersectionDistance = ( zmax - position.z ) / (*newVelocity).z;
		float4 intersection = position + intersectionDistance * *newVelocity;
		float4 normal = (float4)( 0, 0, -1, 0 );
		float4 reflection = *newVelocity - 2.0f * DOT( *newVelocity, normal ) * normal;
		float remaining = timeStep - intersectionDistance;
		position = intersection;
		*newVelocity = reflection;
		*newPosition = intersection + remaining * damping * reflection;
	}

}



__kernel void sortPostPass(
						   __global uint2 * particleIndex,
						   __global uint  * particleIndexBack,
						   __global float4 * position,
						   __global float4 * velocity,
						   __global float4 * sortedPosition,
						   __global float4 * sortedVelocity
						   )
{
	int id = get_global_id( 0 );
	uint2 spi = particleIndex[ id ];//contains id of cell and id of particle it has sorted 
	int serialId = PI_SERIAL_ID( spi );//get a particle Index
	int cellId = PI_CELL_ID( spi );//get a cell Index
	float4 position_ = position[ serialId ];//get position by serialId
	POSITION_CELL_ID( position_ ) = (float)cellId;
	float4 velocity_ = velocity[ serialId ];
	sortedVelocity[ id ] = velocity_;//put velocity to sortedVelocity for right order according to particleIndex
	sortedPosition[ id ] = position_;//put position to sortedVelocity for right order according to particleIndex

	particleIndexBack[ serialId ] = id;
}

//=================================
// PCI SPH KERNELS BELOW
//=================================

__kernel void pcisph_computeDensity(
									 __global float2 * neighborMap,
									 float Wpoly6Coefficient,
									 float gradWspikyCoefficient,
									 float h,
									 float mass,
									 float rho0,
									 float simulationScale,
									 float stiffness,
									 __global float4 * sortedPosition,
									 __global float * pressure,
									 __global float * rho,
									 __global uint * particleIndexBack,
									 float delta									 )
{
	int id = get_global_id( 0 );
	id = particleIndexBack[id];//track selected particle (indices are not shuffled anymore)
	int idx = id * NEIGHBOR_COUNT;
	int nc=0;//neighbor counter
	float density = 0.0f;
	float r_ij2;//squared r_ij
	float hScaled = h * simulationScale;//scaled smoothing radius
	float hScaled2 = hScaled*hScaled;//squared scaled smoothing radius
	float hScaled6 = hScaled2*hScaled2*hScaled2;
	float2 nm;
	int real_nc = 0;

	do// gather density contribution from all neighbors (if they exist)
	{
		if( NEIGHBOR_MAP_ID( neighborMap[ idx + nc ] ) != NO_PARTICLE_ID )
		{
			r_ij2= NEIGHBOR_MAP_DISTANCE( neighborMap[ idx + nc ] );	// distance is already scaled here
			r_ij2 *= r_ij2;
			density += (hScaled2-r_ij2)*(hScaled2-r_ij2)*(hScaled2-r_ij2);
			real_nc++;
		}

	}while( ++nc < NEIGHBOR_COUNT );
	
	//if(density==0.f) density = hScaled2*hScaled2*hScaled2;
	if(density<hScaled6) density = hScaled6;

	density *= mass*Wpoly6Coefficient; // since all particles are same fluid type, factor this out to here
	rho[ id ] = density; 		
}
/*
float4 calcBoundaryForceAcceleration(float4 position,
									 float4 velocity,
									 float xmin,
									 float xmax,
									 float ymin,
									 float ymax,
									 float zmin,
									 float zmax,
									 float h,
									 float simulationScale)
{
    float4 acceleration = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
	float hScaled = h*simulationScale;
    float dist_iw; //i-th particle to wall distance
    float diff;
    float boundaryStiffness = 2000.0f;
    float boundaryDampening = 256.0f;


	float value = 32; //value
    
	//-----------------------------------------------
	if ( ( diff = (position[0]-xmin)*simulationScale ) < hScaled)
    {
        float4 norm =  (float4)( 1.f, 0.f, 0.f, 0.f );
        float adj = boundartStiffness * diff - boundaryDampening * DOT(norm, velocity);
        acceleration +=  norm * adj;
    }

	if ( ( diff = (xmax-position[0])*simulationScale ) < hScaled)
    {
        float4 norm =  (float4)(-1.f, 0.f, 0.f, 0.f );
        float adj = boundartStiffness * diff - boundaryDampening * DOT(norm, velocity);
        acceleration +=  norm * adj;
    }
	//-----------------------------------------------
	if ( ( diff = (position[1]-ymin)*simulationScale ) < hScaled)
    {
        float4 norm =  (float4)( 0.f, 1.f, 0.f, 0.f );
        float adj = boundartStiffness * diff - boundaryDampening * DOT(norm, velocity);
        acceleration +=  norm * adj;
    }

	if ( ( diff = (ymax-position[1])*simulationScale ) < hScaled)
    {
        float4 norm =  (float4)( 0.f,-1.f, 0.f, 0.f );
        float adj = boundartStiffness * diff - boundaryDampening * DOT(norm, velocity);
        acceleration +=  norm * adj;
    }
	//-----------------------------------------------
	if ( ( diff = (position[2]-zmin)*simulationScale ) < hScaled)
    {
        float4 norm =  (float4)( 0.f, 0.f, 1.f, 0.f );
        float adj = boundartStiffness * diff - boundaryDampening * DOT(norm, velocity);
        acceleration +=  norm * adj;
    }

	if ( ( diff = (zmax-position[2])*simulationScale ) < hScaled)
    {
        float4 norm =  (float4)( 0.f, 0.f,-1.f, 0.f );
        float adj = boundartStiffness * diff - boundaryDampening * DOT(norm, velocity);
        acceleration +=  norm * adj;
    }
	//-----------------------------------------------

    return acceleration;
}
*/

__kernel void pcisph_computeForcesAndInitPressure(
								  __global float2 * neighborMap,
								  __global float * rho,
								  __global float  * pressure,
								  __global float4 * sortedPosition,
								  __global float4 * sortedVelocity,
								  __global float4 * acceleration,
								  __global uint * particleIndexBack,
								  float Wpoly6Coefficient,
								  float del2WviscosityCoefficient,
								  float h,
								  float mass,
								  float mu,
								  float simulationScale,
								  float gravity_x,
								  float gravity_y,
								  float gravity_z,
								  __global float4 * position,
								  __global uint2 * particleIndex
								  //int ELASTIC_CONNECTIONS_COUNT,
								  //__global float4 * elasticConnectionsData
								  )
{
	int id = get_global_id( 0 );
	id = particleIndexBack[id];//track selected particle (indices are not shuffled anymore)
	int id_source_particle = PI_SERIAL_ID( particleIndex[id] );
	if((int)(position[ id_source_particle ].w) == BOUNDARY_PARTICLE){
		//FOR BOUNDARY OARTICLE WE SHOULDN'T COMPUDE ACCELERATION BECAUSE THEY DON'T MOVE SOMEWHERE
		acceleration[ id ] = (float4)(0.0f, 0.0f, 0.0f, 0.0f );
		acceleration[ PARTICLE_COUNT+id ] = (float4)(0.0f, 0.0f, 0.0f, 0.0f );
		pressure[id] = 0.f;//initialize pressure with 0
		return;
	}
		
	int idx = id * NEIGHBOR_COUNT;
	float hScaled = h * simulationScale;
	float hScaled2 = hScaled*hScaled;//29aug_A.Palyanov

	float4 acceleration_i;// = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
	float2 nm;
	float r_ij;
	int nc = 0;//neighbor counter
	int jd;
	float4 sum = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
	float4 vi,vj;
	float rho_i,rho_j;
	float4 accel_surfTensForce = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
	//float4 normalVector = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
	//float  nV_length;
	//int neighbor_cnt = 0;


	do{
		if( (jd = NEIGHBOR_MAP_ID(neighborMap[ idx + nc])) != NO_PARTICLE_ID )
		{
			r_ij = NEIGHBOR_MAP_DISTANCE( neighborMap[ idx + nc] );

			if(r_ij<hScaled)
			{
				//neighbor_cnt++;
				rho_i = rho[id];
				rho_j = rho[jd];
				vi = sortedVelocity[id];
				vj = sortedVelocity[jd];
				sum += (sortedVelocity[jd]-sortedVelocity[id])*(hScaled-r_ij)/rho[jd];
				//29aug_A.Palyanov_start_block
				// M.Beckner & M.Teschner / Weakly compressible SPH for free surface flows. 2007.
				//normalVector += sortedPosition[id]-sortedPosition[jd];
				//	-0.3f * (sortedPosition[id]-sortedPosition[jd])*simulationScale * pow(hScaled2/2/*-r_ij*r_ij*/,3);
				//29aug_A.Palyanov_end_block
				//0.09 for experiments with water drops
				//-0.0133
				
				accel_surfTensForce += -0.0013*Wpoly6Coefficient*pow(hScaled2/2/*-r_ij*r_ij*/,3)*(sortedPosition[id]-sortedPosition[jd])*simulationScale;
			}
		}
		
	}while(  ++nc < NEIGHBOR_COUNT );

	accel_surfTensForce.w = 0.f;

	/*
	if(neighbor_cnt)
	{
		normalVector[3] = 0;
		normalVector /= (float)neighbor_cnt;
		nV_length = SQRT( normalVector[0]*normalVector[0]+
						  normalVector[1]*normalVector[1]+
						  normalVector[2]*normalVector[2]);

		if(nV_length>h/10.f)
		{
			accel_surfTensForce = -(normalVector/nV_length)*25;
		}
		}*///TODO REMOVE THIS BLOCK

	float viscosity = 0.3;//0.5f;//0.1f

	sum *= mass*viscosity*del2WviscosityCoefficient/rho[id];

	// apply external forces
	acceleration_i = sum;
	//acceleration_i = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );//sum;
	
	//acceleration_i += calcBoundaryForceAcceleration(sortedPosition[id],sortedVelocity[id],xmin,xmax,ymin,ymax,zmin,zmax,h,simulationScale);
	
	acceleration_i += (float4)( gravity_x, gravity_y, gravity_z, 0.0f );

	acceleration_i +=  accel_surfTensForce; //29aug_A.Palyanov

	acceleration[ id ] = acceleration_i; 

	// 1st half of 'acceleration' array is used to store acceleration corresponding to gravity, visc. force etc.
	acceleration[ PARTICLE_COUNT+id ] = (float4)(0.0f, 0.0f, 0.0f, 0.0f );
	// 2nd half of 'acceleration' array is used to store pressure force

	pressure[id] = 0.f;//initialize pressure with 0

}


__kernel void pcisph_computeElasticForces(
								  __global float2 * neighborMap,
								  __global float4 * sortedPosition,
								  __global float4 * sortedVelocity,
								  __global float4 * acceleration,
								  __global uint * particleIndexBack,
								  float h,
								  float mass,
								  float simulationScale,
								  int ELASTIC_CONNECTIONS_COUNT,
								  __global float4 * elasticConnectionsData
								  )
{
	int index = get_global_id( 0 );
	//return;
	if(index>=ELASTIC_CONNECTIONS_COUNT) {
		return;
	}
	int id = (int)elasticConnectionsData[index].x;
	int jd = (int)elasticConnectionsData[index].y;

	id = particleIndexBack[id];
	jd = particleIndexBack[jd];

	float r_ij_equilibrium = elasticConnectionsData[index].z;

	float4 vect_r_ij = (sortedPosition[id] - sortedPosition[jd])*simulationScale;
	vect_r_ij.w = 0;
	float r_ij = SQRT(vect_r_ij.x*vect_r_ij.x+vect_r_ij.y*vect_r_ij.y+vect_r_ij.z*vect_r_ij.z);
	float delta_r_ij = r_ij - r_ij_equilibrium;
	//if(r_ij!=0.f) vect_r_ij /= r_ij;//normalize TODO: delete this

	float k = 90000.f;// k - coefficient of elasticity
	if(r_ij!=0.f){
		//In this Place can be simultaneously writing to one place in memory
		//Because id can appear more that one time
		//elasticConnectionsData can have more that one cells with 
		//id jd1 rij0 0
		//id jd2 rij1 0
		acceleration[ id ] += -(vect_r_ij/r_ij) * delta_r_ij * k;
	}
		
	float4 centerOfMassVelocity = (sortedVelocity[id]+sortedVelocity[jd])/2.f;
	float4 velocity_i_cm = sortedVelocity[id] - centerOfMassVelocity;
	//float4 velocity_j_cm = sortedVelocity[jd] - centerOfMassVelocity;// always will be equal to (-velocity_i_cm)
	float v_i_cm_length = sqrt(velocity_i_cm.x*velocity_i_cm.x+velocity_i_cm.y*velocity_i_cm.y+velocity_i_cm.z*velocity_i_cm.z);

	if((v_i_cm_length!=0)&&(r_ij!=0))
	{
		float damping_coeff = 0.5f;
		velocity_i_cm.w = 0.f;
		//vect_r_ij.w = 0;
		float4 proj_v_i_cm_on_r_ij = vect_r_ij * DOT(velocity_i_cm,vect_r_ij)/(v_i_cm_length*r_ij);
		float check = DOT(velocity_i_cm,vect_r_ij)/(v_i_cm_length*r_ij);//TODO remove it 

		//In this Place can be simultaneously writing to one place in memory
		//Because id and jd can appear more that one time
		//id jd1 rij0 0
		//id1 jd1 rij1 0
		sortedVelocity[id] -= proj_v_i_cm_on_r_ij * damping_coeff;
		sortedVelocity[jd] += proj_v_i_cm_on_r_ij * damping_coeff;
	}



	return;
}
//boundaryHandling
float w_icb(float r0, float r_ib){
	if((r0-r_ib) >= 0)
		return (r0-r_ib)/r0;
	else
		return 0.f;
}
void calculateBoundaryParticleAffect(
									   int id, 
									   float r0, 
									   __global float2 * neighborMap,
									   __global uint * particleIndexBack,
									   __global uint2 * particleIndex,
									   __global float4 * position,
									   __global float4 * velocity,
									   float4 * pos_,
									   bool tangVel,
									   float4 * vel
									   )
{
	//track selected particle (indices are not shuffled anymore)
	int idx = id * NEIGHBOR_COUNT;
	int id_source_particle, nc = 0;
	float4 n_ci = (float4)(0.f,0.f,0.f,0.f); 
	float4 n_b;
	float w_icb_summ = 0.f, w_icb_second_summ = 0.f,w_icb_current;
	float x_ib_norm;
	float4 dist;
	float val,x_ib_norma;
	int jd;
	//float8 returnedValue;
	do// gather density contribution from all neighbors (if they exist)
	{
		if( (jd = NEIGHBOR_MAP_ID( neighborMap[ idx + nc ])) != NO_PARTICLE_ID )
		{
			//jd = particleIndexBack[jd];
			id_source_particle = PI_SERIAL_ID( particleIndex[jd] );
			if((int)position[id_source_particle].w == 3){
				//dist = pos_ - position[id_source_particle];
				x_ib_norma = ((*pos_).x - position[id_source_particle].x) * ((*pos_).x - position[id_source_particle].x);
				x_ib_norma += ((*pos_).y - position[id_source_particle].y) * ((*pos_).y - position[id_source_particle].y);
				x_ib_norma += ((*pos_).z - position[id_source_particle].z) * ((*pos_).z - position[id_source_particle].z);
				x_ib_norma = SQRT(x_ib_norma);
				w_icb_current = w_icb(r0,x_ib_norma);
				n_b = velocity[id_source_particle];
				n_ci += n_b * w_icb_current;
				w_icb_summ += w_icb_current;
				w_icb_second_summ += w_icb_current * (r0 - x_ib_norma);
			}
		}
	}while( ++nc < NEIGHBOR_COUNT );
	val = DOT(n_ci,n_ci);
	if(val != 0){
		val = sqrt(val);
		//n_ci = n_ci / val;
		dist = n_ci/val * w_icb_second_summ * 1/w_icb_summ;
		(*pos_).x += dist.x;
		(*pos_).y += dist.y;
		(*pos_).z += dist.z;
		if(tangVel){
			float eps = 0.99f;//eps <= 1.0
			float vel_n_len = n_ci.x * (*vel).x + n_ci.y * (*vel).y + n_ci.z * (*vel).z; 
			if(vel_n_len < 0){
				(*vel).x -= n_ci.x * vel_n_len;
				(*vel).y -= n_ci.y * vel_n_len;
				(*vel).z -= n_ci.z * vel_n_len;
				(*vel) = (*vel) * eps;
			}
		}
	}
}
//
__kernel void pcisph_predictPositions(
						__global float4 * acceleration,
						__global float4 * sortedPosition,
						__global float4 * sortedVelocity,
						__global uint2 * particleIndex,
						__global uint * particleIndexBack,
						float gravity_x,
						float gravity_y,
						float gravity_z,
						float simulationScaleInv,
						float timeStep,
						float xmin,
						float xmax,
						float ymin,
						float ymax,
						float zmin,
						float zmax,
						float damping,
						__global float4 * position,
						__global float4 * velocity,
						float r0,
						__global float2 * neighborMap
						)
{
	int id = get_global_id( 0 );
	id = particleIndexBack[id];
	int id_source_particle = PI_SERIAL_ID( particleIndex[id] );
	float4 position_ = sortedPosition[ id ];
	if((int)(position[ id_source_particle ].w) == 3){
		sortedPosition[PARTICLE_COUNT+id] = position_;//this line was missing (absent) and caused serions errors int further program behavior
		return;
	}
	float4 acceleration_ = acceleration[ id ] + acceleration[ PARTICLE_COUNT+id ];
	float4 velocity_ = sortedVelocity[ id ];
	// Semi-implicit Euler integration 
	float4 newVelocity_ = velocity_ + timeStep * acceleration_; //newVelocity_.w = 0.f;
	float posTimeStep = timeStep * simulationScaleInv;			
	float4 newPosition_ = position_ + posTimeStep * newVelocity_; //newPosition_.w = 0.f;


	//handleBoundaryConditions( position_, &newVelocity_, posTimeStep, &newPosition_,
	//	xmin, xmax, ymin, ymax, zmin, zmax, damping ); //TODO: remove this

	//sortedVelocity[id] = newVelocity_;// sorted position, as well as velocity, 
	calculateBoundaryParticleAffect(id,r0,neighborMap,particleIndexBack,particleIndex,position,velocity,&newPosition_,false, &newVelocity_);
	sortedPosition[PARTICLE_COUNT+id] = newPosition_;// in current version sortedPosition array has double size, 
													 // PARTICLE_COUNT*2, to store both x(t) and x*(t+1)
}


__kernel void pcisph_predictDensity(
									 __global float2 * neighborMap,
									 __global uint * particleIndexBack,
									 float Wpoly6Coefficient,
									 float gradWspikyCoefficient,
									 float h,
									 float mass,
									 float rho0,
									 float simulationScale,
									 float stiffness,
									 __global float4 * sortedPosition,
									 __global float * pressure,
									 __global float * rho,
									 float delta
									 )
{
	int id = get_global_id( 0 );
	id = particleIndexBack[id];//track selected particle (indices are not shuffled anymore)
	int idx = id * NEIGHBOR_COUNT;
	int nc=0;//neighbor counter
	float density = 0.0f;
	float4 r_ij;
	float r_ij2;//squared r_ij
	float hScaled = h * simulationScale;//scaled smoothing radius
	float hScaled2 = hScaled*hScaled;//squared scaled smoothing radius
	float hScaled6 = hScaled2*hScaled2*hScaled2;
	//float2 nm;
	int jd;
	int real_nc = 0;

	do// gather density contribution from all neighbors (if they exist)
	{
		if( (jd = NEIGHBOR_MAP_ID( neighborMap[ idx + nc ])) != NO_PARTICLE_ID )
		{
			r_ij = sortedPosition[PARTICLE_COUNT+id]-sortedPosition[PARTICLE_COUNT+jd];
			r_ij2 = (r_ij.x*r_ij.x+r_ij.y*r_ij.y+r_ij.z*r_ij.z)*simulationScale*simulationScale;

			if(r_ij2<hScaled2)
			{
				density += (hScaled2-r_ij2)*(hScaled2-r_ij2)*(hScaled2-r_ij2);
				real_nc++;
			}
		}

	}while( ++nc < NEIGHBOR_COUNT );
	
	//if(density==0.f) 
	if(density<hScaled6)
	{
		//density += hScaled6;
		density = hScaled6;
	}


	density *= mass*Wpoly6Coefficient; // since all particles are same fluid type, factor this out to here
	rho[ PARTICLE_COUNT+id ] = density; 
}


__kernel void pcisph_correctPressure(
									 __global float2 * neighborMap,
									  __global uint * particleIndexBack,
									 float Wpoly6Coefficient,
									 float gradWspikyCoefficient,
									 float h,
									 float mass,
									 float rho0,
									 float simulationScale,
									 float stiffness,
									 __global float4 * sortedPosition,
									 __global float * pressure,
									 __global float * rho,
									 float delta,
									 __global float4 * position,
									 __global uint2 * particleIndex
									 )
{
	
	int id = get_global_id( 0 );
	id = particleIndexBack[id];//track selected particle (indices are not shuffled anymore)
	/*int id_source_particle = PI_SERIAL_ID( particleIndex[id] );
	if((int)(position[ id_source_particle ].w) == 3)
		return;*/
	int idx = id * NEIGHBOR_COUNT;
	int nc = 0;// neigbor counter
	float rho_err;
	float p_corr;


	rho_err = rho[PARTICLE_COUNT+id] - rho0;
	p_corr = rho_err*delta;
	if(p_corr < 0) p_corr = 0;//non-negative pressure
	pressure[ id ] += p_corr;

	//just to view the variable value;
	//p_corr = pressure[ id ];
	//p_corr = 0.f;
}


__kernel void pcisph_computePressureForceAcceleration(
								  __global float2 * neighborMap,
								  __global float * pressure,
								  __global float * rho,
								  __global float4 * sortedPosition,
								  __global float4 * sortedVelocity,
								  __global uint * particleIndexBack,
								  float CFLLimit,
								  float del2WviscosityCoefficient,
								  float gradWspikyCoefficient,
								  float h,
								  float mass,
								  float mu,
								  float simulationScale,
								  __global float4 * acceleration,
								  float rho0,
								  __global float4 * position,
								  __global uint2 * particleIndex
								  )
{
	int id = get_global_id( 0 );
	id = particleIndexBack[id];//track selected particle (indices are not mixed anymore)
	int id_source_particle = PI_SERIAL_ID( particleIndex[id] );
	if((int)(position[ id_source_particle ].w) == 3){//TODO change all number const to symbol
		acceleration[ PARTICLE_COUNT+id ] = 0.f;
		return;
	}
	int idx = id * NEIGHBOR_COUNT;
	float hScaled = h * simulationScale;

	//float4 position_i = sortedPosition[ id ];
	//float4 velocity_i = sortedVelocity[ id ];
	float pressure_i  = pressure[ id ]; 
	float rho_i		  = rho[ PARTICLE_COUNT+id ];

	float4 result = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
	//float2 nm;

	int nc=0;
	float4 gradW_ij;
	//float4 rj,rij,ri = sortedPosition[id];//x_i(t) // (not x_i*(t+1)) 
	//ri.w = 0.f;
	float r_ij,rho_err;
	float4 vr_ij;
	int jd;
	float value;
	int real_neighbors = 0;
	int total_neighbors = 0;
	
	do
	{
		if( (jd = NEIGHBOR_MAP_ID( neighborMap[ idx + nc ])) != NO_PARTICLE_ID)
		{
			r_ij = NEIGHBOR_MAP_DISTANCE( neighborMap[ idx + nc] );

			if(r_ij<hScaled)
			{
				value = -(hScaled-r_ij)*(hScaled-r_ij)*0.5f*(pressure[id]+pressure[jd])/rho[PARTICLE_COUNT+jd];
				/*2*///value = -(hScaled-r_ij)*(hScaled-r_ij)*( pressure[id]/(rho[PARTICLE_COUNT+id]*rho[PARTICLE_COUNT+id])
				/*2*///										+pressure[jd]/(rho[PARTICLE_COUNT+id]*rho[PARTICLE_COUNT+id]) );
				vr_ij = (sortedPosition[id]-sortedPosition[jd])*simulationScale; vr_ij.w = 0;
				result += value*vr_ij/r_ij;
				//result = result;

				// according to formula (3.3) in B. Solenthaler's dissertation "Incompressible Fluid Simulation and Advanced Surface Handling with SPH"
				// http://www.ifi.uzh.ch/pax/uploads/pdf/publication/1299/Solenthaler.pdf

				//result += -mass*(pressure[id]/(rho[PARTICLE_COUNT+id]*rho[PARTICLE_COUNT+id])
				//			   + pressure[jd]/(rho[PARTICLE_COUNT+jd]*rho[PARTICLE_COUNT+jd]))*gradW_ij;
				real_neighbors++;
			}

			total_neighbors++;
		}

	}while( ++nc < NEIGHBOR_COUNT );

	/*1*/result *= mass*gradWspikyCoefficient/rho[PARTICLE_COUNT+id];
	/*2*///result *= mass*gradWspikyCoefficient;
	//
	//result = -2.f*mass*pressure[id]*sum_gradW/(rho0*rho0);
	//result.w = 0.0f;
	acceleration[ PARTICLE_COUNT+id ] = result; // pressureForceAcceleration "=" or "+=" ???

}

__kernel void pcisph_integrate(
						__global float4 * acceleration,
						__global float4 * sortedPosition,
						__global float4 * sortedVelocity,
						__global uint2 * particleIndex,
						__global uint * particleIndexBack,
						float gravity_x,
						float gravity_y,
						float gravity_z,
						float simulationScaleInv,
						float timeStep,
						float xmin,
						float xmax,
						float ymin,
						float ymax,
						float zmin,
						float zmax,
						float damping,
						__global float4 * position,
						__global float4 * velocity,
						__global float * rho,
						float r0,
						__global float2 * neighborMap
						)
{
	int id = get_global_id( 0 ); if(id>=PARTICLE_COUNT) return;
	id = particleIndexBack[id]; if(id>=PARTICLE_COUNT) return;
	int id_source_particle = PI_SERIAL_ID( particleIndex[id] );
	float4 position_ = sortedPosition[ id ];

	if((int)(position[ id_source_particle ].w) == 3){
		return;
	}
	float4  accelOld = acceleration[ id ];
	float4  accelT = acceleration[ PARTICLE_COUNT+id ];
	float4 acceleration_ = acceleration[ id ] + acceleration[ PARTICLE_COUNT+id ]; acceleration_.w = 0.f;
	float4 velocity_ = sortedVelocity[ id ];
	//float4 acceleration_Fp;
	//acceleration_Fp = acceleration[ PARTICLE_COUNT+id ]; acceleration_Fp.w = 0.f;
	//acceleration_ += acceleration[ PARTICLE_COUNT+id ]; 
	//float speedLimit = 100.f;

	// Semi-implicit Euler integration 
	//float nV_length;
	//float vel_limit = 100.f;
	float4 newVelocity_ = velocity_ + timeStep * acceleration_  ; //newVelocity_.w = 0.f;
		
	//newVelocity_[3] = 0;
	//nV_length = SQRT(DOT(newVelocity_,newVelocity_));
	//if(nV_length>vel_limit)
	//	newVelocity_ = newVelocity_*vel_limit/nV_length;

	float posTimeStep = timeStep * simulationScaleInv;			
	float4 newPosition_ = position_ + posTimeStep * newVelocity_; //newPosition_.w = 0.f;


	//handleBoundaryConditions( position_, &newVelocity_, posTimeStep, &newPosition_,
	//	xmin, xmax, ymin, ymax, zmin, zmax, damping );

	//if(mode==0)
	/*{
		//sortedVelocity[PARTICLE_COUNT+id] = newVelocity_;// sorted position, as well as velocity, in current version has double size, PARTICLE_COUNT*2, to store both x(t) and x*(t+1)
		sortedPosition[PARTICLE_COUNT+id] = newPosition_;
	}*/
	//else//if mode==1

	// in Chao Fang version here is also acceleration 'speed limit' applied

	//newPosition_.w = rho[id];

	//rho[ id_source_particle ] = rho[id];

	if(newPosition_.x<xmin) newPosition_.x = xmin;//A.Palyanov 30.08.2012
	if(newPosition_.y<ymin) newPosition_.y = ymin;//A.Palyanov 30.08.2012
	if(newPosition_.z<zmin) newPosition_.z = zmin;//A.Palyanov 30.08.2012
	if(newPosition_.x>xmax-0.000001f) newPosition_.x = xmax-0.000001f;//A.Palyanov 30.08.2012
	if(newPosition_.y>ymax-0.000001f) newPosition_.y = ymax-0.000001f;//A.Palyanov 30.08.2012
	if(newPosition_.z>zmax-0.000001f) newPosition_.z = zmax-0.000001f;//A.Palyanov 30.08.2012
	// better replace 0.0000001 with smoothingRadius*0.001 or smth like this to keep this

	float particleType = position[ id_source_particle ].w;
	newVelocity_ = (velocity_ + newVelocity_) * 0.5f ;
	calculateBoundaryParticleAffect(id,r0,neighborMap,particleIndexBack,particleIndex,position,velocity,&newPosition_, true, &newVelocity_);
	velocity[ id_source_particle ] = newVelocity_;//newVelocity_;
	position[ id_source_particle ] = newPosition_;
	position[ id_source_particle ].w = particleType;
	// position[0..2] stores x,y,z; position[3] - for particle type
}
