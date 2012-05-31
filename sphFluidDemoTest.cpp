
#ifndef NDEBUG

#include <assert.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.hpp>
#endif


#define EPSILON 0.1
#define MAX_ERRORS 10
#define NO_PARTICLE_ID -1

#define NEIGHBOR_COUNT 32
#define PARTICLE_COUNT ( 32 * 1024 )

static float * readField( 
						 int particleId,
						 float * _field
						 ){
							 return _field + 4 * particleId;
}


bool testBeginClearBuffers( cl::CommandQueue queue ){
	return true;
}

bool testEndClearBuffers( cl::CommandQueue queue, int _NK, cl::Buffer neighborMap ){
	bool result = true;
	int messageCount = 0;

	float * _neighborMap = new float[ ( _NK ) * 2 ];
	queue.enqueueReadBuffer( neighborMap, CL_TRUE, 0, ( _NK ) * sizeof( float ) * 2, _neighborMap );
	queue.finish();

	for( int i = 0; i < PARTICLE_COUNT; ++i ){
		float neighborId = _neighborMap[ i * 2 ];
		if( !( neighborId == -1 )){
			std::cerr << "Error: failed to initialize neighborMap at location " << i << std::endl;
			result = false;
		}
		messageCount += ( result == false );
		if( messageCount > MAX_ERRORS ) break;
	}//for
	delete [] _neighborMap;

	return result;
}



// test kernel HashParticles

bool testBeginHashParticles( cl::CommandQueue queue ){
	bool result = true;
	return result;
}

bool testEndHashParticles( cl::CommandQueue queue, float hashGridCellSize, int gridCellCount,
						  int gridCellsX, int gridCellsY, int gridCellsZ,
						  float xmin, float ymin, float zmin,
						  cl::Buffer position, cl::Buffer particleIndex ){
	bool result = true;
	int messageCount = 0;

	float * _position = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( position, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _position );
	queue.finish();

	unsigned int * _particleIndex = new unsigned int[ PARTICLE_COUNT * 2 ];
	queue.enqueueReadBuffer( particleIndex, CL_TRUE, 0, PARTICLE_COUNT * sizeof( unsigned int ) * 2, _particleIndex );
	queue.finish();

	for( int i = 0; i < PARTICLE_COUNT; ++i ){

		int cellId = _particleIndex[ i * 2 ];
		int particleId = _particleIndex[ i * 2 + 1 ];
		if ( cellId < 0 || cellId >= gridCellCount ){
			std::cerr << "Error: bad cellId " << cellId << " at particle " << particleId << std::endl;
			result = false;
			messageCount++;
		}

		float * p = readField( particleId, _position );
		int ix = (int)(( p[ 0 ] - xmin ) / hashGridCellSize ) % gridCellsX;
		int iy = (int)(( p[ 1 ] - ymin ) / hashGridCellSize ) % gridCellsY;
		int iz = (int)(( p[ 2 ] - zmin ) / hashGridCellSize ) % gridCellsZ;
		if( ix < 0 ) ix += gridCellsX;
		if( iy < 0 ) iy += gridCellsY;
		if( iz < 0 ) iz += gridCellsZ;
		assert( ix >= 0 && iy >= 0 && iz >= 0 );
		int computedCellId = ix + iy * gridCellsX + iz * gridCellsX * gridCellsY;
		computedCellId = computedCellId & 0xffff; // truncate to low 16 bits
		if( computedCellId != cellId ){
			std::cerr << "Error: inconsistent cellId " << cellId << " at particle " << particleId
				<< " expected " << computedCellId << std::endl;
			result = false;
			messageCount++;
		}

		if( messageCount > MAX_ERRORS ) break;
	}//for

	delete [] _position;
	delete [] _particleIndex;
	return result;
}

// test kernel Sort

bool testBeginSort( cl::CommandQueue queue, cl::Buffer particleIndex, int gridCellCount ){
	int messageCount = 0;
	bool result = true;

	unsigned int * _particleIndex = new unsigned int[ PARTICLE_COUNT * 2 ];
	queue.enqueueReadBuffer( particleIndex, CL_TRUE, 0, PARTICLE_COUNT * sizeof( unsigned int ) * 2, _particleIndex );
	queue.finish();

	for( int i = 0; i < PARTICLE_COUNT; ++i ){
		int cellId = _particleIndex[ i * 2 ];
		int particleId = _particleIndex[ i * 2 + 1 ];
		if( cellId < 0 || cellId > gridCellCount ){
			std::cerr << "Error: bad cellId " << cellId << " input data to sort, element " << i << std::endl;
			result = false;
			if( ++messageCount > MAX_ERRORS ) break;
			continue;
		}
		if( particleId != i ){
			std::cerr << "Error: bad particleId " << particleId << " input data to sort, element " << i << std::endl;
			result = false;
			if( ++messageCount > MAX_ERRORS ) break;
			continue;
		}
	}//for

	delete [] _particleIndex;
	return result;
}

bool testEndSort( cl::CommandQueue queue, cl::Buffer particleIndex ){
	bool result = true;
	int messageCount = 0;

	unsigned int * _particleIndex = new unsigned int[ PARTICLE_COUNT * 2 ];
	queue.enqueueReadBuffer( particleIndex, CL_TRUE, 0, PARTICLE_COUNT * sizeof( unsigned int ) * 2, _particleIndex );
	queue.finish();

	int lastCellId = 0;
	for( int i = 0; i < PARTICLE_COUNT; ++i ){
		int cellId = _particleIndex[ i * 2 ];
		int particleId = _particleIndex[ i * 2 + 1 ];
		if ( cellId < lastCellId || cellId < 0 ){
			std::cerr << "Error: data out of order " << cellId << " at position " << i << std::endl;
			result = false;
		}
		lastCellId = cellId;
		messageCount += ( result == false );
		if( messageCount > MAX_ERRORS ) break;
	}//for
	delete [] _particleIndex;
	return result;
}


// test kernel SortPostPass

bool testBeginSortPostPass( cl::CommandQueue queue ){
	return true;
}


bool testEndSortPostPass( cl::CommandQueue queue, cl::Buffer particleIndex, cl::Buffer position,
						 cl::Buffer velocity, cl::Buffer sortedPosition, cl::Buffer sortedVelocity,
						 int gridCellCount ){
	bool result = true;
	int messageCount = 0;

	unsigned int * _particleIndex = new unsigned int[ PARTICLE_COUNT * 2 ];
	queue.enqueueReadBuffer( particleIndex, CL_TRUE, 0, PARTICLE_COUNT * sizeof( unsigned int ) * 2, _particleIndex );
	queue.finish();

	float * _position = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( position, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _position );
	queue.finish();

	float * _velocity = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( velocity, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _velocity );
	queue.finish();

	float * _sortedPosition = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( sortedPosition, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _sortedPosition );
	queue.finish();

	float * _sortedVelocity = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( sortedVelocity, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _sortedVelocity );
	queue.finish();

	for( int i = 0; i < PARTICLE_COUNT; ++i ){
		unsigned int * spi = _particleIndex + 2 * i;	
		int cellId = spi[ 0 ];
		int serialId = spi[ 1 ];
		if( cellId < 0 || cellId > gridCellCount
			|| serialId < 0 || serialId > PARTICLE_COUNT ){
				std::cerr << "Error: impossible results from sort, particle " << i 
					<< " cell " << cellId << " serialId " << serialId << std::endl;
				result = false;
				if( ++messageCount > MAX_ERRORS ) break;
				continue;
		}
		float * oldPosition = readField( serialId, _position );
		float * newPosition = readField( i, _sortedPosition );
		if( oldPosition[ 0 ] != newPosition[ 0 ] ||
			oldPosition[ 1 ] != newPosition[ 1 ] ||
			oldPosition[ 2 ] != newPosition[ 2 ] ){
				std::cerr << "Error: particle " << i << " position does not match sorting element "
					<< serialId << std::endl;
				result = false;
				if( ++messageCount > MAX_ERRORS ) break;
				continue;
		}
		float * oldVelocity = readField( serialId, _velocity );
		float * newVelocity = readField( i, _sortedVelocity );
		if( oldVelocity[ 0 ] != newVelocity[ 0 ] ||
			oldVelocity[ 1 ] != newVelocity[ 1 ] ||
			oldVelocity[ 2 ] != newVelocity[ 2 ] ){
				std::cerr << "Error: particle " << i << " velocity does not match sorting element "
					<< serialId << std::endl;
				result = false;
				if( ++messageCount > MAX_ERRORS ) break;
				continue;
		}
	}//for

	delete [] _particleIndex;
	delete [] _position;
	delete [] _velocity;
	delete [] _sortedPosition;
	delete [] _sortedVelocity;
	return result;
}


// test kernel Index

bool testBeginIndexx( cl::CommandQueue queue ){
	return true;
}

bool testEndIndexx( cl::CommandQueue queue, int gridCellCount, cl::Buffer particleIndex,
				  cl::Buffer gridCellIndex ){
	bool result = true;
	int messageCount = 0;

	unsigned int * _particleIndex = new unsigned int[ PARTICLE_COUNT * 2 ];
	queue.enqueueReadBuffer( particleIndex, CL_TRUE, 0, PARTICLE_COUNT * sizeof( unsigned int ) * 2, _particleIndex );
	queue.finish();

	unsigned int * _gridCellIndex = new unsigned int[ ( (gridCellCount) + (1) ) * 1 ];
	queue.enqueueReadBuffer( gridCellIndex, CL_TRUE, 0, ( (gridCellCount) + (1) ) * sizeof( unsigned int ) * 1, _gridCellIndex );
	queue.finish();

	// test that all values are increasing except for NO_CELL_ID entries
	int lastParticleIndex = 0;
	for( int i = 0; i < gridCellCount + 1; ++i ){
		int particleIndex = _gridCellIndex[ i ];
		if( particleIndex < 0 ) continue;
		if( particleIndex < lastParticleIndex ){
			std::cerr << "Error: particle index out of order " << particleIndex << " at grid cell " << i << std::endl;
			result = false;
		}
		lastParticleIndex = particleIndex;
		messageCount += ( result == false );
		if( messageCount > MAX_ERRORS ) break;
	}//for

	// compute the expected result and compare
	int lastCellId = -1;
	for( int i = 0; i < PARTICLE_COUNT; ++i ){
		int cellId = _particleIndex[ i * 2 ];
		int particleId = _particleIndex[ i * 2 + 1 ];
		if( cellId < 0 || cellId > gridCellCount 
			|| particleId < 0 || particleId > PARTICLE_COUNT ){
				std::cerr << "Error: impossible particle index values, particle " << i
					<< " cell " << cellId << " particleId " << particleId << std::endl;
				result = false;
				if( ++messageCount > MAX_ERRORS ) break;
				continue;
		}
		if( cellId > lastCellId ){
			int indexValue = _gridCellIndex[ cellId ];
			if( indexValue >= 0 && indexValue != i ){
				std::cerr << "Error: grid cell index mismatch at cellId " << cellId << " = " 
					<< _gridCellIndex[ cellId ] << " expected " << i << std::endl;
				result = false;
			}
			lastCellId = cellId;
		}
		messageCount += ( result == false );
		if( messageCount > MAX_ERRORS ) break;
	}//for

	// verify the last+1 grid cell index
	if( _gridCellIndex[ gridCellCount ] != PARTICLE_COUNT ){
		std::cerr << "Error: last+1 grid cell index is wrong " 
			<< _gridCellIndex[ gridCellCount ] << " should be " << PARTICLE_COUNT << std::endl;
		result = false;
	}

	delete [] _particleIndex;
	delete [] _gridCellIndex;
	return result;
}

// test kernel IndexPostPass

bool testBeginIndexPostPass( cl::CommandQueue queue ){
	return true;
}


bool testEndIndexPostPass( cl::CommandQueue queue, int gridCellCount, cl::Buffer gridCellIndex,
						  cl::Buffer gridCellIndexFixedUp, cl::Buffer particleIndex ){
	bool result = true;
	int messageCount = 0;

	unsigned int * _gridCellIndex = new unsigned int[ ( (gridCellCount) + (1) ) * 1 ];
	queue.enqueueReadBuffer( gridCellIndex, CL_TRUE, 0, ( (gridCellCount) + (1) ) * sizeof( unsigned int ) * 1, _gridCellIndex );
	queue.finish();

	unsigned int * _gridCellIndexFixedUp = new unsigned int[ ( (gridCellCount) + (1) ) * 1 ];
	queue.enqueueReadBuffer( gridCellIndexFixedUp, CL_TRUE, 0, ( (gridCellCount) + (1) ) * sizeof( unsigned int ) * 1, _gridCellIndexFixedUp );
	queue.finish();

	// test that all values are increasing except for NO_CELL_ID entries
	int * _particlesPerCell = new int[ gridCellCount ];
	int lastParticleIndex = 0;
	for( int i = 0; i < gridCellCount + 1; ++i ){
		int particleIndex = _gridCellIndexFixedUp[ i ];
		if( particleIndex < lastParticleIndex ){
			std::cerr << "Error: gridCellIndexFixedUp out of order " << particleIndex << " at grid cell " << i << std::endl;
			result = false;
		}
		if( i ){
			_particlesPerCell[ i - 1 ] = particleIndex - lastParticleIndex;
		}
		lastParticleIndex = particleIndex;
		messageCount += ( result == false );
		if( messageCount > MAX_ERRORS ) break;
	}//for

	float sumParticlesPerCell = 0.0f;
	for( int i = 0; i < gridCellCount; ++i ){
		sumParticlesPerCell += _particlesPerCell[ i ];
	}//for
	float meanParticlesPerCell = sumParticlesPerCell / (float)gridCellCount;
	float totalVariance = 0.0f;
	for( int i = 0; i < gridCellCount; ++i ){
		float varianceThisCell = ( _particlesPerCell[ i ] - meanParticlesPerCell ) * ( _particlesPerCell[ i ] - meanParticlesPerCell );
		totalVariance += varianceThisCell;
	}//for
	float stdDev = sqrt( totalVariance / (float)gridCellCount );

	std::cout << "Particles per cell, mean " << meanParticlesPerCell << " standard deviation " << stdDev << std::endl;

	float expectedParticlesPerCell = (float)PARTICLE_COUNT / (float)gridCellCount;
	if( meanParticlesPerCell != expectedParticlesPerCell ){
		std::cerr << "Error: wrong mean particles per cell " << meanParticlesPerCell 
			<< " should be " << expectedParticlesPerCell << std::endl;
		result = false;
	}

	float maxExpectedStandardDeviation = (float)expectedParticlesPerCell * 0.06f;//observed empirically 
	// to do - predict expected std dev based on grid size and particle count
	if( stdDev > maxExpectedStandardDeviation ){
		std::cerr << "Warning: nonequilibrium standard deviation " << stdDev 
			<< " of particles per cell, equilbrium value would be " 
			<< maxExpectedStandardDeviation << std::endl;
	}

	delete [] _particlesPerCell;

	// compute the expected result and compare

	unsigned int * _particleIndex = new unsigned int[ PARTICLE_COUNT * 2 ];
	queue.enqueueReadBuffer( particleIndex, CL_TRUE, 0, PARTICLE_COUNT * sizeof( unsigned int ) * 2, _particleIndex );
	queue.finish();

	int lastCellId = -1;
	for( int i = 0; i < PARTICLE_COUNT; ++i ){
		int cellId = _particleIndex[ i * 2 ];
		int particleId = _particleIndex[ i * 2 + 1 ];
		if( cellId > lastCellId ){
			int indexValue = _gridCellIndex[ cellId ];
			if( indexValue != i ){
				std::cerr << "Error: grid cell index mismatch at cellId " << cellId << " = " 
					<< _gridCellIndex[ cellId ] << " expected " << i << std::endl;
				result = false;
			}
			lastCellId = cellId;
		}
		messageCount += ( result == false );
		if( messageCount > MAX_ERRORS ) break;
	}//for

	// verify the last+1 grid cell index
	if( _gridCellIndex[ gridCellCount ] != PARTICLE_COUNT ){
		std::cerr << "Error: last+1 grid cell index is wrong " << _gridCellIndex[ gridCellCount ] 
		<< " should be " << PARTICLE_COUNT << std::endl;
		result = false;
	}

	delete [] _particleIndex;
	delete [] _gridCellIndex;
	delete [] _gridCellIndexFixedUp;
	return result;
}

// test kernel FindNeighbors


bool testBeginFindNeighbors( cl::CommandQueue queue ){
	return true;
}



bool testEndFindNeighbors( cl::CommandQueue queue, int gridCellCount, float hashGridCellSize, float hashGridCellSizeInv,
						  cl::Buffer gridCellIndexFixedUp, cl::Buffer sortedPosition, cl::Buffer neighborMap, int _NK,
						  float simulationScale, float h ){
	bool result = true;
	int messageCount = 0;

	unsigned int * _gridCellIndexFixedUp = new unsigned int[ ( (gridCellCount) + (1) ) * 1 ];
	queue.enqueueReadBuffer( gridCellIndexFixedUp, CL_TRUE, 0, ( (gridCellCount) + (1) ) * sizeof( unsigned int ) * 1, _gridCellIndexFixedUp );
	queue.finish();

	float * _sortedPosition = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( sortedPosition, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _sortedPosition );
	queue.finish();

	float * _neighborMap = new float[ ( _NK ) * 2 ];
	queue.enqueueReadBuffer( neighborMap, CL_TRUE, 0, ( _NK ) * sizeof( float ) * 2, _neighborMap );
	queue.finish();

	// make sure neighbors are symmetric and agree with distance

	for( int particleId = 0; particleId < PARTICLE_COUNT; ++particleId ){
		float * nm = _neighborMap + particleId * NEIGHBOR_COUNT * 2;
		float * myPosition = readField( particleId, _sortedPosition );

		for( int neighborNum = 0; neighborNum < NEIGHBOR_COUNT; ++neighborNum ){
			float neighborParticleId = nm[ neighborNum * 2 ];
			if( (int)neighborParticleId == NO_PARTICLE_ID ) continue;

			// make sure this particle is listed in the neighbors neighbor map

			float *nnm = _neighborMap + (int)neighborParticleId * NEIGHBOR_COUNT * 2;
			bool isNeighbor = false;
			for( int j = 0; j < NEIGHBOR_COUNT; ++j ){
				int nnNum = (int)nnm[ j * 2 ];
				isNeighbor |= nnNum == particleId;
			}//for
			if( !isNeighbor ){
				std::cerr << "Warning: particle " << particleId << " has neighbor particle "
					<< neighborParticleId << " but is not listed in the neighbor map for that particle\n";
				messageCount++;
				continue;
			}

			// make sure the computed distance is correct

			float distance = nm[ neighborNum * 2 + 1 ];
			float * neighborPosition = readField( (int)neighborParticleId, _sortedPosition );

			float dx = myPosition[ 0 ] - neighborPosition[ 0 ];
			float dy = myPosition[ 1 ] - neighborPosition[ 1 ];
			float dz = myPosition[ 2 ] - neighborPosition[ 2 ];
			float distanceSquared = dx * dx + dy * dy + dz * dz;
			float computedDistance = sqrt( distanceSquared ) * simulationScale;

			if( fabs( 1.0f - computedDistance / distance ) > EPSILON ){
				std::cerr << "Error: particle " << particleId << " has neighbor particle "
					<< neighborParticleId << " but distance " << distance
					<< " does not agree with computed distance "
					<< computedDistance << std::endl;
				result = false;
				messageCount++;
			}	
			if( messageCount > MAX_ERRORS ) break;
		}//for
		if( messageCount > MAX_ERRORS) break;
	}//for

	float totalNeighbors = 0.0f;
	for( int particleId = 0; particleId < PARTICLE_COUNT; ++particleId ){
		float * nm = _neighborMap + particleId * NEIGHBOR_COUNT * 2;
		float * myPosition = readField( particleId, _sortedPosition );

		for( int neighborNum = 0; neighborNum < NEIGHBOR_COUNT; ++neighborNum ){
			float neighborParticleId = nm[ neighborNum * 2 ];
			if( (int)neighborParticleId == NO_PARTICLE_ID ) continue;
			totalNeighbors += 1.0f;
		}//for
	}//for
	float meanNeighbors = totalNeighbors / PARTICLE_COUNT;
	if( meanNeighbors < 0.7f * NEIGHBOR_COUNT ){
		std::cerr << "Error: unexpectedly low meanNeighbors " << meanNeighbors << std::endl;
		result = false;
	}else{
		std::cout << "mean neighbors " << meanNeighbors << " NEIGHBOR_COUNT=" << NEIGHBOR_COUNT << std::endl;
	}

	// truly randomize sampling, assume data is same from run to run
	time_t t;
	t = time( NULL );
	srand( (unsigned int)t );

	for( int i = 0; i < PARTICLE_COUNT; ++i ){
		int particleId = i;
		float * myPosition = readField( particleId, _sortedPosition );
		float * nm = _neighborMap + particleId * NEIGHBOR_COUNT * 2;
		for( int neighborParticleId = 0; neighborParticleId < PARTICLE_COUNT; ++neighborParticleId ){
			if( neighborParticleId == particleId ) continue;
			float * neighborPosition = readField( neighborParticleId, _sortedPosition );
			float dx = myPosition[ 0 ] - neighborPosition[ 0 ];
			float dy = myPosition[ 1 ] - neighborPosition[ 1 ];
			float dz = myPosition[ 2 ] - neighborPosition[ 2 ];
			float distance = sqrt( dx * dx + dy * dy + dz * dz );
			if( distance < h - EPSILON ){
				bool found = false;
				for( int i = 0; i < NEIGHBOR_COUNT; ++i ){
					float npid = nm[ i * 2 ];
					if( (int)npid == neighborParticleId ) found = true;
				}//for
				if( !found ){
					std::cerr << "Warning: particle " << particleId << " at ("
						<< myPosition[ 0 ] << ", "
						<< myPosition[ 1 ] << ", "
						<< myPosition[ 2 ] <<") did not find neighbor "
						<< neighborParticleId << " at ("
						<< neighborPosition[ 0 ] << ", "
						<< neighborPosition[ 1 ] << ", "
						<< neighborPosition[ 2 ] << ") distance " 
						<< distance << std::endl;
					messageCount++;
				}
			}
			if( messageCount > MAX_ERRORS ) break;
		}//for
		if( messageCount > MAX_ERRORS ) break;
	}//for

	delete [] _gridCellIndexFixedUp;
	delete [] _sortedPosition;
	delete [] _neighborMap;
	return result;
}


// test kernel ComputeDensityPressure

bool testBeginComputeDensityPressure( cl::CommandQueue queue ){
	return true;
}


bool testEndComputeDensityPressure( cl::CommandQueue queue, int _NK, float mass, float h,
								   float Wpoly6Coefficient, float rho0, float stiffness,
								   cl::Buffer neighborMap, cl::Buffer pressure, cl::Buffer rho,
								   cl::Buffer rhoInv, float simulationScale ){

	bool result = true;
	int messageCount = 0;

	float * _neighborMap = new float[ ( _NK ) * 2 ];
	queue.enqueueReadBuffer( neighborMap, CL_TRUE, 0, ( _NK ) * sizeof( float ) * 2, _neighborMap );
	queue.finish();

	float * _pressure = new float[ PARTICLE_COUNT * 1 ];
	queue.enqueueReadBuffer( pressure, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 1, _pressure );
	queue.finish();

	float * _rho = new float[ PARTICLE_COUNT * 1 ];
	queue.enqueueReadBuffer( rho, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 1, _rho );
	queue.finish();

	float * _rhoInv = new float[ PARTICLE_COUNT * 1 ];
	queue.enqueueReadBuffer( rhoInv, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 1, _rhoInv );
	queue.finish();

	float hScaled = h * simulationScale;

	for( int particleId = 0; particleId < PARTICLE_COUNT; ++particleId ){
		float * nm = _neighborMap + particleId * NEIGHBOR_COUNT * 2;
		float computedDensity = 0.0f;

		for( int neighborNum = 0; neighborNum < NEIGHBOR_COUNT; ++neighborNum ){
			float r = 0.0f;
			float neighborParticleId = nm[ neighborNum * 2 ];
			if( neighborParticleId == NO_PARTICLE_ID ) continue;
			float distance = nm[ neighborNum * 2 + 1 ];
			r = distance;
			float x = hScaled * hScaled - r * r;
			float Wpoly6 = x * x * x * Wpoly6Coefficient;
			float contribution = mass * Wpoly6;
			computedDensity += contribution;
		}//for

		float drho = computedDensity - rho0;
		float k = stiffness;
		float computedPressure = k * drho;

		float density = _rho[ particleId ];
		float p = _pressure[ particleId ];

		float densityError = 1.0f - density / computedDensity;
		//float densityError = fabs( density - computedDensity );
		if( densityError > EPSILON ){
			std::cerr << "Error: particle " << particleId << " expected density "
				<< computedDensity << " but kernel computed " << density 
				<< " error " << densityError << "%" << std::endl;
			result = false;
			messageCount++;
		}
		float pressureError = 1.0f - p / computedPressure;
		//float pressureError = fabs( p - computedPressure );
		if( pressureError > EPSILON ){
			std::cerr << "Error: particle " << particleId << " expected pressure "
				<< computedPressure << " but kernel computed " << p 
				<< " error " << pressureError << "%" << std::endl;
			result = false;
			messageCount++;
		}
		if( messageCount > MAX_ERRORS ) break;
	}//for

	float totalRho = 0.0;
	float totalP = 0.0;

	for( int i = 0; i < PARTICLE_COUNT; ++i ){
		totalRho += _rho[ i ];
		totalP += _pressure[ i ];
	}//for

	float meanRho = totalRho / PARTICLE_COUNT;
	float meanP = totalP / PARTICLE_COUNT;
	std::cout << "mean rho " << meanRho << " mean pressure " << meanP << std::endl;

	delete [] _neighborMap;
	delete [] _pressure;
	delete [] _rho;
	delete [] _rhoInv;

	return result;
}

// test kernel ComputeAcceleration

bool testBeginComputeAcceleration( cl::CommandQueue queue ){
	return true;
}

bool testEndComputeAcceleration( cl::CommandQueue queue, int _NK,
								float simulationScale, float gradWspikyCoefficient,
								float del2WviscosityCoefficient, float CFLLimit,
								cl::Buffer neighborMap, cl::Buffer pressure, cl::Buffer rho,
								cl::Buffer rhoInv, cl::Buffer sortedPosition, cl::Buffer sortedVelocity,
								cl::Buffer acceleration,
								float gravity_x, float gravity_y, float gravity_z,
								float mass, float h, float mu ){

	bool result = true;
	int messageCount = 0;

	float * _neighborMap = new float[ ( _NK ) * 2 ];
	queue.enqueueReadBuffer( neighborMap, CL_TRUE, 0, ( _NK ) * sizeof( float ) * 2, _neighborMap );
	queue.finish();

	float * _pressure = new float[ PARTICLE_COUNT * 1 ];
	queue.enqueueReadBuffer( pressure, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 1, _pressure );
	queue.finish();

	float * _rho = new float[ PARTICLE_COUNT * 1 ];
	queue.enqueueReadBuffer( rho, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 1, _rho );
	queue.finish();

	float * _rhoInv = new float[ PARTICLE_COUNT * 1 ];
	queue.enqueueReadBuffer( rhoInv, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 1, _rhoInv );
	queue.finish();

	float * _sortedPosition = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( sortedPosition, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _sortedPosition );
	queue.finish();

	float * _sortedVelocity = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( sortedVelocity, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _sortedVelocity );
	queue.finish();

	float * _acceleration = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( acceleration, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _acceleration );
	queue.finish();

	float gravity[ 3 ] = { gravity_x, gravity_y, gravity_z };
	float totalMagnitude = 0.0f;
	float hScaled = h * simulationScale;

	for( int particleId = 0; particleId < PARTICLE_COUNT; ++particleId ){
		float rho_i = _rho[ particleId ];
		float rho_i_inv = 1.0f / rho_i;
		float p_i = _pressure[ particleId ];
		float * nm = _neighborMap + particleId * NEIGHBOR_COUNT * 2;
		float computedAcceleration[ 3 ] = { 0.0f, 0.0f, 0.0f }; 
		float * position_i = readField( particleId, _sortedPosition );
		float * v_i = readField( particleId, _sortedVelocity );
		float gradP[ 3 ] = { 0.0f, 0.0f, 0.0f };
		float del2V[ 3 ] = { 0.0f, 0.0f, 0.0f };

		// integrate the gradP and del^2 V terms over all particles

		for( int neighborNum = 0 /* -1 */; neighborNum < NEIGHBOR_COUNT; ++neighborNum ){
			int neighborParticleId;
			float r;
			float rho_j_inv;
			float p_j;
			float * neighborPosition;
			float d[ 3 ];
			float * v_j;


			if( neighborNum < 0 ){
				r = 0.0f;
				rho_j_inv = 1.0f / rho_i;
				p_j = p_i;
				neighborPosition = position_i;
				d[ 0 ] = d[ 1 ] = d[ 2 ] = 0.0f;
				v_j = v_i;
			}else{
				neighborParticleId = (int)nm[ neighborNum * 2 ];
				if( neighborParticleId == NO_PARTICLE_ID ) continue;
				r = nm[ neighborNum * 2 + 1 ];
				float rhoP_j = _rho[ neighborParticleId ];
				rho_j_inv = 1.0f / rhoP_j;
				p_j = _pressure[ neighborParticleId ];
				neighborPosition = readField( neighborParticleId, _sortedPosition );
				d[ 0 ] = position_i[ 0 ] - neighborPosition[ 0 ];
				d[ 1 ] = position_i[ 1 ] - neighborPosition[ 1 ];
				d[ 2 ] = position_i[ 2 ] - neighborPosition[ 2 ];
				d[ 0 ] *= simulationScale;
				d[ 1 ] *= simulationScale;
				d[ 2 ] *= simulationScale;
				d[ 0 ] /= r;
				d[ 1 ] /= r;
				d[ 2 ] /= r;
				v_j =  readField( neighborParticleId, _sortedVelocity );		  
			}

			// gradP

			float x = hScaled - r;
			float gradWspiky[ 3 ];
			gradWspiky[ 0 ] = x * x * d[ 0 ] * gradWspikyCoefficient;
			gradWspiky[ 1 ] = x * x * d[ 1 ] * gradWspikyCoefficient;
			gradWspiky[ 2 ] = x * x * d[ 2 ] * gradWspikyCoefficient;
			for( int j = 0; j < 3; ++j ){
				gradP[ j ] += mass * ( p_i + p_j ) * 0.5f * rho_j_inv * gradWspiky[ j ];
			}//for

			// del^2 V

			float del2Wviscosity = ( hScaled - r ) * del2WviscosityCoefficient;
			d[ 0 ] = ( v_j[ 0 ] - v_i[ 0 ] );
			d[ 1 ] = ( v_j[ 1 ] - v_i[ 1 ] );
			d[ 2 ] = ( v_j[ 2 ] - v_i[ 2 ] );
			for( int j = 0; j < 3; ++j ){
				del2V[ j ] += mass * d[ j ] * rho_j_inv * del2Wviscosity;
			}//for
		}//for

		for( int j = 0; j < 3; ++j ){
			computedAcceleration[ j ] = rho_i_inv * ( mu * del2V[ j ] - gradP[ j ] );
		}//for

		// apply CFL limiting
		float magnitudeSquared = computedAcceleration[ 0 ] * computedAcceleration[ 0 ]
		+ computedAcceleration[ 1 ] * computedAcceleration[ 1 ]
		+ computedAcceleration[ 2 ] * computedAcceleration[ 2 ];
		float magnitude = sqrt( magnitudeSquared );
		if( magnitude > CFLLimit ){
			float scale = CFLLimit / magnitude;
			for( int j = 0; j < 3; ++j ){
				computedAcceleration[ j ] *= scale;
			}//for
		}
		totalMagnitude += magnitude;

		float * acceleration = readField( particleId, _acceleration );
		float d[ 3 ];
		d[ 0 ] = computedAcceleration[ 0 ] - acceleration[ 0 ];
		d[ 1 ] = computedAcceleration[ 1 ] - acceleration[ 1 ];
		d[ 2 ] = computedAcceleration[ 2 ] - acceleration[ 2 ];
		float l2Norm = sqrt( d[ 0 ] * d[ 0 ] + d[ 1 ] * d[ 1 ] + d[ 2 ] * d[ 2 ] );
		if( l2Norm > EPSILON ){
			std::cerr << "Error: particle " << particleId << " has unexpected acceleration ("
				<< acceleration[ 0 ] << ", " << acceleration[ 1 ] << ", " << acceleration[ 2 ]
			<< ") \texpected (" << computedAcceleration[ 0 ] << ", "
				<< computedAcceleration[ 1 ] << ", "
				<< computedAcceleration[ 2 ] << ")" << std::endl;
			result = false;
			messageCount++;
		}
		if( messageCount > MAX_ERRORS ) break;
	}//for

	float totalAcceleration[ 3 ] = { 0.0f, 0.0f, 0.0f };

	for( int i = 0; i < PARTICLE_COUNT; ++i ){
		float * a = _acceleration + i * 4;
		totalAcceleration[ 0 ] += a[ 0 ];
		totalAcceleration[ 1 ] += a[ 1 ];
		totalAcceleration[ 2 ] += a[ 2 ];
	}//for

	float meanAcceleration[ 3 ];
	for( int j = 0; j < 3; ++j ){
		meanAcceleration[ j ] = totalAcceleration[ j ] / PARTICLE_COUNT;
	}//for
	float magnitude = sqrt( 
		meanAcceleration[ 0 ] * meanAcceleration[ 0 ] +
		meanAcceleration[ 1 ] * meanAcceleration[ 1 ] +
		meanAcceleration[ 2 ] * meanAcceleration[ 2 ]
	);
	float meanMagnitude = totalMagnitude / PARTICLE_COUNT;
	std::cout << "mean acceleration (" <<
		meanAcceleration[ 0 ] << "," <<
		meanAcceleration[ 1 ] << "," <<
		meanAcceleration[ 2 ] << ") mean magnitude " << meanMagnitude << std::endl;

	delete [] _neighborMap;
	delete [] _pressure;
	delete [] _rho;
	delete [] _rhoInv;
	delete [] _sortedPosition;
	delete [] _sortedVelocity;
	delete [] _acceleration;

	return result;

}



bool
testBeginHandleBoundaryConditions( cl::CommandQueue queue )
{
	return true;
}

bool
testEndHandleBoundaryConditions( cl::CommandQueue queue, cl::Buffer acceleration, cl::Buffer sortedPosition,
								cl::Buffer sortedVelocity )
{
	float * _acceleration = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( acceleration, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _acceleration );
	queue.finish();

	float * _sortedPosition = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( sortedPosition, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _sortedPosition );
	queue.finish();

	float * _sortedVelocity = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( sortedVelocity, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _sortedVelocity );
	queue.finish();


	// to be completed

	delete [] _acceleration;
	delete [] _sortedPosition;
	delete [] _sortedVelocity;
	return true;
}


// test kernel Integrate

static float *_position;
static float *_velocity;

bool testBeginIntegrate( cl::CommandQueue queue ){
	bool result = true;
	return result;
}

bool testEndIntegrate( cl::CommandQueue queue, cl::Buffer acceleration, cl::Buffer sortedPosition,
					  cl::Buffer sortedVelocity, cl::Buffer position, cl::Buffer velocity,
					  float xmin, float ymin, float zmin, float xmax, float ymax, float zmax ){
	bool result = true;
	int messageCount = 0;

	float * _acceleration = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( acceleration, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _acceleration );
	queue.finish();

	float * _sortedPosition = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( sortedPosition, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _sortedPosition );
	queue.finish();

	float * _sortedVelocity = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( sortedVelocity, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _sortedVelocity );
	queue.finish();

	float * _position = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( position, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _position );
	queue.finish();

	float * _velocity = new float[ PARTICLE_COUNT * 4 ];
	queue.enqueueReadBuffer( velocity, CL_TRUE, 0, PARTICLE_COUNT * sizeof( float ) * 4, _velocity );
	queue.finish();


	for( int i = 0; i < PARTICLE_COUNT; ++i ){
		float * sortedPosition_ = readField( i, _sortedPosition );
		float * position_ = readField( i, _position );
		for( int j = 0; j < 3; ++j ){
			if( position_[ 0 ] < xmin || position_[ 0 ] > xmax ||
				position_[ 1 ] < ymin || position_[ 1 ] > ymax ||
				position_[ 2 ] < zmin || position_[ 2 ] > zmax ){
					std::cerr << "Error: particle " << i << " escaped the environment at ("
						<< position_[ 0 ] << "," << position_[ 1 ]
					<< "," << position_[ 2 ] << ")" << std::endl;
					result = false;
					break;
			}
		}//for
		messageCount += ( result == false );
		if( messageCount > MAX_ERRORS ) break;
	}//for

	delete [] _acceleration;
	delete [] _sortedPosition;
	delete [] _sortedVelocity;
	delete [] _position;
	delete [] _velocity;
	return result;
}


#endif//NDEBUG