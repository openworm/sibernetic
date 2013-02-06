#ifndef OW_PHYSICS_SIMULATOR_H
#define OW_PHYSICS_SIMULATOR_H

#include "owPhysicsConstant.h"
#include "owHelper.h"
#include "owOpenCLSolver.h"

class owPhysicsFluidSimulator
{
public:
	owPhysicsFluidSimulator(void);
	owPhysicsFluidSimulator(owHelper * helper);
	~owPhysicsFluidSimulator(void);
	float * getPositionBuffer() { return positionBuffer; };
	float * getVelocityBuffer() { return velocityBuffer; };
	float * getDensityBuffer() { ocl_solver->read_density_b( densityBuffer ); return densityBuffer; };
	unsigned int * getParticleIndexBuffer() { ocl_solver->read_particleIndex_b( particleIndexBuffer ); return particleIndexBuffer; };
	double simulationStep();
	int get_numOfLiquidP() { return numOfLiquidP; };
	int get_numOfElasticP() { return numOfElasticP; };
	int get_numOfBoundaryP() { return numOfBoundaryP; };
	
private:
	int numOfLiquidP;
	int numOfElasticP;
	int numOfBoundaryP;
	owOpenCLSolver * ocl_solver;
	float * positionBuffer;
	float * velocityBuffer;
	float * elasticConnections;
	//Helper buffers
	float * densityBuffer;
	unsigned int * particleIndexBuffer;
	owHelper * helper;
};

#endif //OW_PHYSICS_SIMULATOR_H