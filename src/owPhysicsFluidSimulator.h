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
	//TODO helper functions delete after fix!!
	float * getElasticConnections() { return elasticConnections; };
	//
	double simulationStep();
private:
	owOpenCLSolver * ocl_solver;
	float * positionBuffer;
	float * velocityBuffer;
	float * elasticConnections;
	//Helper buffers
	float * densityBuffer;
	unsigned int * particleIndexBuffer;
	float * accelerationBuffer;//TODO REMOVE after fixing
	owHelper * helper;
};

#endif //OW_PHYSICS_SIMULATOR_H