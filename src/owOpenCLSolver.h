#ifndef OW_OPENCL_SOLVER_H
#define OW_OPENCL_SOLVER_H

#pragma comment( lib, "opencl.lib" )							// opencl.lib

#if defined(__APPLE__) || defined(__MACOSX)
	#include <OpenCL/cl.hpp>
	#include <OpenCL/cl_d3d10.h>
#else
	#include <CL/cl.hpp>
#endif
#include "owPhysicsConstant.h"

extern int PARTICLE_COUNT;
extern int PARTICLE_COUNT_RoundedUp;
extern int local_NDRange_size;

#if INTEL_OPENCL_DEBUG
	//#define  OPENCL_DEBUG_PROGRAM_PATH "-g -s \"C:\\Users\\Serg\\Documents\\GitHub\\Smoothed-Particle-Hydrodynamics\\src\\sphFluid.cl\"" // if you debuging with intel opencl debuger you need past here full path to you opencl program
#define  OPENCL_DEBUG_PROGRAM_PATH "-g -s \"C:\\GitHub\\Smoothed-Particle-Hydrodynamics\\src\\sphFluid.cl\"" // if you debuging with intel opencl debuger you need past here full path to you opencl program
#endif
#define OPENCL_PROGRAM_PATH "src/sphFluid.cl"

class owOpenCLSolver
{
public:
	owOpenCLSolver(const float * positionBuffer, const float * velocityBuffer, const float * elasticConnections = NULL);
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
	unsigned int _run_pcisph_integrate();
	
	void read_position_b( float * positionBuffer ) { copy_buffer_from_device( positionBuffer, position, PARTICLE_COUNT * sizeof( float ) * 4 ); };
	void read_density_b( float * densityBuffer ) { copy_buffer_from_device( densityBuffer, rho, PARTICLE_COUNT * sizeof( float ) * 1 ); }; // This need only for visualization current density of particle (graphic effect)
	void read_particleIndex_b( unsigned int * particeleIndexBuffer ) { copy_buffer_from_device( particeleIndexBuffer, particleIndex, PARTICLE_COUNT * sizeof( unsigned int ) * 2 ); }; // This need only for visualization current density of particle (graphic effect)
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
	cl::Buffer 	acceleration;				// forceAcceleration and pressureForceAcceleration
	cl::Buffer 	gridCellIndex;
	cl::Buffer 	gridCellIndexFixedUp;
	cl::Buffer 	neighborMap;
	cl::Buffer 	particleIndex;				// list of pairs [CellIndex, particleIndex]
	cl::Buffer 	particleIndexBack;			// list of indexes of particles before sort 
	cl::Buffer 	position;
	cl::Buffer 	pressure;
	cl::Buffer 	rho;						// size * 2
	cl::Buffer 	sortedPosition;				// size * 2
	cl::Buffer 	sortedVelocity;
	cl::Buffer 	velocity;
	cl::Buffer 	elasticConnectionsData;		//list of particle pairs connected with springs and rest distance between them
	// Kernels
	cl::Kernel clearBuffers;
	cl::Kernel computeAcceleration;
	cl::Kernel computeDensityPressure;
	cl::Kernel findNeighbors;
	cl::Kernel hashParticles;
	cl::Kernel indexx;
	cl::Kernel integrate;
	cl::Kernel sortPostPass;
	// Additional kernels for PCISPH and for calculation ellastic forces
	cl::Kernel preElasticMatterPass;
	cl::Kernel pcisph_computeDensity;
	cl::Kernel pcisph_computeForcesAndInitPressure;
	cl::Kernel pcisph_integrate;
	cl::Kernel pcisph_predictPositions;
	cl::Kernel pcisph_predictDensity;
	cl::Kernel pcisph_correctPressure;
	cl::Kernel pcisph_computePressureForceAcceleration;
	cl::Kernel pcisph_computeElasticForces;
};

#endif //OW_OPENCL_SOLVER_H