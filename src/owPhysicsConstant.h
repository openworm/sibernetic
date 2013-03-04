#ifndef OW_PHYSICS_CONSTANT_H
#define OW_PHYSICS_CONSTANT_H

#include "owOpenCLConstant.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415927f
#endif
//Sizes of the box containing simulated 'world'
//Sizes choice is realized this way because it should be proportional to smoothing radius h
#define XMIN 0
#define XMAX 120.24f//120//200//50//100
#define YMIN 0
#define YMAX 80.16f//80//330//40//80
#define ZMIN 0
#define ZMAX 182.03f//180//10//30

const float rho0 = 1000.0f;
const float stiffness = 0.75f;
const float h = 3.34f;
const float hashGridCellSize = 2.0f * h;
const float hashGridCellSizeInv = 1.0f / hashGridCellSize;
const float mass = 0.0003f;//0.0003
const float simulationScale = 0.004f;
const float simulationScaleInv = 1.0f / simulationScale;
const float mu = 10.0f;//why this value? Dynamic viscosity of water at 25 C = 0.89e-3 Pa*s
const float timeStep = 0.001f;//0.0005f;//0.0042f;// ATTENTION you should remember about time step: if it is larger 0.001 it can lead to 'explosion' of elastic matter objects
const float CFLLimit = 100.0f;
const int NK = NEIGHBOR_COUNT * PARTICLE_COUNT;
const float damping = 0.75f;
const float r0 = 0.5f * h; // distance between two boundary particle
const float beta = timeStep*timeStep*mass*mass*2/(rho0*rho0);// B. Solenthaler's dissertation, formula 3.6 (end of page 30)
const float betaInv = 1.f/beta;
const float Wpoly6Coefficient = 315.0f / ( 64.0f * M_PI * pow( h * simulationScale, 9.0f ) );
const float gradWspikyCoefficient= -45.0f / ( M_PI * pow( h * simulationScale, 6.0f ) );
const float del2WviscosityCoefficient = -gradWspikyCoefficient;
const float gravity_x = 0.0f;
const float gravity_y = -9.8f;
const float gravity_z = 0.0f;
extern const float delta;
const int maxIteration = 3;
const int ELASTIC_CONNECTIONS_COUNT = 0;

#endif // #ifndef OW_PHYSICS_CONSTANT_H