#ifndef OW_WORLD_SIMULATION_H
#define OW_WORLD_SIMULATION_H


#if defined(_WIN32) || defined (_WIN64)
	#pragma comment( lib, "gl//glut32.lib" )
	#include <windows.h>
	#include "../gl/glut.h"
#else
	#include <string.h>
	#include <GL/glut.h>
#endif

#include "owPhysicsFluidSimulator.h"
#include "VectorMath.h"

#define TIMER_INTERVAL 30  //this is the interval between calls to timer func (in milliseconds)
#define ROTATION_STEP_ANGLE 1      //this is the step angle that the mesh will rotate every SOME_INTERVAL milliseconds
void run(int argc, char** argv, const bool with_graphics = true);
#endif //OW_WORLD_SIMULATION_H