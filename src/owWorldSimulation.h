#ifndef OW_WORLD_SIMULATION_H
#define OW_WORLD_SIMULATION_H

#if defined(_WIN32) || defined (_WIN64)
	#include <GL/glew.h>
	#include <GL/wglew.h>
#else
	#include <string.h>
#endif
#if defined(__APPLE__) || defined(MACOSX)
    #include <GLUT/glut.h>
#else
    #include <GL/freeglut.h>
#endif

#include "owPhysicsFluidSimulator.h"
#include "VectorMath.h"

#define TIMER_INTERVAL 30  //this is the interval between calls to timer func (in milliseconds)
#define ROTATION_STEP_ANGLE 1      //this is the step angle that the mesh will rotate every SOME_INTERVAL milliseconds
void run(int argc, char** argv, const bool with_graphics = true);
#endif //OW_WORLD_SIMULATION_H
