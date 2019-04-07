//
// Created by serg on 07.04.19.
//

#ifndef SIBERNETIC_GRAPH_H
#define SIBERNETIC_GRAPH_H

#if defined(_WIN32) || defined (_WIN64)
#include <GL/glew.h>
	#include <GL/wglew.h>
#else
	#include <cstring>
#endif

#if defined(__APPLE__) || defined(MACOSX)
	#include <GLUT/glut.h>
#else
	#include <GL/freeglut.h>
#endif


class graph {
public:
	void draw_scene();
	void run();
};


#endif //SIBERNETIC_GRAPH_H
