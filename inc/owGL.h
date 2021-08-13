#ifndef __OW_GL_H__
#define __OW_GL_H__
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

#endif
