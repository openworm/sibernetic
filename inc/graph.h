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

#include <memory>
#include "sph_model.hpp"
#include "solver_container.hpp"

#define TIMER_INTERVAL 30  			//this is the interval between calls to timer func (in milliseconds)
#define ROTATION_STEP_ANGLE 1 //this is the step angle that the mesh will rotate every SOME_INTERVAL milliseconds

namespace sibernetic {
	namespace graphics {
		struct g_config{
			float xmin;
			float ymin;
			float zmin;
			float xmax;
			float ymax;
			float zmax;
		};
		class graph {
		public:
			static void run(int, char **);
			static g_config * config;
			static std::shared_ptr<sibernetic::model::sph_model<float>> model;
			static sibernetic::solver::solver_container<float> *s_container;
		private:
			static int old_x;
			static int old_y; // Used for mouse event
			static float camera_trans[];
			static float camera_rot[]; // camera rotation settings at start
			static float camera_trans_lag[];
			static float camera_rot_lag[];
			static int button_state;
			static float sc; // 0.0145;//0.045;//0.07
			static double total_time;
			static int frames_counter;

			static void draw_partition();
			static void draw_model();
			static void draw_scene();
			static void init();
			static void resize_callback(GLsizei width, GLsizei height);
			static void display();
			static void idle();
			static void timer(int);
			static void key_pressed_callback(unsigned char key, int x, int y);
			static void mouse_motion_callback(int x, int y);
			static void respond_mouse_callback(int button, int state, int x, int y);
		};
	}
}
#endif //SIBERNETIC_GRAPH_H
