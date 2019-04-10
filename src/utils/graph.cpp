//
// Created by serg on 07.04.19.
//
#include <cmath>
#include <csignal>
#include <iostream>
#include <sstream>

#include "graph.h"
#include "VectorMath.h"

using sibernetic::graphics::graph;
using sibernetic::graphics::g_config;

float simulationScale = 0.5f;
inline void beginWinCoords() {
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(0.0f, (GLfloat)glutGet(GLUT_WINDOW_HEIGHT) - 10, 0.0f);
	glScalef(8.0f, -1.0f, 1.0f);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, glutGet(GLUT_WINDOW_WIDTH), 0, glutGet(GLUT_WINDOW_HEIGHT), -1, 1);

	glMatrixMode(GL_MODELVIEW);
}

inline void endWinCoords() {
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}
void glPrint(float x, float y, const char *s, void *font) {
	glRasterPos2f((GLfloat)x, (GLfloat)y);
	int len = (int)strlen(s);
	for (int i = 0; i < len; ++i) {
		glutBitmapCharacter(font, s[i]);
	}
}
void glPrint3D(float x, float y, float z, const char *s, void *font) {
	glRasterPos3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
	int len = (int)strlen(s);
	for (int i = 0; i < len; ++i) {
		glutBitmapCharacter(font, s[i]);
	}
}

void graph::draw_scene(){
	//       [7]----[6]
	//      / |     /|
	//    [3]----[2] |
	//     | [4]--|-[5]
	//     | /    | /
	//    [0]----[1]
	Vector3D vcenter(0, 0, 0);
	Vector3D vbox[8];
	float s_v = 1 / (simulationScale); // = 1 m in simulation
	float order = 0;
	while (s_v >= 1) {
		s_v /= 10;
		if (s_v < 1) {
			s_v *= 10;
			break;
		}
		++order;
	}
	vbox[0] = Vector3D(config->xmin, config->ymin, config->zmin);
	vbox[1] = Vector3D(config->xmax, config->ymin, config->zmin);
	vbox[2] = Vector3D(config->xmax, config->ymax, config->zmin);
	vbox[3] = Vector3D(config->xmin, config->ymax, config->zmin);
	vbox[4] = Vector3D(config->xmin, config->ymin, config->zmax);
	vbox[5] = Vector3D(config->xmax, config->ymin, config->zmax);
	vbox[6] = Vector3D(config->xmax, config->ymax, config->zmax);
	vbox[7] = Vector3D(config->xmin, config->ymax, config->zmax);

	// Display user interface if enabled
//	bool displayInfos = true;
//	if (displayInfos) {
//		glDisable(GL_DEPTH_TEST);
//		glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ZERO); // invert color
//		glEnable(GL_BLEND);
//		renderInfo(0, 0);
//		glDisable(GL_BLEND);
//		glEnable(GL_DEPTH_TEST);
//	}
	glBegin(GL_LINES);
	sc *= 10;
	glColor3ub(255, 0, 0);
	glVertex3d(vcenter.x, vcenter.y, vcenter.z);
	glVertex3d(vcenter.x + sc, vcenter.y, vcenter.z);
	glColor3ub(0, 255, 0);
	glVertex3d(vcenter.x, vcenter.y, vcenter.z);
	glVertex3d(vcenter.x, vcenter.y + sc, vcenter.z);
	glColor3ub(0, 0, 255);
	glVertex3d(vcenter.x, vcenter.y, vcenter.z);
	glVertex3d(vcenter.x, vcenter.y, vcenter.z + sc);
	sc /= 10;
	vcenter = Vector3D(-(config->xmin + config->xmax) / 2,
	                   -(config->ymin + config->ymax) / 2,
	                   -(config->zmin + config->zmax) / 2);
	vcenter *= sc;
	Vector3D v1, v2, v3, v4, v5, v6, v7, v8;
	v1 = Vector3D(-config->xmax / 2, -config->ymax / 2,
	              -config->zmax / 2) * sc;
	v2 = Vector3D(config->xmax / 2, -config->ymax / 2,
	              -config->zmax / 2) * sc;
	v3 = Vector3D(config->xmax / 2, config->ymax / 2,
	              -config->zmax / 2) * sc;
	v4 = Vector3D(-config->xmax / 2, config->ymax / 2,
	              -config->zmax / 2) * sc;
	v5 = Vector3D(-config->xmax / 2, -config->ymax / 2,
	              config->zmax / 2) * sc;
	v6 = Vector3D(config->xmax / 2, -config->ymax / 2,
	              config->zmax / 2) * sc;
	v7 = Vector3D(config->xmax / 2, config->ymax / 2,
	              config->zmax / 2) * sc;
	v8 = Vector3D(-config->xmax / 2, config->ymax / 2,
	              config->zmax / 2) * sc;
	glColor3ub(255, 255, 255); // yellow
	glVertex3d(v1.x, v1.y, v1.z);
	glVertex3d(v2.x, v2.y, v2.z);

	glColor3ub(255, 255, 255); // yellow
	glVertex3d(v2.x, v2.y, v2.z);
	glVertex3d(v3.x, v3.y, v3.z);

	glVertex3d(v3.x, v3.y, v3.z);
	glVertex3d(v4.x, v4.y, v4.z);

	glVertex3d(v4.x, v4.y, v4.z); // glColor3ub(0,255,0);//green
	glVertex3d(v1.x, v1.y, v1.z);

	// glColor3ub(0,0,255);//blue
	glVertex3d(v1.x, v1.y, v1.z); // glColor3ub(255,255,0);//yellow
	glVertex3d(v5.x, v5.y, v5.z);

	glVertex3d(v2.x, v2.y, v2.z);
	glVertex3d(v6.x, v6.y, v6.z);

	glVertex3d(v3.x, v3.y, v3.z);
	glVertex3d(v7.x, v7.y, v7.z);

	glVertex3d(v4.x, v4.y, v4.z);
	glVertex3d(v8.x, v8.y, v8.z);

	glVertex3d(v5.x, v5.y, v5.z);
	glVertex3d(v6.x, v6.y, v6.z);

	glVertex3d(v6.x, v6.y, v6.z);
	glVertex3d(v7.x, v7.y, v7.z);

	glVertex3d(v7.x, v7.y, v7.z);
	glVertex3d(v8.x + s_v * sc, v8.y, v8.z);

	glVertex3d(v8.x, v8.y, v8.z);
	glVertex3d(v5.x, v5.y, v5.z);
	glEnd();
	//
	glBegin(GL_LINES);
	glColor3ub(0, 0, 0); // black

	Vector3D v_s = Vector3D(-config->xmax / 2 + s_v, config->ymax / 2,
	                        config->zmax / 2) * sc;
	glVertex3d(v_s.x, v_s.y, v_s.z);
	glVertex3d(v_s.x, v_s.y - 0.5f * sc, v_s.z);
	glLineWidth((GLfloat)10.0);
	glVertex3d(v8.x, v8.y, v8.z);
	glVertex3d(v_s.x, v_s.y, v_s.z);

	glEnd();
	glLineWidth((GLfloat)1.0);
	void *m_font = GLUT_BITMAP_8_BY_13;
	std::stringstream ss;
	std::string s;
	ss << order;
	s = "1E-" + ss.str() + "m";
	glPrint3D((float)v8.x + 0.4f * sc, (float)v8.y - 2.f * sc, (float)v8.z, "0",
	          m_font);
	glPrint3D((float)v_s.x, (float)v_s.y - 2.f * sc, (float)v_s.z, s.c_str(),
	          m_font);
	ss.str("");
	while (v_s.x < config->xmax / 2 * sc) {
		v_s.x += s_v * sc;
		if (v_s.x < config->xmax / 2 * sc) {
			glBegin(GL_LINES);
			glVertex3d(v_s.x, v_s.y, v_s.z);
			glVertex3d(v_s.x, v_s.y - 0.5f * sc, v_s.z);
			glEnd();
		}
	}
}

void graph::display(){
	// Update Scene if not paused
	int i, j, k;
	int err_coord_cnt = 0;
	double calculationTime;
	double renderTime;
	void *m_font = GLUT_BITMAP_8_BY_13;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	draw_scene();
	draw_model();
//	float dc, rho;
//	glLineWidth((GLfloat)1.0);
	glutSwapBuffers();
}

void graph::draw_model() {
	glBegin(GL_POINTS);
	glPointSize(1.3f * sqrt(sc / 0.025f));
	for(auto p :model->get_particles()){
		if (1/*p.type != 3*/) {
			glColor4f(0, 0, 0, 1.0f); // color of elastic particles
			glPointSize(1.3f * sqrt(sc / 0.025f));
			glVertex3f((p.pos[0] - config->xmax / 2) * sc,
					   (p.pos[1] - config->ymax / 2) * sc,
			           (p.pos[2] - config->zmax / 2) * sc);
		}
	}
	glEnd();
}

void graph::respond_mouse_callback(int button, int state, int x, int y) {
	if (button == GLUT_RIGHT_BUTTON)
		button_state = 3;
	if (button == GLUT_LEFT_BUTTON)
		button_state = 1;
	int mods = glutGetModifiers();
	if (mods & GLUT_ACTIVE_CTRL) {
		button_state = 2;
	}
	if (state == GLUT_UP)
		button_state = 0;
	old_x = x;
	old_y = y;
	if (button == 3) // mouse wheel up
	{
		sc *= 1.1f; // Zoom in
	}
	if (button == 4) // mouse wheel down
	{
		sc /= 1.1f; // Zoom out
	}
}

// GLUT callback
// called on mouse movement

void graph::mouse_motion_callback(int x, int y) {
	float dx, dy;
	dy = static_cast<float>(y - old_y);
	dx = static_cast<float>(x - old_x);

	if (button_state == 1) {
		camera_rot[0] += dy / 5.0f;
		camera_rot[1] += dx / 5.0f;
	}
	if (button_state == 3) {

		camera_trans[0] += dx / 100.0f;
		camera_trans[1] -= dy / 100.0f;
	}
	if (button_state == 2) {
		camera_trans[0] += dx / 100.0f;
		camera_trans[1] -= dy / 100.0f;
	}
	old_x = x;
	old_y = y;
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	for (int c = 0; c < 3; ++c) {
		camera_trans_lag[c] += (camera_trans[c] - camera_trans_lag[c]);
		camera_rot_lag[c] += (camera_rot[c] - camera_rot_lag[c]);
	}
	glTranslatef(camera_trans_lag[0], camera_trans_lag[1], camera_trans_lag[2]);
	glRotatef(camera_rot_lag[0], 1.0, 0.0, 0.0);
	glRotatef(camera_rot_lag[1], 0.0, 1.0, 0.0);
}

void graph::key_pressed_callback(unsigned char key, int x, int y) {
	switch (key) {
		case '1':
			//config->setConfigFileName("demo1");
			//helper->refreshTime();
//			fluid_simulation->reset();
//			sPause = false;
			break;
		case '2':
			// owHelper::configFileName = "demo2";
			//config->setConfigFileName("demo2");
			//helper->refreshTime();
			//fluid_simulation->reset();
//			sPause = false;
			break;
		case '\033': // Escape quits
		case 'Q':    // Q quits
		case 'q':    // q quits
//			cleanupSimulation();
			// break;
			exit(EXIT_SUCCESS);
		case ' ':
//			sPause = !sPause;
			std::cout << "\nSimulation Is Paused" << std::endl;
			break;
		case 's':
//			fluid_simulation->makeSnapshot();
			break;
		case 'r': // reset simulation
			//helper->refreshTime();
			try {
//				fluid_simulation->reset();
			} catch (std::runtime_error &ex) {
//				cleanupSimulation();
				std::cout << "ERROR: " << ex.what() << std::endl;
				exit(EXIT_FAILURE);
			}
			break;
		case 'i':
//			showInfo = !showInfo;
			break;
	}
	glutPostRedisplay();
}

void graph::idle() { glutPostRedisplay(); }

void graph::timer(int value) {
	// Re-register for next callback
	glutPostRedisplay();
	glutTimerFunc(TIMER_INTERVAL * 0, &timer, 0);
}

void graph::resize_callback(GLsizei width, GLsizei height) {
	if (height == 0) {
		height = 1;
	}
	if (width == 0) {
		width = 1;
	}

	glViewport(0, 0, width, height); // Set view area
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	float aspectRatio = (GLfloat)width / (GLfloat)height;
	if (aspectRatio > 1.f)
		glFrustum(-1 * aspectRatio, 1 * aspectRatio, -1, 1, 3, 45);
	else
		glFrustum(-1, 1, -1 / aspectRatio, 1 / aspectRatio, 3, 45);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	for (int c = 0; c < 3; ++c) {
		camera_trans_lag[c] += (camera_trans[c] - camera_trans_lag[c]);
		camera_rot_lag[c] += (camera_rot[c] - camera_rot_lag[c]);
	}
	glTranslatef(camera_trans_lag[0], camera_trans_lag[1], camera_trans_lag[2]);
	glRotatef(camera_rot_lag[0], 1.0, 0.0, 0.0);
	glRotatef(camera_rot_lag[1], 0.0, 1.0, 0.0);
}

void graph::init() {
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
	glEnable(GL_AUTO_NORMAL);
	float ambient[4] = {1.0, 1.0, 1.0, 1};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
	glClearColor(0.7f, 0.7f, 0.7f, 1.0f); // background color
	glClearDepth(1.0f);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

int graph::old_x = 0;
int graph::old_y = 0; // Used for mouse event
float graph::camera_trans[] = {0, 0, -8.f};
float graph::camera_rot[] = {60, -90, 0}; // camera rotation settings at start
float graph::camera_trans_lag[] = {0, 0, -8.f};
float graph::camera_rot_lag[] = {0, 0, 0};
int graph::button_state = 0;
float graph::sc = 0.025f; // 0.0145;//0.045;//0.07
double graph::total_time = 0;
int graph::frames_counter = 0;
std::shared_ptr<sibernetic::model::sph_model<float>> graph::model = nullptr;

g_config * graph::config = new g_config({0.0,0.0,0.0, 100.0, 100.0, 100.0});

void graph::run(int argc, char **argv){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(1200, 800);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("SIBERNETIC (2011-2017) by Andrey Palyanov and Sergey "
	                 "Khayrulin. Build from 28/10/2017 sources (development "
	                 "branch)");
	glutIdleFunc(&idle);
	init();
	glutDisplayFunc(&display);
	glutReshapeFunc(&resize_callback);
	glutMouseFunc(&respond_mouse_callback);
	glutMotionFunc(&mouse_motion_callback); // process movement in case if the mouse is clicked,
	glutKeyboardFunc(&key_pressed_callback);
	glutTimerFunc(TIMER_INTERVAL * 0, &timer, 0);
	glutMainLoop();
}
