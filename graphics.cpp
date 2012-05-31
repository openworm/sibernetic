//WM_MOUSEWHEEL is only defined in later versions of windows. 
//To have the identifier defined you'll need to put the line 
#define _WIN32_WINDOWS 0x501

#define TIMER_INTERVAL 30  //this is the interval between calls to timr func (in milliseconds)
#define ROTATION_STEP_ANGLE 1      //this is the step angle that the mesh will rotate every SOME_INTERVAL milliseconds

#include "engine.h"
#include "VectorMath.h"
#include "sph.h"
#include <stdio.h>
#include <time.h>
Engine *engine;
bool rotate = false;
int old_x=0, old_y=0;     // Used for mouse event
float rotX = 0.0f;    // Rotate screen on x axis 
float rotY = 0.0f;    // Rotate screen on y axis
float rotZ = 0.0f;    // Rotate screen on z axis
bool lbutton = false;
float sc = 0.04;
Vector3D ort1(1,0,0),ort2(0,1,0),ort3(0,0,1);
//bool mouse_event = false; // need to reDraw
//extern int particleCount;
extern float * positionBuffer;

int frames_counter = 0;

extern int sph_fluid_main_start ();
extern void sph_fluid_main_step ();
extern void sph_fluid_main_stop ();

 
GLvoid display(GLvoid)
{
	frames_counter++;
	int c = clock();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Очищается буфер кадра и буфер глубины

	
	Vector3D vcenter(0,0,0);
	Vector3D vbox[8];

	//       [7]----[6]
	//      / |     /| 
	//    [3]----[2] | 
	//     | [4]--|-[5]   
	//     | /    | /
	//    [0]----[1]  

	vbox[0] = Vector3D(XMIN,YMIN,ZMIN);
	vbox[1] = Vector3D(XMAX,YMIN,ZMIN);
	vbox[2] = Vector3D(XMAX,YMAX,ZMIN);
	vbox[3] = Vector3D(XMIN,YMAX,ZMIN);
	vbox[4] = Vector3D(XMIN,YMIN,ZMAX);
	vbox[5] = Vector3D(XMAX,YMIN,ZMAX);
	vbox[6] = Vector3D(XMAX,YMAX,ZMAX);
	vbox[7] = Vector3D(XMIN,YMAX,ZMAX);
	//sph_fluid_main_step();
	glBegin(GL_LINES);

	sc *=10;

	glColor3ub(255, 0, 0);
	glVertex3d(vcenter.x,vcenter.y,vcenter.z);
	glVertex3d(vcenter.x+sc,vcenter.y,vcenter.z);

	glColor3ub(0,255, 0);
	glVertex3d(vcenter.x,vcenter.y,vcenter.z);
	glVertex3d(vcenter.x,vcenter.y+sc,vcenter.z);

	glColor3ub(0, 0, 255);
	glVertex3d(vcenter.x,vcenter.y,vcenter.z);
	glVertex3d(vcenter.x,vcenter.y,vcenter.z+sc);

	sc /=10;

	vcenter = Vector3D(-(XMIN+XMAX)/2,-(YMIN+YMAX)/2,-(ZMIN+ZMAX)/2);
	vcenter *= sc;

	Vector3D v1,v2,v3,v4,v5,v6,v7,v8;

	v1 = Vector3D( -XMAX/2, -YMAX/2, -ZMAX/2)*sc;
	v2 = Vector3D(  XMAX/2, -YMAX/2, -ZMAX/2)*sc;
	v3 = Vector3D(  XMAX/2,  YMAX/2, -ZMAX/2)*sc;
	v4 = Vector3D( -XMAX/2,  YMAX/2, -ZMAX/2)*sc;
	v5 = Vector3D( -XMAX/2, -YMAX/2,  ZMAX/2)*sc;
	v6 = Vector3D(  XMAX/2, -YMAX/2,  ZMAX/2)*sc;
	v7 = Vector3D(  XMAX/2,  YMAX/2,  ZMAX/2)*sc;
	v8 = Vector3D( -XMAX/2,  YMAX/2,  ZMAX/2)*sc;

	glColor3ub(255,0,0);//red
	glVertex3d(v1.x,v1.y,v1.z); glColor3ub(255,255,0);//yellow
	glVertex3d(v2.x,v2.y,v2.z);

	glVertex3d(v2.x,v2.y,v2.z);
	glVertex3d(v3.x,v3.y,v3.z);

	glVertex3d(v3.x,v3.y,v3.z);
	glVertex3d(v4.x,v4.y,v4.z);

	glVertex3d(v4.x,v4.y,v4.z); //glColor3ub(0,255,0);//green
	glVertex3d(v1.x,v1.y,v1.z);

	//glColor3ub(0,0,255);//blue
	glVertex3d(v1.x,v1.y,v1.z); //glColor3ub(255,255,0);//yellow
	glVertex3d(v5.x,v5.y,v5.z);

	glVertex3d(v2.x,v2.y,v2.z);
	glVertex3d(v6.x,v6.y,v6.z);

	glVertex3d(v3.x,v3.y,v3.z);
	glVertex3d(v7.x,v7.y,v7.z);

	glVertex3d(v4.x,v4.y,v4.z);
	glVertex3d(v8.x,v8.y,v8.z);

	glVertex3d(v5.x,v5.y,v5.z);
	glVertex3d(v6.x,v6.y,v6.z);

	glVertex3d(v6.x,v6.y,v6.z);
	glVertex3d(v7.x,v7.y,v7.z);

	glVertex3d(v7.x,v7.y,v7.z);
	glVertex3d(v8.x,v8.y,v8.z);

	glVertex3d(v8.x,v8.y,v8.z);
	glVertex3d(v5.x,v5.y,v5.z);
	
	glEnd();
/*
    GLfloat material_diffuse[] = {1.0, 1.0, 1.0, 1.0};
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);
    // установка источников света
    // направленный источник света

    GLfloat light0_diffuse[] = {0.7, 0.7, 1.0};
    GLfloat light0_direction[] = {0.0, 0.0, 1.0, 0.0};
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_direction);

*/

	//glPopMatrix();	

	glColor3ub(75, 135, 195);

	Vector3D v;

	//int x = positionBuffer[0];
	//float distrib[301];
	//for(int j=0;j<301;j++) distrib[j] = 0;
	//glPushMatrix();
	glBegin(GL_POINTS);
	glColor3f(1.0,1.0,1.0);
	glPointSize(1.f);
	//glEnable(GL_POINT_SMOOTH);
	for(int i = 0; i<PARTICLE_COUNT; i ++)
	{
		//GLfloat vert[] = {0.f,0.f,0.f};
		glPushMatrix();
		//glVertex3fv(vert);
		//glTranslated( (positionBuffer[i*4]-XMAX/2)*sc , (positionBuffer[i*4+1]-YMAX/2)*sc, (positionBuffer[i*4+2]-ZMAX/2)*sc );
		//distrib[(int)(positionBuffer[i*4]*300/XMAX)]+=1.f;
		//glutWireSphere( 1.0, 8, 8 );
		glVertex3f((positionBuffer[i*4]-XMAX/2)*sc , (positionBuffer[i*4+1]-YMAX/2)*sc, (positionBuffer[i*4+2]-ZMAX/2)*sc );
//		glutSolidSphere( 0.3*sc, 4, 2 );

		glPopMatrix();
	}
	glEnd();
	//glPopMatrix();
	//printf("\ntime:%d",clock()-c);
	c = clock();
	/*FILE *f = fopen("distrib.txt","wt");
	for(int j=0;j<301;j++) fprintf(f,"%.3f\n",distrib[j]);
	fclose(f);*/


    //engine->Draw();

	/**///glPopMatrix();

    glutSwapBuffers();
	//printf("\ntime:%d",clock()-c);
}
 
GLvoid reshape(GLsizei width, GLsizei height)
{
      engine->Resize(width, height);

}

void respond_mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON)
	{
		lbutton = true;
	}
	else
		lbutton = false;

	old_x=x;
	old_y=y;

	if (button == 3)// mouse wheel up
    {
        sc *= 1.1;// Zoom in
    }
    else
	if (button == 4)// mouse wheel down
    {
        sc /= 1.1;// Zoom out
    }


  // Respond to mouse button presses.
  // If button1 pressed, mark this state so we know in motion function.

	/*
	if (button == GLUT_LEFT_BUTTON)
    {
      //g_bButton1Down = (state == GLUT_DOWN) ? TRUE : FALSE;
      //g_yClick = y - 3 * g_fViewDistance;
    }*/

	/*
	if( (button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN) )
        rotate = true;
    else
        rotate = false;
	*/

	/*
            if ((mousePressed == 0))    // If left mouse button is pressed
            {
                X = (x - old_x) / 15;       // I did divide by 15 to adjust 
                                            // for a nice translation 
                Y = -(y - old_y) / 15;
            }
	*/
}

// GLUT callback
//    called on mouse movement
void mouse_motion (int x, int y) 
{
	//if(lbutton)
	{
		int rx,ry;

		ry = (GLfloat)(y - old_y)/2;	//Изменение угола поворота
		rx = (GLfloat)(x - old_x)/2;

		old_x=x;
		old_y=y;
		
		if(rx)
		{
			ort1 = Vector3D::RotateVector1AroundVector2(ort1,Vector3D(0,1,0),rx);
			ort2 = Vector3D::RotateVector1AroundVector2(ort2,Vector3D(0,1,0),rx);
			ort3 = Vector3D::RotateVector1AroundVector2(ort3,Vector3D(0,1,0),rx);
		}
		
		if(ry)
		{
			ort1 = Vector3D::RotateVector1AroundVector2(ort1,Vector3D(0,0,1),ry);
			ort2 = Vector3D::RotateVector1AroundVector2(ort2,Vector3D(0,0,1),ry);
			ort3 = Vector3D::RotateVector1AroundVector2(ort3,Vector3D(0,0,1),ry);
		}
	}
}



void Timer(int value)
{
	int c = clock();
	int work_time;
	sph_fluid_main_step();
	// Re-register for next callback
	//c = clock();
    glutTimerFunc(TIMER_INTERVAL*0, Timer, 0);
	//work_time = clock() - c;
	//printf("glutTimerFunc work:%d\n",work_time);
	//c = clock();
	glutPostRedisplay();
	//work_time = clock() - c;
	//printf("glutPostRedisplay work:%d\n",work_time);
}


 
int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 600);
    glutInitWindowPosition(100, 100);
	glutCreateWindow("Palyanov Andrey for OpenWorm: OpenCL SPH fluid + OpenGL(GLUT) visualization [2012]");
	
    /*
    glClearColor (0.3, 0.3, 0.3, 0.0); // цвет фона
    glEnable(GL_LIGHTING);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glEnable(GL_NORMALIZE);
	*/
	//===============================
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);

	glEnable(GL_NORMALIZE);
	glEnable(GL_AUTO_NORMAL);

	float ambient[4] = {1.0, 1.0, 1.0, 1};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

	//===============================

	engine = new Engine();
	engine->Init();

	sph_fluid_main_start();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMouseFunc(respond_mouse);
	glutMotionFunc(mouse_motion);	// The former handles movement while the mouse is clicked, 
	//glutPassiveMotionFunc			// and the latter while no button is clicked

	glutTimerFunc(TIMER_INTERVAL*0, Timer, 0);
	

	

	/*
	while(1)
	{
		sph_fluid_main_step();
	}
	/**/

	glutMainLoop();

	sph_fluid_main_stop();

	//exit(1);
	return 0;
}