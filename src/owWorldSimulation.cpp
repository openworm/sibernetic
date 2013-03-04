#include "owWorldSimulation.h"
#include <stdio.h>

extern int numOfLiquidP;
extern int numOfElasticP;
extern int numOfBoundaryP;
extern int iterationCount;

bool rotate = false;
int old_x=0, old_y=0;	// Used for mouse event
float rotX = 0.0f;		// Rotate screen on x axis 
float rotY = 0.0f;		// Rotate screen on y axis
float rotZ = 0.0f;		// Rotate screen on z axis
bool lbutton = false;
float sc = 0.025;		//0.0145;//0.045;//0.07
Vector3D ort1(1,0,0),ort2(0,1,0),ort3(0,0,1);
GLsizei Height, Width;
int winIdMain;
int winIdSub;
owPhysicsFluidSimulator * fluid_simulation;
owHelper * helper;
double calculationTime;
double renderTime;
double totalTime = 0;
int frames_counter = 0;
double fps;
char device_full_name [1000];
double prevTime;
int frameCount = 0;
void calculateFPS();
unsigned int * p_indexb;
float * d_b;
float * p_b;
void display(void)
{
	helper->refreshTime();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
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
	glColor3ub(255,255,255);//yellow
	glVertex3d(v1.x,v1.y,v1.z); 
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

	//glColor3ub(255,255,255);//yellow
	p_indexb = fluid_simulation->getParticleIndexBuffer();
	int pib;
	for(int i=0;i<PARTICLE_COUNT;i++)
	{
		pib = p_indexb[2*i + 1];
		p_indexb[2*pib + 0] = i;
	}
	glPointSize(3.f);
	glBegin(GL_POINTS);
	p_b = fluid_simulation->getPositionBuffer();
	d_b = fluid_simulation->getDensityBuffer();
	float dc, rho;
	for(int i = 0; i<PARTICLE_COUNT; i++)
	{
		rho = d_b[ p_indexb[ i * 2 + 0 ] ];
		if( rho < 0 ) rho = 0;
		if( rho > 2 * rho0) rho = 2 * rho0;
		dc = 100.0 * ( rho - rho0 ) / rho0 ;
		if(dc>1.f) dc = 1.f;
		//  R   G   B
		glColor4f(  0,  0,  1, 1.0f);//blue
		if( (dc=100*(rho-rho0*1.00f)/rho0) >0 )	glColor4f(   0,  dc,   1,1.0f);//cyan
		if( (dc=100*(rho-rho0*1.01f)/rho0) >0 )	glColor4f(   0,   1,1-dc,1.0f);//green
		if( (dc=100*(rho-rho0*1.02f)/rho0) >0 )	glColor4f(  dc,   1,   0,1.0f);//yellow
		if( (dc=100*(rho-rho0*1.03f)/rho0) >0 )	glColor4f(   1,1-dc,   0,1.0f);//red
		if( (dc=100*(rho-rho0*1.04f)/rho0) >0 )	glColor4f(   1,   0,   0,1.0f);
		if((int)p_b[i*4 + 3] != BOUNDARY_PARTICLE /*&& (int)p_b[i*4 + 3] != ELASTIC_PARTICLE*/){
			glBegin(GL_POINTS);
			if((int)p_b[i*4+3]==2) glColor4f(   1,   1,   0,  1.0f);
			glVertex3f( (p_b[i*4]-XMAX/2)*sc , (p_b[i*4+1]-YMAX/2)*sc, (p_b[i*4+2]-ZMAX/2)*sc );
			glEnd();
		}
	}
	glEnd();
	glutSwapBuffers();
	helper->watch_report("graphics: \t\t%9.3f ms\n====================================\n");
	renderTime = helper->get_elapsedTime();
	totalTime += calculationTime + renderTime;
	calculateFPS();
}
void calculateFPS()
{
    //  Increase frame count
	frames_counter++;
    int timeInterval = totalTime - prevTime;
    if(timeInterval > 1000)
    {
		fps = frames_counter / (totalTime / 1000.0f);
        prevTime = totalTime;
        frameCount = 0;
		printf("FPS: \t\t%9.3f fps\n====================================\n",	fps );
    }
}/**/
void respond_mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON)
		lbutton = true;
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
}

// GLUT callback
// called on mouse movement
void mouse_motion (int x, int y) 
{
	if(lbutton)
	{
		int rx,ry;
		ry = (GLfloat)(y - old_y)/2;	
		rx = (GLfloat)(x - old_x)/2;

		old_x=x;
		old_y=y;
		
		if(rx)
		{
			rotX += rx;
			glRotatef(rx, 0.0, 1.0, 0.0);          
		}
		if(ry)
		{
		}
	}
}
//Auxiliary function
/* There can be only one idle() callback function. In an 
   animation, this idle() function must update not only the 
   main window but also all derived subwindows */ 
void idle (void) 
{ 
  glutSetWindow (winIdMain); 
  glutPostRedisplay (); 
  glutSetWindow (winIdSub); 
  glutPostRedisplay (); 
}; 
void drawString (char *s) 
{ 
  unsigned int i; 
  for (i = 0; i < strlen (s); i++) 
    glutBitmapCharacter (GLUT_BITMAP_HELVETICA_10, s[i]); 
}; 
void drawStringBig (char *s) 
{ 
  unsigned int i; 
  for (i = 0; i < strlen (s); i++) 
	  glutBitmapCharacter (GLUT_BITMAP_HELVETICA_18, s[i]); 
}; 
static char label[100];                            /* Storage for current string   */

void subMenuDisplay() 
{ 
	/* Clear subwindow */ 
	glutSetWindow (winIdSub); 
	glClearColor(0.7f, 0.7f, 0.7f, 1.0f);
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	/* Write State Variables */ 
	glColor3f (1.0F, 1.0F, 1.0F);
	sprintf(label,"Liquid particles: %d, elastic matter particles: %d, boundary particles: %d; total count: %d", numOfLiquidP,
																												 numOfElasticP,
																												 numOfBoundaryP,PARTICLE_COUNT); 
	glRasterPos2f (0.01F, 0.65F); 
	drawStringBig (label); 
	glColor3f (1.0F, 1.0F, 1.0F); 
	sprintf(label,"Selected device: %s     FPS = %.2f, time step: %d", device_full_name, fps, iterationCount); 
	glRasterPos2f (0.01F, 0.20F); 
	drawStringBig (label); 

	glutSwapBuffers (); 
}; 
void subMenuReshape (int w, int h) 
{ 
  glViewport (0, 0, w, h); 
  glMatrixMode (GL_PROJECTION); 
  glLoadIdentity (); 
  gluOrtho2D (0.0F, 1.0F, 0.0F, 1.0F); 
}; 
void Timer(int value)
{
	calculationTime = fluid_simulation->simulationStep();
	// Re-register for next callback
    glutTimerFunc(TIMER_INTERVAL*0, Timer, 0);
	glutPostRedisplay();
}
void SetProjectionMatrix(void){
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();									// Current projection matrix is dropped to identity matrix 
	glFrustum(-1, 1, -1, 1, 3, 15);						// Set up perspective projection
}
void SetModelviewMatrix(void){
     glMatrixMode(GL_MODELVIEW);                                   
     glLoadIdentity();                                             
     glTranslatef(0.0, 0.0, -8.0);                              
     glRotatef(10.0, 1.0, 0.0, 0.0);
     glRotatef(rotX, 0.0, 1.0, 0.0);                              
}
GLvoid resize(GLsizei width, GLsizei height){
	if(height == 0)
	{
		height = 1;										
	}

	glViewport(0, 0, width, height);					// Sev view area

	Height = height;
	Width = width;

	SetProjectionMatrix();
	SetModelviewMatrix();
}
void init(void){
	glClearColor(0.7f, 0.7f, 0.7f, 1.0f);
	glClearDepth(1.0f);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

}
void draw(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	glColor3f(1.0f, 1.0f, 1.0f);
	glPopMatrix();
}

void run(int argc, char** argv, const bool with_graphics)
{
	helper = new owHelper();
	fluid_simulation = new owPhysicsFluidSimulator(helper);
	if(with_graphics){
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
		glutInitWindowSize(800, 600);
		glutInitWindowPosition(100, 100);
		winIdMain = glutCreateWindow("Palyanov Andrey for OpenWorm: OpenCL PCISPH fluid + elastic matter demo [2013]");
		glutIdleFunc (idle); 
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_COLOR_MATERIAL);
		glEnable(GL_NORMALIZE);
		glEnable(GL_AUTO_NORMAL);
		float ambient[4] = {1.0, 1.0, 1.0, 1};
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
		//Init physic Simulation
		init();
		glutDisplayFunc(display);
		glutReshapeFunc(resize);
		glutMouseFunc(respond_mouse);
		glutMotionFunc(mouse_motion);	// The former handles movement while the mouse is clicked, 
		//Create sub window which contains information about simulation: FPS, and particles count
		winIdSub = glutCreateSubWindow (winIdMain, 5, 5, 1000 - 10, 600 / 10); 
		glutDisplayFunc (subMenuDisplay); 
		glutReshapeFunc (subMenuReshape); 
		glutTimerFunc(TIMER_INTERVAL*0, Timer, 0);
		glutMainLoop();
		fluid_simulation->~owPhysicsFluidSimulator();
	}else{
		while(1){
			fluid_simulation->simulationStep();
			helper->refreshTime();
		}
	}
}
