/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/

#include <stdio.h>
#include <sstream>
#include <csignal>

#include "owWorldSimulation.h"

#include "owWorldSimulation.h"

extern int numOfLiquidP;
extern int numOfElasticP;
extern int numOfBoundaryP;
extern int numOfMembranes;
extern bool load_from_file;

int old_x=0, old_y=0;	// Used for mouse event
float camera_trans[] = {0, 0, -8.0};
float camera_rot[]   = {0, 0, 0};
float camera_trans_lag[] = {0, 0, -8.0};
float camera_rot_lag[] = {0, 0, 0};
const float inertia = 1.0f;
float modelView[16];
int buttonState = 0;

float sc = 0.045f;		//0.0145;//0.045;//0.07

float sc_scale = 1.0f;

Vector3D ort1(1,0,0),ort2(0,1,0),ort3(0,0,1);
GLsizei viewHeight, viewWidth;
int winIdMain;
int MUSCLE_COUNT = 100;//increase this value and modify corresponding code if you plan to add more than 10 muscles
double totalTime = 0;
int frames_counter = 0;
double calculationTime;
double renderTime;
double fps;
char device_full_name [1000];
double prevTime;
unsigned int * p_indexb;
float * d_cpp;
float * p_cpp;
float * ec_cpp;
float * muscle_activation_signal_cpp;
int   * md_cpp;// pointer to membraneData_cpp
owPhysicsFluidSimulator * fluid_simulation;
owHelper * helper;
owConfigProrerty * loacalConfig;
bool flag = false;
bool sPause = false;
void * m_font = (void *) GLUT_BITMAP_8_BY_13;
int iteration = 0;

void calculateFPS();
void drawScene();
void renderInfo(int,int);
void glPrint(float,float,const char *, void*);
void glPrint3D(float,float,float,const char *, void*);
//float muscle_activation_signal [10] = {0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f};
void beginWinCoords(void)
{
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

void endWinCoords(void)
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}
void glPrint(float x, float y, const char *s, void *font)
{
	glRasterPos2f((GLfloat)x, (GLfloat)y);
    int len = (int) strlen(s);
    for (int i = 0; i < len; i++) {
        glutBitmapCharacter(font, s[i]);
    }
}
void glPrint3D(float x, float y, float z, const char *s, void *font)
{
	glRasterPos3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
    int len = (int) strlen(s);
    for (int i = 0; i < len; i++) {
        glutBitmapCharacter(font, s[i]);
    }
}
/** Main displaying function
 */
void display(void)
{
	//Update Scene if not paused
	int i,j,k;
	int err_coord_cnt = 0;
	if(!sPause){
		if(load_from_file){
			owHelper::loadConfigurationFromFile_experemental(p_cpp,ec_cpp,md_cpp, loacalConfig,iteration);
			iteration++;
		}else{
			calculationTime = fluid_simulation->simulationStep(); // Run one simulation step
			int pib;
			p_indexb = fluid_simulation->getParticleIndex_cpp();
			for(i=0;i<loacalConfig->getParticleCount();i++)
			{
				pib = p_indexb[2*i + 1];
				p_indexb[2*pib + 0] = i;
			}
			p_cpp = fluid_simulation->getPosition_cpp();
			d_cpp = fluid_simulation->getDensity_cpp();
			ec_cpp = fluid_simulation->getElasticConnectionsData_cpp();
		}
		helper->refreshTime();
	}

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
	drawScene();
	glPointSize(3.f);
	glBegin(GL_POINTS);
	float dc, rho;
	//Display all particles
	for(i = 0; i<loacalConfig->getParticleCount(); i++)
	{
		if(!load_from_file){
			rho = d_cpp[ p_indexb[ i * 2 + 0 ] ];
			if( rho < 0 ) rho = 0;
			if( rho > 2 * rho0) rho = 2 * rho0;
			dc = 100.f * ( rho - rho0 ) / rho0 ;
			if(dc>1.f) dc = 1.f;
			//  R   G   B
			glColor4f(  0,  0,  1, 1.0f);//blue
			if(!load_from_file){
				if( (dc=100*(rho-rho0*1.00f)/rho0) >0 )	glColor4f(   0,  dc,   1,1.0f);//cyan
				if( (dc=100*(rho-rho0*1.01f)/rho0) >0 )	glColor4f(   0,   1,1-dc,1.0f);//green
				if( (dc=100*(rho-rho0*1.02f)/rho0) >0 )	glColor4f(  dc,   1,   0,1.0f);//yellow
				if( (dc=100*(rho-rho0*1.03f)/rho0) >0 )	glColor4f(   1,1-dc,   0,1.0f);//red
				if( (dc=100*(rho-rho0*1.04f)/rho0) >0 )	glColor4f(   1,   0,   0,1.0f);
			}
		}
		else
			glColor4f(  0,  0,  1, 1.0f);//blue
		if((int)p_cpp[i*4 + 3] != BOUNDARY_PARTICLE /*&& (int)p_cpp[i*4 + 3] != ELASTIC_PARTICLE*/)
		{
			glBegin(GL_POINTS);
			if((int)p_cpp[i*4+3]==2) 
			{ 
				glColor4f(   0,   0,   0,  1.0f);// color of elastic particles
				glPointSize(6.f);
			}
			glVertex3f( (p_cpp[i*4]-loacalConfig->xmax/2)*sc , (p_cpp[i*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[i*4+2]-loacalConfig->zmax/2)*sc );
			glPointSize(3.f);
			glEnd();

			if(!((p_cpp[i*4  ]>=0)&&(p_cpp[i*4  ]<=loacalConfig->xmax)&&
				(p_cpp[i*4+1]>=0)&&(p_cpp[i*4+1]<=loacalConfig->ymax)&&
				(p_cpp[i*4+2]>=0)&&(p_cpp[i*4+2]<=loacalConfig->zmax) ))
			{
				char label[50];
				beginWinCoords();
				glRasterPos2f (0.01F, 0.05F);
				if(err_coord_cnt<50){
				sprintf(label,"%d: %f , %f , %f",i,p_cpp[i*4  ],p_cpp[i*4+1],p_cpp[i*4+2]);
				glPrint( 0.f, (float)(50+err_coord_cnt*11), label, m_font);}
				if(err_coord_cnt==50) {
				glPrint( 0, (float)(50+err_coord_cnt*11), "............", m_font);}
				err_coord_cnt++;
				endWinCoords();
			}
		}
	}
	glLineWidth((GLfloat)0.1);
	int ecc=0;//elastic connections counter;
	//Display elastic connections
	for(int i_ec=0; i_ec < numOfElasticP * MAX_NEIGHBOR_COUNT; i_ec++)
	{
		//offset = 0
		if((j=(int)ec_cpp[ 4 * i_ec + 0 ])>=0)
		{
			i = (i_ec / MAX_NEIGHBOR_COUNT);// + (generateInitialConfiguration!=1)*numOfBoundaryP;
			if(i<j)	
			{
				glColor4b(150/2, 125/2, 0, 100/2/*alpha*/);
				if(ec_cpp[ 4 * i_ec + 2 ]>1.f)//muscles 
				{
					glLineWidth((GLfloat)1.0);
					if(ec_cpp[4*i_ec+2]-floor(ec_cpp[4*i_ec+2])>0.45f) 
					{
						if(muscle_activation_signal_cpp[ (int)(floor( ec_cpp[4*i_ec+2])-1) ]>0.1)
						glLineWidth((GLfloat)6.0); else glLineWidth((GLfloat)2.0);
						glColor4b(127/2, 0, 255/2, 255/2);/* muscle_number+0.5 <--> violet*/  
						glBegin(GL_LINES);
						glVertex3f( (p_cpp[i*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[i*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[i*4+2]-loacalConfig->zmax/2)*sc );
						glColor4b(255/2, 255/2, 255/2, 255/2);
						glVertex3f( (p_cpp[j*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[j*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[j*4+2]-loacalConfig->zmax/2)*sc );
						glEnd();
					}
					else
					if(ec_cpp[4*i_ec+2]-floor(ec_cpp[4*i_ec+2])>0.35f) 
					{ 
						if(muscle_activation_signal_cpp[ (int)(floor( ec_cpp[4*i_ec+2])-1) ]>0.1)
						glLineWidth((GLfloat)6.0); else glLineWidth((GLfloat)2.0);
						glColor4b(255/2, 0, 255/2, 255/2);/* muscle_number+0.4 <--> magenta*/  
						glBegin(GL_LINES);
						glVertex3f( (p_cpp[i*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[i*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[i*4+2]-loacalConfig->zmax/2)*sc );
						glColor4b(255/2, 255/2, 255/2, 255/2);
						glVertex3f( (p_cpp[j*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[j*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[j*4+2]-loacalConfig->zmax/2)*sc );
						glEnd();
					}
					else
					if(ec_cpp[4*i_ec+2]-floor(ec_cpp[4*i_ec+2])>0.25f) 
					{ 
						if(muscle_activation_signal_cpp[ (int)(floor( ec_cpp[4*i_ec+2])-1) ]>0.1)
						glLineWidth((GLfloat)6.0); else glLineWidth((GLfloat)2.0);
						glColor4b(255/2, 127/2, 0, 255/2);/* muscle_number+0.3 <--> orange*/  
						glBegin(GL_LINES);
						glVertex3f( (p_cpp[i*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[i*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[i*4+2]-loacalConfig->zmax/2)*sc );
						glColor4b(255/2, 255/2, 255/2, 255/2);
						glVertex3f( (p_cpp[j*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[j*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[j*4+2]-loacalConfig->zmax/2)*sc );
						glEnd();
					}
					else
					if(ec_cpp[4*i_ec+2]-floor(ec_cpp[4*i_ec+2])>0.15f) 
					{ 
						if(muscle_activation_signal_cpp[ (int)(floor( ec_cpp[4*i_ec+2])-1) ]>0.1)
						glLineWidth((GLfloat)6.0); else glLineWidth((GLfloat)2.0);
						glColor4b(255/2, 0, 0, 255/2);/* muscle_number+0.2 <--> red*/  
						glBegin(GL_LINES);
						glVertex3f( (p_cpp[i*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[i*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[i*4+2]-loacalConfig->zmax/2)*sc );
						glColor4b(255/2, 255/2, 255/2, 255/2);
						glVertex3f( (p_cpp[j*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[j*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[j*4+2]-loacalConfig->zmax/2)*sc );
						glEnd();
					}
					else
					{
						glColor4b(255/2, 0,     0, 255/2);/* muscle_number+0.1 <--> red */

						glBegin(GL_LINES);
						glVertex3f( (p_cpp[i*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[i*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[i*4+2]-loacalConfig->zmax/2)*sc );
						glVertex3f( (p_cpp[j*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[j*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[j*4+2]-loacalConfig->zmax/2)*sc );
						glEnd();
					}
				}
				else
				{//ordinary springs
					glLineWidth((GLfloat)0.1);
					glBegin(GL_LINES);
											glColor4b(150/2, 125/2, 0, 100/2);
					if(p_cpp[i*4+3]>2.15)	glColor4b( 50/2, 125/2, 0, 100/2);
					glVertex3f( (p_cpp[i*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[i*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[i*4+2]-loacalConfig->zmax/2)*sc );
											glColor4b(150/2, 125/2, 0, 100/2);
					if(p_cpp[j*4+3]>2.15)	glColor4b( 50/2, 125/2, 0, 100/2);
					glVertex3f( (p_cpp[j*4+0]-loacalConfig->xmax/2)*sc , (p_cpp[j*4+1]-loacalConfig->ymax/2)*sc, (p_cpp[j*4+2]-loacalConfig->zmax/2)*sc );
					glEnd();
				}
				
				ecc++;
			}
		}
	}
	// Draw membranes
	if(!load_from_file)
		md_cpp = fluid_simulation->getMembraneData_cpp();
	glColor4b(0, 200/2, 150/2, 255/2/*alpha*/);
	for(int i_m = 0; i_m < numOfMembranes; i_m++)
	{
		i = md_cpp [i_m*3+0];
		j = md_cpp [i_m*3+1];
		k = md_cpp [i_m*3+2];

		glBegin(GL_LINES);
		glVertex3f( ((p_cpp[i*4]+p_cpp[j*4]+4*p_cpp[k*4])/6-loacalConfig->xmax/2)*sc , ((p_cpp[i*4+1]+p_cpp[j*4+1]+4*p_cpp[k*4+1])/6-loacalConfig->ymax/2)*sc, ((p_cpp[i*4+2]+p_cpp[j*4+2]+4*p_cpp[k*4+2])/6-loacalConfig->zmax/2)*sc );
		glVertex3f( ((p_cpp[i*4]+p_cpp[k*4]+4*p_cpp[j*4])/6-loacalConfig->xmax/2)*sc , ((p_cpp[i*4+1]+p_cpp[k*4+1]+4*p_cpp[j*4+1])/6-loacalConfig->ymax/2)*sc, ((p_cpp[i*4+2]+p_cpp[k*4+2]+4*p_cpp[j*4+2])/6-loacalConfig->zmax/2)*sc );

		glVertex3f( ((p_cpp[i*4]+p_cpp[k*4]+4*p_cpp[j*4])/6-loacalConfig->xmax/2)*sc , ((p_cpp[i*4+1]+p_cpp[k*4+1]+4*p_cpp[j*4+1])/6-loacalConfig->ymax/2)*sc, ((p_cpp[i*4+2]+p_cpp[k*4+2]+4*p_cpp[j*4+2])/6-loacalConfig->zmax/2)*sc );
		glVertex3f( ((p_cpp[j*4]+p_cpp[k*4]+4*p_cpp[i*4])/6-loacalConfig->xmax/2)*sc , ((p_cpp[j*4+1]+p_cpp[k*4+1]+4*p_cpp[i*4+1])/6-loacalConfig->ymax/2)*sc, ((p_cpp[j*4+2]+p_cpp[k*4+2]+4*p_cpp[i*4+2])/6-loacalConfig->zmax/2)*sc );

		glVertex3f( ((p_cpp[j*4]+p_cpp[k*4]+4*p_cpp[i*4])/6-loacalConfig->xmax/2)*sc , ((p_cpp[j*4+1]+p_cpp[k*4+1]+4*p_cpp[i*4+1])/6-loacalConfig->ymax/2)*sc, ((p_cpp[j*4+2]+p_cpp[k*4+2]+4*p_cpp[i*4+2])/6-loacalConfig->zmax/2)*sc );
		glVertex3f( ((p_cpp[i*4]+p_cpp[j*4]+4*p_cpp[k*4])/6-loacalConfig->xmax/2)*sc , ((p_cpp[i*4+1]+p_cpp[j*4+1]+4*p_cpp[k*4+1])/6-loacalConfig->ymax/2)*sc, ((p_cpp[i*4+2]+p_cpp[j*4+2]+4*p_cpp[k*4+2])/6-loacalConfig->zmax/2)*sc );
		glEnd();
	}
	glLineWidth((GLfloat)1.0);
	glutSwapBuffers();
	helper->watch_report("graphics: \t\t%9.3f ms\n====================================\n");
	renderTime = helper->get_elapsedTime();
	totalTime += calculationTime + renderTime;
	calculateFPS();
}
/** Drawing main scene and bounding box
 */
inline void drawScene()
{
	//       [7]----[6]
	//      / |     /| 
	//    [3]----[2] | 
	//     | [4]--|-[5]   
	//     | /    | /
	//    [0]----[1]  
	Vector3D vcenter(0,0,0);
	Vector3D vbox[8];
	float s_v = 1 /(simulationScale);// = 1 m in simulation
	float order = 0;
	while(s_v >= 1){
		s_v /= 10;
		if(s_v < 1)
		{
			s_v *= 10;
			break;
		}
		++order;
	}
	vbox[0] = Vector3D(loacalConfig->xmin,loacalConfig->ymin,loacalConfig->zmin);
	vbox[1] = Vector3D(loacalConfig->xmax,loacalConfig->ymin,loacalConfig->zmin);
	vbox[2] = Vector3D(loacalConfig->xmax,loacalConfig->ymax,loacalConfig->zmin);
	vbox[3] = Vector3D(loacalConfig->xmin,loacalConfig->ymax,loacalConfig->zmin);
	vbox[4] = Vector3D(loacalConfig->xmin,loacalConfig->ymin,loacalConfig->zmax);
	vbox[5] = Vector3D(loacalConfig->xmax,loacalConfig->ymin,loacalConfig->zmax);
	vbox[6] = Vector3D(loacalConfig->xmax,loacalConfig->ymax,loacalConfig->zmax);
	vbox[7] = Vector3D(loacalConfig->xmin,loacalConfig->ymax,loacalConfig->zmax);
	// Display user interface if enabled
	bool displayInfos = true;
    if (displayInfos) 
    {
        glDisable(GL_DEPTH_TEST);
        glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ZERO); // invert color
        glEnable(GL_BLEND);
        renderInfo(0, 0);
        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
    }
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
	vcenter = Vector3D(-(loacalConfig->xmin+loacalConfig->xmax)/2,-(loacalConfig->ymin+loacalConfig->ymax)/2,-(loacalConfig->zmin+loacalConfig->zmax)/2);
	vcenter *= sc;
	Vector3D v1,v2,v3,v4,v5,v6,v7,v8;
	v1 = Vector3D( -loacalConfig->xmax/2, -loacalConfig->ymax/2, -loacalConfig->zmax/2)*sc;
	v2 = Vector3D(  loacalConfig->xmax/2, -loacalConfig->ymax/2, -loacalConfig->zmax/2)*sc;
	v3 = Vector3D(  loacalConfig->xmax/2,  loacalConfig->ymax/2, -loacalConfig->zmax/2)*sc;
	v4 = Vector3D( -loacalConfig->xmax/2,  loacalConfig->ymax/2, -loacalConfig->zmax/2)*sc;
	v5 = Vector3D( -loacalConfig->xmax/2, -loacalConfig->ymax/2,  loacalConfig->zmax/2)*sc;
	v6 = Vector3D(  loacalConfig->xmax/2, -loacalConfig->ymax/2,  loacalConfig->zmax/2)*sc;
	v7 = Vector3D(  loacalConfig->xmax/2,  loacalConfig->ymax/2,  loacalConfig->zmax/2)*sc;
	v8 = Vector3D( -loacalConfig->xmax/2,  loacalConfig->ymax/2,  loacalConfig->zmax/2)*sc;
	glColor3ub(255,255,255);//yellow
	glVertex3d(v1.x,v1.y,v1.z); 
	glVertex3d(v2.x,v2.y,v2.z);
	

	glColor3ub(255,255,255);//yellow
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
	glVertex3d(v8.x + s_v*sc,v8.y,v8.z);

	glVertex3d(v8.x,v8.y,v8.z);
	glVertex3d(v5.x,v5.y,v5.z);
	glEnd();
	//
	glBegin(GL_LINES);
	glColor3ub(0,0,0);//black
	

	Vector3D v_s = Vector3D(  -loacalConfig->xmax/2 + s_v,  loacalConfig->ymax/2,  loacalConfig->zmax/2)*sc;
	glVertex3d(v_s.x, v_s.y, v_s.z);
	glVertex3d(v_s.x, v_s.y - 0.5f * sc , v_s.z);
	glLineWidth((GLfloat)10.0);
	glVertex3d( v8.x,  v8.y,  v8.z);
	glVertex3d(v_s.x, v_s.y, v_s.z);

	glEnd();
	glLineWidth((GLfloat)1.0);
	std::stringstream ss;
	std::string s;
	ss << order;
	s = "1E-" + ss.str() + "m";
	glPrint3D( (float)v8.x + 0.4f*sc , (float)v8.y - 2.f * sc, (float)v8.z, "0", m_font);
	glPrint3D( (float)v_s.x , (float)v_s.y - 2.f * sc, (float)v_s.z, s.c_str(), m_font);
	ss.str("");
	while(v_s.x < loacalConfig->xmax/2*sc){
		v_s.x += s_v * sc;
		if(v_s.x < loacalConfig->xmax/2*sc){
			glBegin(GL_LINES);
				glVertex3d(v_s.x, v_s.y, v_s.z);
				glVertex3d(v_s.x, v_s.y - 0.5f * sc , v_s.z);
			glEnd();
		}
	}
}

int count_s = 0;
float current_sv ;
static char label[1000];                            /* Storage for current string   */
bool showInfo = true;
bool showRuler = false;
/** Render addition test information
 */
void renderInfo(int x, int y)
{
	beginWinCoords();
	int y_m = y;
	if(showInfo){
		glColor3f (0.5F, 1.0F, 1.0F);
		sprintf(label,"Liquid particles: %d, elastic matter particles: %d, boundary particles: %d; total count: %d", numOfLiquidP,
																													 numOfElasticP,
																													 numOfBoundaryP,loacalConfig->getParticleCount());
		glPrint( 0 , 2 , label, m_font);
		glColor3f (1.0F, 1.0F, 1.0F); 
		sprintf(label,"Selected device: %s FPS = %.2f, time step: %d (%f s)", device_full_name+7, fps, fluid_simulation->getIteration(),((float)fluid_simulation->getIteration())*timeStep);
		glPrint( 0 , 17 , label, m_font);
	}
	if(showRuler){
		glColor3ub(255, 0, 0);
		float accuracy = 100;//what it it?
		float s_v = 1 * sc_scale * (1 /( accuracy * simulationScale));
		float s_v_10 = s_v / 10;
		std::stringstream ss;
		std::string s;
		glBegin(GL_LINES);
			glColor3f(1.0f,0.0f,0.0f);
			glVertex2f((GLfloat) 0.f,(GLfloat)y_m );
			glVertex2f((GLfloat) s_v,(GLfloat)y_m );
			glVertex2f((GLfloat) s_v,(GLfloat)y_m );
			glVertex2f((GLfloat) s_v,(GLfloat)y_m + 5.f );
		glEnd();
			glPrint( s_v , y_m + 15.f , "1E-02 m", m_font);
		glBegin(GL_LINES);		
			glVertex2f((GLfloat) s_v_10,(GLfloat)y_m + 0.f);
			glVertex2f((GLfloat) s_v_10,(GLfloat)y_m + 5.f);
		glEnd();
		if( 8 * s_v/pow(10.f,count_s) >= glutGet(GLUT_WINDOW_WIDTH)/2 ){
			count_s++;
			flag = true;
		}else{
			if(count_s != 0 && 8 * s_v/pow(10.f,count_s - 1) < glutGet(GLUT_WINDOW_WIDTH)/2){
				//flag = false;
				count_s--;
			}
		}
		if(flag){
			for(int i = 1;i <= count_s; i++){
				glBegin(GL_LINES);		
					glVertex2f((GLfloat) s_v/pow(10.f,i + 1),(GLfloat)y_m + 0.f);
					glVertex2f((GLfloat) s_v/pow(10.f,i + 1),(GLfloat)y_m + 5.f);
				glEnd();
				ss << i + 1 + 2;
				s = "1E-" + ss.str() + "m";
				ss.str("");
				glPrint( s_v/pow(10.f,i + 1) , y_m + 15.f , s.c_str(), m_font);
			}
		}
	}
	endWinCoords();
}
/** Calculation of FPS
 */
void calculateFPS()
{
    //  Increase frame count
	frames_counter++;
    int timeInterval = (int)(totalTime - prevTime);
    if(timeInterval >= 1000)
    {
		fps = frames_counter / (timeInterval / 1000.0f);
        prevTime = totalTime;
        frames_counter = 0;
		printf("FPS: \t\t%9.3f fps\n====================================\n",	fps );
    }
}
void respond_mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON)
		buttonState = 1;
	if (button == GLUT_RIGHT_BUTTON)
		buttonState = 3;
	int mods;
	mods = glutGetModifiers();
    if (mods & GLUT_ACTIVE_CTRL) 
    {
        buttonState = 2;
    } 
	if(state == GLUT_UP)
		buttonState = 0;
	old_x=x;
	old_y=y;
	if (button == 3)     // mouse wheel up
    {
        sc *= 1.1f;		 // Zoom in
		sc_scale *= 1.1f;// Zoom in
    }
    else
	if (button == 4)	 // mouse wheel down
    {
        sc /= 1.1f;		 // Zoom out
		sc_scale /= 1.1f;// Zoom out
    }
}

// GLUT callback
// called on mouse movement

void mouse_motion (int x, int y) 
{
	float dx,dy;
	dy = (float)(y - old_y);	
	dx = (float)(x - old_x);
	
	if(buttonState == 1)
	{
		camera_rot[0] += dy / 5.0f;
		camera_rot[1] += dx / 5.0f;
	}
	if(buttonState == 3){
		// middle = translate
		camera_trans[0] += dx / 100.0f;
		camera_trans[1] -= dy / 100.0f;
		//camera_trans[2] += (dy / 100.0f) * 0.5f * fabs(camera_trans[2]);
	}
	if(buttonState == 2){
		// middle = translate
		//camera_trans[0] += dx / 100.0f;
		//camera_trans[1] -= dy / 100.0f;
		camera_trans[2] += (dy / 100.0f) * 0.5f * fabs(camera_trans[2]);
	}
	old_x=x;
	old_y=y;
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	for (int c = 0; c < 3; ++c)
	{
	    camera_trans_lag[c] += (camera_trans[c] - camera_trans_lag[c]) * inertia;
		camera_rot_lag[c] += (camera_rot[c] - camera_rot_lag[c]) * inertia;
	}
    glTranslatef(camera_trans_lag[0], camera_trans_lag[1], camera_trans_lag[2]);
	glRotatef(camera_rot_lag[0], 1.0, 0.0, 0.0);
	glRotatef(camera_rot_lag[1], 0.0, 1.0, 0.0);
	glGetFloatv(GL_MODELVIEW_MATRIX, modelView);
}

extern float *muscle_activation_signal_cpp;
void Cleanup(int exitCode){
	delete fluid_simulation;
	delete helper;
	exit (exitCode);
}
void RespondKey(unsigned char key, int x, int y)
{
	switch(key)
	{
	case '1':
		owHelper::suffix = "";
		helper->refreshTime();
		fluid_simulation->reset();
		sPause = false;
		break;
	case '2':
		owHelper::suffix = "_membranes_demo";
		helper->refreshTime();
		fluid_simulation->reset();
		sPause = false;
		break;
	case '\033':// Escape quits
	case 'Q':   // Q quits
	case 'q':   // q quits
		Cleanup(EXIT_SUCCESS);
		break;
	case ' ':
		sPause = !sPause;
		std::cout << "\nSimulation Is Paused" << std::endl;
		break;
	}

	if(key == 'i')
	{
		showInfo = !showInfo;
	}
	if(key == 'r')
	{
		showRuler = !showRuler;
	}
	glutPostRedisplay();
}

//Auxiliary function
/** There can be only one idle() callback function. In an
 *  animation, this idle() function must update not only the
 *  main window but also all derived subwindows
 */
void idle (void) 
{ 
  glutSetWindow (winIdMain); 
  glutPostRedisplay (); 
} 
//static char label[1000];                            /* Storage for current string   */

void Timer(int value)
{
	// Re-register for next callback
    glutPostRedisplay();
    glutTimerFunc(TIMER_INTERVAL*0, Timer, 0);
}

GLvoid resize(GLsizei width, GLsizei height){

	if(height == 0) { height = 1; }										 
	if(width == 0) { width = 1; }										 
	
	glViewport(0, 0, width, height);					// Set view area
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	float aspectRatio = (GLfloat)width / (GLfloat)height;
	if (aspectRatio>1.f)
		glFrustum(-1*aspectRatio, 1*aspectRatio, -1, 1, 3, 45);
	else
		glFrustum(-1, 1, -1/aspectRatio, 1/aspectRatio, 3, 45);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	for (int c = 0; c < 3; ++c)
	{
	    camera_trans_lag[c] += (camera_trans[c] - camera_trans_lag[c]) * inertia;
		camera_rot_lag[c] += (camera_rot[c] - camera_rot_lag[c]) * inertia;
	}
    glTranslatef(camera_trans_lag[0], camera_trans_lag[1], camera_trans_lag[2]);
	glRotatef(camera_rot_lag[0], 1.0, 0.0, 0.0);
	glRotatef(camera_rot_lag[1], 0.0, 1.0, 0.0);
	glGetFloatv(GL_MODELVIEW_MATRIX, modelView);
}
inline void init(void){
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
	glEnable(GL_AUTO_NORMAL);
	float ambient[4] = {1.0, 1.0, 1.0, 1};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
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
void sighandler(int s){
	std::cerr << "\nCaught signal CTRL+C. Exit Simulation..." << "\n"; // this is undefined behaviour should check signal value
	Cleanup(EXIT_SUCCESS);
}
/** Init & start simulation and graphic component if with_graphics==true
 *
 * 	@param argc
 * 	command line arguments going throw the main function
 * 	@param with_graphics
 * 	Flag indicates that simulation will be run with graphic or not
 * 	@param load_to
 * 	Flag indicates that simulation will in "load configuration to file" mode
 */
void run(int argc, char** argv, const bool with_graphics, const bool load_to)
{
	helper = new owHelper();
	if(!load_from_file){
		DEVICE dev_type = CPU;
		for(int i = 1; i<argc; i++){
			if(strncmp(argv[i], "device=", 7) == 0){
				if(strstr(argv[i], "gpu") != NULL || strstr(argv[i], "GPU") != NULL)
					dev_type = GPU;
			}
		}
		fluid_simulation = new owPhysicsFluidSimulator(helper, dev_type);
		loacalConfig = fluid_simulation->getConfig();
	}
	else{
		loacalConfig = new owConfigProrerty();
		muscle_activation_signal_cpp = new float [MUSCLE_COUNT];
		for(int i=0;i<MUSCLE_COUNT;i++)
		{
			muscle_activation_signal_cpp[i] = 0.f;
		}
	}
	std::signal(SIGINT,sighandler);
	if(with_graphics){
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
		glutInitWindowSize(1024, 1024);
		glutInitWindowPosition(100, 100);
		winIdMain = glutCreateWindow("Palyanov Andrey for OpenWorm: OpenCL PCISPH fluid + elastic matter + membranes [2013]: C.elegans body generator demo");
		glutIdleFunc (idle); 
		//Init physic Simulation
		init();
		glutDisplayFunc(display);
		glutReshapeFunc(resize);
		glutMouseFunc(respond_mouse);
		glutMotionFunc(mouse_motion);	//process movement in case if the mouse is clicked,
		glutKeyboardFunc(RespondKey);
		glutTimerFunc(TIMER_INTERVAL * 0, Timer, 0);
		glutMainLoop();
		if(!load_from_file){
			Cleanup(EXIT_SUCCESS);
		}
	}else{
		while(1){
			fluid_simulation->simulationStep(load_to);
			helper->refreshTime();
		}
	}
    exit(EXIT_SUCCESS);
}
