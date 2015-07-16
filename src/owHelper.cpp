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

#include <string>
#include <stdio.h>
#include <sstream>
#include <string>
#include <vector>


#if defined(__APPLE__) || defined (__MACOSX)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

#define EXPEREMENTAL_WRITE 0
#if EXPEREMENTAL_WRITE
#define EXPEREMENTAL_READ
#endif
#include "owHelper.h"
#include "owPhysicsConstant.h"

using namespace std;

extern int numOfMembranes;
extern int numOfElasticP;
extern int numOfLiquidP;


/** owHelpre class constructor
 */
owHelper::owHelper(void)
{
	refreshTime();
}
/** owHelpre class destructor
 */
owHelper::~owHelper(void)
{
}
/** Update time variable
 *
 *  Using precision system functions for getting exact system time.
 *  Initialization of start time value.
 *  For more info:
 *  Windows QueryPerformanceFrequency - http://msdn.microsoft.com/ru-ru/library/windows/desktop/ms644905(v=vs.85).aspx
 *  Linux clock_gettime - http://man7.org/linux/man-pages/man2/clock_gettime.2.html
 *  MacOS - https://developer.apple.com/library/mac/documentation/Darwin/Conceptual/KernelProgramming/services/services.html
 */
void owHelper::refreshTime()
{
#if defined(_WIN32) || defined (_WIN64)
    QueryPerformanceFrequency(&frequency);	// get ticks per second
    QueryPerformanceCounter(&t1);			// start timer
	t0 = t1;
#elif defined(__linux__)
	clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
	t0 = t1;
#elif defined(__APPLE__)
    t1 = mach_absolute_time();
    t0 = t1;
#endif
}

int generateWormShell(int stage, int i_start,float *position_cpp, float *velocity_cpp, int &numOfMembranes, int *membraneData_cpp, owConfigProrerty * config)
{
	//return 0;
	if( !((stage==0)||(stage==1)) ) return 0;

	float alpha;
	float wormBodyRadius;
	int pCount = 0;//total counter of particles being created within this function
	int i,j;
	float *positionVector;
	float *velocityVector;
	float xc = config->xmax*0.5f;
	float yc = config->ymax*0.3f;
	float zc = config->zmax*0.5f;
	int elasticLayers = 1;//starting value
	float PI = 3.1415926536f;
	int currSlice_pCount;
	int prevSlice_pCount = 0;
	int currSlice_start;
	int prevSlice_start = 0;
	int mc = 0;
	float angle;
	int tip = 0;

	int jmin=-100,jmax=98;

	if(stage==0) numOfMembranes=0;

	//return 0;

	for(j=jmin;j<=jmax;j++)	// number of cross-slices of the worm, in direction from head to tail. -
							// 1..-1 for example will give 3 slices in the middle of the worm
	{////////////////////////////////////////////////////
		currSlice_pCount = 0;
		currSlice_start = pCount;
		wormBodyRadius = 6.0f*r0*sqrt(max(1.f-(1.0e-4f)*j*j,0.f));
		tip = 0;

		if((wormBodyRadius>0.707f*r0)&&
		   (wormBodyRadius<1.000f*r0)) wormBodyRadius = 1.000f*r0;

		if(wormBodyRadius<0.707f*r0) { tip = 1; wormBodyRadius = 0.707f*r0; }//0.707 = sqrt(2)/2

		//alpha = 2*asin(0.5*r0/wormBodyRadius);//in radians
		//angle = alpha;

		elasticLayers = 1;

		if(stage==1)
		{
			positionVector = position_cpp + 4 * (pCount+i_start);
			positionVector[ 0 ] = xc + wormBodyRadius*cos(0.0f);
			positionVector[ 1 ] = yc + wormBodyRadius*sin(0.0f);
			positionVector[ 2 ] = zc + r0*j;
			positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow

			positionVector = position_cpp + 4 * (pCount+1+i_start);
			positionVector[ 0 ] = xc - wormBodyRadius*cos(0.0f);
			positionVector[ 1 ] = yc - wormBodyRadius*sin(0.0f);
			positionVector[ 2 ] = zc + r0*j;
			positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow
		}
		pCount+=2;
		currSlice_pCount+=2;

		if(tip==1)
		{
			if(stage==1)
			{
				positionVector = position_cpp + 4 * (pCount+i_start);
				positionVector[ 0 ] = xc + wormBodyRadius*sin(0.0f);
				positionVector[ 1 ] = yc + wormBodyRadius*cos(0.0f);
				positionVector[ 2 ] = zc + r0*j;
				positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow

				positionVector = position_cpp + 4 * (pCount+1+i_start);
				positionVector[ 0 ] = xc - wormBodyRadius*sin(0.0f);
				positionVector[ 1 ] = yc - wormBodyRadius*cos(0.0f);
				positionVector[ 2 ] = zc + r0*j;
				positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow
			}
			pCount+=2;
			currSlice_pCount+=2;

		}

		while(elasticLayers<=2)//this number defines number of radial muscle layers. can be 1, 2, 3 or more
		{
			if((elasticLayers==2)&&(j==jmin))
			{
				if(stage==1)
				{
					positionVector = position_cpp + 4 * (pCount+i_start);
					positionVector[ 0 ] = xc;
					positionVector[ 1 ] = yc;
					positionVector[ 2 ] = zc + r0*(j-1);
					positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow
				}
				pCount++;
				currSlice_pCount++;
			}

			if(wormBodyRadius>0)
			if(elasticLayers>=2)
			{
				if(wormBodyRadius>r0*(1.00f))
				{
					if(stage==1)
					{
						positionVector = position_cpp + 4 * (pCount+i_start);
						positionVector[ 0 ] = xc + wormBodyRadius*cos(0.0f);
						positionVector[ 1 ] = yc + wormBodyRadius*sin(0.0f);
						positionVector[ 2 ] = zc + r0*j;
						positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow

						positionVector = position_cpp + 4 * (pCount+1+i_start);
						positionVector[ 0 ] = xc - wormBodyRadius*cos(0.0f);
						positionVector[ 1 ] = yc - wormBodyRadius*sin(0.0f);
						positionVector[ 2 ] = zc + r0*j;
						positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow
					}
					pCount+=2;
					currSlice_pCount+=2;
				}
				else
				if(wormBodyRadius<r0*(1.00-0.707))
				{
					if(stage==1)
					{
						positionVector = position_cpp + 4 * (pCount+i_start);
						positionVector[ 0 ] = xc;
						positionVector[ 1 ] = yc;
						positionVector[ 2 ] = zc + r0*j;
						positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow
					}
					pCount++;
					currSlice_pCount++;
				}
			}

			if(wormBodyRadius<r0*0.707f) break;
			alpha = 2*asin(0.5f*r0/wormBodyRadius);//in radians//recalculate -- wormBodyRadius changed
			angle = alpha;

			while(angle<0.89/*radians = less or equal to 51 degrees*/)
			{
				if(stage==1)
				{	
					positionVector = position_cpp + 4 * (pCount+i_start);
					positionVector[ 0 ] = xc + wormBodyRadius*cos(angle);
					positionVector[ 1 ] = yc + wormBodyRadius*sin(angle);
					positionVector[ 2 ] = zc + r0*j;
					positionVector[ 3 ] = 2.2f;/* 2 = elastic matter, green */ 

					positionVector = position_cpp + 4 * (pCount+1+i_start);
					positionVector[ 0 ] = xc + wormBodyRadius*cos(angle);
					positionVector[ 1 ] = yc - wormBodyRadius*sin(angle);
					positionVector[ 2 ] = zc + r0*j;
					positionVector[ 3 ] = 2.2f;/* 2 = elastic matter, green*/ 

					positionVector = position_cpp + 4 * (pCount+2+i_start);
					positionVector[ 0 ] = xc - wormBodyRadius*cos(angle);
					positionVector[ 1 ] = yc + wormBodyRadius*sin(angle);
					positionVector[ 2 ] = zc + r0*j;
					positionVector[ 3 ] = 2.2f;/* 2 = elastic matter, green */ 

					positionVector = position_cpp + 4 * (pCount+3+i_start);
					positionVector[ 0 ] = xc - wormBodyRadius*cos(angle);
					positionVector[ 1 ] = yc - wormBodyRadius*sin(angle);
					positionVector[ 2 ] = zc + r0*j;
					positionVector[ 3 ] = 2.2f;/* 2 = elastic matter, green*/ 
				}
				pCount += 4;
				currSlice_pCount +=4;
				angle+= alpha;
			}

			//if(elasticLayers==1)
			{
				angle-= alpha;//step back for 1 radial segment
				float non_muscle_angle = PI - 2.f*angle;
				int n_non_muscle_particles = (int)(floor(non_muscle_angle / alpha)-1);// distance between each 2 radially adjacent particles will be r0 or more (not less); alpha corresponds to r0
				if(n_non_muscle_particles>0)
				{
					float beta = non_muscle_angle / (n_non_muscle_particles+1);
					int nmp_counter = 0;//non muscle particles counter

					for(i=0;i<(n_non_muscle_particles+1)/2;i++)
					{
						angle+= beta;

						if(stage==1)
						{	
							positionVector = position_cpp + 4 * (pCount+i_start);
							positionVector[ 0 ] = xc + wormBodyRadius*cos(angle);
							positionVector[ 1 ] = yc + wormBodyRadius*sin(angle);
							positionVector[ 2 ] = zc + r0*j;
							positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow
						
							positionVector = position_cpp + 4 * (pCount+1+i_start);
							positionVector[ 0 ] = xc + wormBodyRadius*cos(angle);
							positionVector[ 1 ] = yc - wormBodyRadius*sin(angle);
							positionVector[ 2 ] = zc + r0*j;
							positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow
						}

						pCount += 2;
						currSlice_pCount += 2;
						nmp_counter += 2;

						if(nmp_counter/2==n_non_muscle_particles)
						{
							break;
						}

						if(stage==1)
						{	
							positionVector = position_cpp + 4 * (pCount+i_start);
							positionVector[ 0 ] = xc - wormBodyRadius*cos(angle);
							positionVector[ 1 ] = yc + wormBodyRadius*sin(angle);
							positionVector[ 2 ] = zc + r0*j;
							positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow
						
							positionVector = position_cpp + 4 * (pCount+1+i_start);
							positionVector[ 0 ] = xc - wormBodyRadius*cos(angle);
							positionVector[ 1 ] = yc - wormBodyRadius*sin(angle);
							positionVector[ 2 ] = zc + r0*j;
							positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow
						}

						pCount += 2;
						currSlice_pCount += 2;
						nmp_counter += 2;
					}
				}

				//1-st, outer layer (outer worm shell) was just created (I mean mass points composing it)
				//now it's time for membrane

				if(elasticLayers==1)
				{///////////
					if(stage==0)
					{
						if((j==jmin)||(j==jmax))
						{
							numOfMembranes+= currSlice_pCount;//currently we plan to build membrane over only outer worm shell surface
							if((j==jmin)&&(currSlice_pCount==4)) numOfMembranes+=2;
							if((j==jmax)&&(currSlice_pCount==6)) numOfMembranes+=6;
						}
						else
							numOfMembranes+= currSlice_pCount*2;//all inner structures of the worm will not be covered by membranes in actual version
					}
					else
					{
						//here, at stage = 1, memory for membraneData_cpp (numOfMembranes) is already allocated,
						//and we fill it with i-j-k indexes of membranes which we build here:
						if((j==jmin)&&(currSlice_pCount==4))
						{
							membraneData_cpp [mc*3+0] = 0+i_start;
							membraneData_cpp [mc*3+1] = 1+i_start;
							membraneData_cpp [mc*3+2] = 2+i_start;
							mc++;
							membraneData_cpp [mc*3+0] = 0+i_start;
							membraneData_cpp [mc*3+1] = 1+i_start;
							membraneData_cpp [mc*3+2] = 3+i_start;
							mc++;
						}

						if((j==jmax)&&(currSlice_pCount==6))
						{
							membraneData_cpp [mc*3+0] = currSlice_start+0+i_start;
							membraneData_cpp [mc*3+1] = currSlice_start+2+i_start;
							membraneData_cpp [mc*3+2] = currSlice_start+6+i_start;
							mc++;
							membraneData_cpp [mc*3+0] = currSlice_start+0+i_start;
							membraneData_cpp [mc*3+1] = currSlice_start+3+i_start;
							membraneData_cpp [mc*3+2] = currSlice_start+6+i_start;
							mc++;

							membraneData_cpp [mc*3+0] = currSlice_start+2+i_start;
							membraneData_cpp [mc*3+1] = currSlice_start+4+i_start;
							membraneData_cpp [mc*3+2] = currSlice_start+6+i_start;
							mc++;
							membraneData_cpp [mc*3+0] = currSlice_start+3+i_start;
							membraneData_cpp [mc*3+1] = currSlice_start+5+i_start;
							membraneData_cpp [mc*3+2] = currSlice_start+6+i_start;
							mc++;

							membraneData_cpp [mc*3+0] = currSlice_start+1+i_start;
							membraneData_cpp [mc*3+1] = currSlice_start+4+i_start;
							membraneData_cpp [mc*3+2] = currSlice_start+6+i_start;
							mc++;
							membraneData_cpp [mc*3+0] = currSlice_start+1+i_start;
							membraneData_cpp [mc*3+1] = currSlice_start+5+i_start;
							membraneData_cpp [mc*3+2] = currSlice_start+6+i_start;
							mc++;
						}

						if(j>jmin)
						{
							int ii,jj,kk=0;
							float middle_pos_x;
							float middle_pos_y;
							float middle_pos_z;
							float dist;//distance between middle of [ii]-[jj] and [kk]
							float dist_min;
							int q,w;
							float *pii;
							float *pjj;
							float *pkk;

							// ii and jj on prevSlice, kk - on currSlice
							// 11111111111111111111111111111111111111111111111111111111111111111
							for(q=0;q<prevSlice_pCount;q++)
							{
								if(q==0) { ii=prevSlice_start+0; jj=prevSlice_start+2; } else
								if(q==1) { ii=prevSlice_start+0; jj=prevSlice_start+3; } else
								if(q==2) { ii=prevSlice_start+1; jj=prevSlice_start+4; } else
								if(q==3) { ii=prevSlice_start+1; jj=prevSlice_start+5; } else
										 { ii=prevSlice_start+q-2; jj=prevSlice_start+q+2*(q+2<prevSlice_pCount); }

								if(prevSlice_pCount==4)//head or tail tip
								{
									if(q==0) { ii=prevSlice_start+0; jj=prevSlice_start+2; } else
									if(q==1) { ii=prevSlice_start+0; jj=prevSlice_start+3; } else
									if(q==2) { ii=prevSlice_start+1; jj=prevSlice_start+2; } else
									if(q==3) { ii=prevSlice_start+1; jj=prevSlice_start+3; }
								}

								ii+= i_start;
								jj+= i_start;
								
								pii = position_cpp + 4 * (ii);
								pjj = position_cpp + 4 * (jj);

								middle_pos_x = (pii[0]+pjj[0])/2.f;
								middle_pos_y = (pii[1]+pjj[1])/2.f;
								middle_pos_z = (pii[2]+pjj[2])/2.f;

								dist_min = 10*r0;
									
								for(w=0;w<currSlice_pCount;w++)
								{
									pkk = position_cpp + 4 * (currSlice_start+w+i_start);
									dist = sqrt( (middle_pos_x - pkk[0])*(middle_pos_x - pkk[0])+
												 (middle_pos_y - pkk[1])*(middle_pos_y - pkk[1])+
												 (middle_pos_z - pkk[2])*(middle_pos_z - pkk[2]) );
									if(dist<=dist_min)//!!! "<=" here and "<" at the similar place below is necessary
									{
										dist_min = dist;
										kk = currSlice_start+w+i_start;
									}
								}

								membraneData_cpp [mc*3+0] = ii;//i;
								membraneData_cpp [mc*3+1] = jj;//array_j[j];
								membraneData_cpp [mc*3+2] = kk;//array_j[k];
								mc++;
							}
							// 11111111111111111111111111111111111111111111111111111111111111111

							
							// ii and jj on currSlice, kk - on prevSlice
							// 22222222222222222222222222222222222222222222222222222222222222222
							for(q=0;q<currSlice_pCount;q++)
							{
								if(q==0) { ii=currSlice_start+0; jj=currSlice_start+2; } else
								if(q==1) { ii=currSlice_start+0; jj=currSlice_start+3; } else
								if(q==2) { ii=currSlice_start+1; jj=currSlice_start+4; } else
								if(q==3) { ii=currSlice_start+1; jj=currSlice_start+5; } else
										 { ii=currSlice_start+q-2; jj=currSlice_start+q+2*(q+2<currSlice_pCount); }

								if(currSlice_pCount==4)//head or tail tip
								{
									if(q==0) { ii=currSlice_start+0; jj=currSlice_start+2; } else
									if(q==1) { ii=currSlice_start+0; jj=currSlice_start+3; } else
									if(q==2) { ii=currSlice_start+1; jj=currSlice_start+2; } else
									if(q==3) { ii=currSlice_start+1; jj=currSlice_start+3; }
								}

								ii+= i_start;
								jj+= i_start;
								
								pii = position_cpp + 4 * (ii);
								pjj = position_cpp + 4 * (jj);

								middle_pos_x = (pii[0]+pjj[0])/2.f;
								middle_pos_y = (pii[1]+pjj[1])/2.f;
								middle_pos_z = (pii[2]+pjj[2])/2.f;

								dist_min = 10*r0;
									
								for(w=0;w<prevSlice_pCount;w++)
								{
									pkk = position_cpp + 4 * (prevSlice_start+w+i_start);
									dist = sqrt( (middle_pos_x - pkk[0])*(middle_pos_x - pkk[0])+
												 (middle_pos_y - pkk[1])*(middle_pos_y - pkk[1])+
												 (middle_pos_z - pkk[2])*(middle_pos_z - pkk[2]) );
									if(dist<dist_min)//!!! 
									{
										dist_min = dist;
										kk = prevSlice_start+w+i_start;
									}
								}

								membraneData_cpp [mc*3+0] = ii;//i;
								membraneData_cpp [mc*3+1] = jj;//array_j[j];
								membraneData_cpp [mc*3+2] = kk;//array_j[k];
								mc++;
							}
							// 22222222222222222222222222222222222222222222222222222222222222222
						}
					}

					prevSlice_pCount = currSlice_pCount;
					prevSlice_start = currSlice_start;
				}//////////
			}

			wormBodyRadius -= r0;
			elasticLayers++;
		}
	}////////////////////////////////////////////////////

	if(stage==1)
	{
		for(i=0;i<pCount;i++)
		{
			velocityVector = velocity_cpp + 4 * (i+i_start);	
			velocityVector[ 0 ] = 0;
			velocityVector[ 1 ] = 0;
			velocityVector[ 2 ] = 0;
			velocityVector[ 3 ] = 0;
		}
	}

	//if(stage==0) numOfMembranes /= 2;//remove this line when development of this function is complete

	return pCount;

	// makeworm
}

int generateInnerWormLiquid(int stage, int i_start,float *position_cpp, float *velocity_cpp, owConfigProrerty * config)
{
	//return 0;
	float alpha;// = 2.f*3.14159f/segmentsCount;
	float wormBodyRadius;// = h*coeff / sin(alpha/2);
	int pCount = 0;//particle counter
	int i;

	float *positionVector;
	float *velocityVector;
	int elasticLayers;//starting from 2, because 1 is for outer shell and doesn't contain liquid particles
	float xc = config->xmax*0.5f;
	float yc = config->ymax*0.3f;
	float zc = config->zmax*0.5f;
	float PI = 3.1415926536f;
	float angle;
	float x,y,z;

	// makeworm
	// outer worm shell elastic cylinder generation

	pCount = 0;
	
	float jmin=-100.f,jmax=100.f;
	//float jmin=-10.f,jmax=10.f;

	float j;


	for(j=jmin;j<=jmax;j+=0.85f)	// number of cross-slices of the worm, in direction from head to tail. -
							// 1..-1 for example will give 3 slices in the middle of the worm
	{////////////////////////////////////////////////////

		elasticLayers = 2;
		wormBodyRadius = 6.0f*r0*sqrt(max(1.f-(1.0e-4f)*j*j,0.f)) - r0*(1+0.85f);

		while(1)
		{
			if(wormBodyRadius>0.707f*r0)
			{
				if(stage==1)
				{
					positionVector = position_cpp + 4 * (pCount+i_start);
					positionVector[ 0 ] = xc + wormBodyRadius*sin(0.0f);
					positionVector[ 1 ] = yc + wormBodyRadius*cos(0.0f);
					positionVector[ 2 ] = zc + r0*j;
					positionVector[ 3 ] = 1.1f;// liquid

					positionVector = position_cpp + 4 * (pCount+1+i_start);
					positionVector[ 0 ] = xc - wormBodyRadius*sin(0.0f);
					positionVector[ 1 ] = yc - wormBodyRadius*cos(0.0f);
					positionVector[ 2 ] = zc + r0*j;
					positionVector[ 3 ] = 1.1f;// liquid
				}
				pCount+=2;
			}
			else 
			{
				/*
				if(wormBodyRadius>0.3*r0)
				{
					if(stage==1)
					{
						positionVector = position_cpp + 4 * (pCount+i_start);
						positionVector[ 0 ] = xc;
						positionVector[ 1 ] = yc;
						positionVector[ 2 ] = zc + r0*j;
						positionVector[ 3 ] = 1.1f;// liquid
					}
					pCount++;
				}*/
				break;
			}

			alpha = 2*asin(0.5f*r0/wormBodyRadius);//in radians//recalculate -- wormBodyRadius changed
			
			angle = 0;

			if(elasticLayers == 1) 
			{
				angle = alpha;

				if(alpha>0)
				while(angle<0.89)
				{
					angle+= alpha;
				}
				angle-= alpha;//step back for 1 radial segment
			}

			float non_muscle_angle = PI - 2.f*angle;
			int n_non_muscle_particles = (int)(floor(non_muscle_angle / (alpha*0.85f) )-1);
			float beta = non_muscle_angle / (n_non_muscle_particles+1);

			for(i=0;i<n_non_muscle_particles;i++)
			{
				angle+= beta;

				if(stage==1)
				{	
					positionVector = position_cpp + 4 * (pCount+i_start);
					positionVector[ 0 ] = xc + wormBodyRadius*sin(angle);
					positionVector[ 1 ] = yc + wormBodyRadius*cos(angle);
					positionVector[ 2 ] = zc + r0*j;
					positionVector[ 3 ] = 1.1f;// liquid
				
					positionVector = position_cpp + 4 * (pCount+1+i_start);
					positionVector[ 0 ] = xc - wormBodyRadius*sin(angle);
					positionVector[ 1 ] = yc + wormBodyRadius*cos(angle);
					positionVector[ 2 ] = zc + r0*j;
					positionVector[ 3 ] = 1.1f;// liquid
				}

				pCount += 2;
			}
			
			elasticLayers++;
			wormBodyRadius -= r0*0.85f;
		}
	}

	//and here we add outer liquid for worm swimming
/**/
	for(x=3*r0;x<config->xmax-3*r0;x+=r0)
	{
		for(y=3*r0;y<config->ymax*0.15;y+=r0)
		{
			for(z=3*r0;z<config->zmax-3*r0;z+=r0)
			{
				if(stage==1)
				{	
					positionVector = position_cpp + 4 * (pCount+i_start);
					positionVector[ 0 ] = x;
					positionVector[ 1 ] = y;
					positionVector[ 2 ] = z;
					positionVector[ 3 ] = 1.1f;// liquid
				}

				pCount++;
			}
		}
	}
/**/

	if(stage==1)
	{
		for(i=0;i<pCount;i++)
		{
			velocityVector = velocity_cpp + 4 * (i+i_start);	
			velocityVector[ 0 ] = 0;
			velocityVector[ 1 ] = 0;
			velocityVector[ 2 ] = 0;
			velocityVector[ 3 ] = 0;
		}
	}

	return pCount;
}


void owHelper::generateConfiguration(int stage, float *position_cpp, float *velocity_cpp, float *& elasticConnectionsData_cpp, int *membraneData_cpp, int & numOfLiquidP, int & numOfElasticP, int & numOfBoundaryP, int & numOfElasticConnections, int & numOfMembranes, int *particleMembranesList_cpp, owConfigProrerty * config)
{
	float p_type = LIQUID_PARTICLE;
	int i = 0;// particle counter
	int ix,iy,iz;
	int ecc = 0;//elastic connections counter

	int nx = (int)( ( config->xmax - config->xmin ) / r0 ); //X
	int ny = (int)( ( config->ymax - config->ymin ) / r0 ); //Y
	int nz = (int)( ( config->zmax - config->zmin ) / r0 ); //Z

	int numOfMembraneParticles = generateWormShell(0,0,position_cpp,velocity_cpp, numOfMembranes, membraneData_cpp, config);

	if(stage==0)
	{
		numOfLiquidP = generateInnerWormLiquid(0,0,position_cpp,velocity_cpp, config);
		numOfElasticP = numOfMembraneParticles;
		numOfBoundaryP = 0;

		if(numOfElasticP<=0) elasticConnectionsData_cpp = NULL; else elasticConnectionsData_cpp = new float[ 4 * numOfElasticP * MAX_NEIGHBOR_COUNT ];
	}

	//=============== create worm body (elastic parts) ==================================================
	if(stage==1)
	{
		i += generateWormShell(1/*stage*/,i,position_cpp,velocity_cpp, numOfMembranes,membraneData_cpp, config);
		//initialize elastic connections data structure (with NO_PARTICLE_ID values)
		for(int ii = 0; ii < numOfElasticP * MAX_NEIGHBOR_COUNT; ii++)
		{
			ecc = 0;

			elasticConnectionsData_cpp[ 4 * ii + 0 ] = NO_PARTICLE_ID;
			elasticConnectionsData_cpp[ 4 * ii + 1 ] = 0;
			elasticConnectionsData_cpp[ 4 * ii + 2 ] = 0;
			elasticConnectionsData_cpp[ 4 * ii + 3 ] = 0; 
		}
	}


	//=============== create worm body (inner liquid) ==================================================
	if(stage==1)
	{
		i += generateInnerWormLiquid(stage,i,position_cpp,velocity_cpp, config);
	}


	if(stage==0) 
	{
		numOfBoundaryP = 2 * ( nx*ny + (nx+ny-2)*(nz-2) ); 
	}
	else
	if(stage==1)
	{
		//===================== create boundary particles ==========================================================
		p_type = BOUNDARY_PARTICLE;
		
		// 1 - top and bottom 
		for(ix=0;ix<nx;ix++)
		{
			for(iy=0;iy<ny;iy++)
			{
				if( ((ix==0)||(ix==nx-1)) || ((iy==0)||(iy==ny-1)) )
				{
					if( ((ix==0)||(ix==nx-1)) && ((iy==0)||(iy==ny-1)) )//corners
					{
						position_cpp[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position_cpp[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position_cpp[ 4 * i + 2 ] =  0*r0 + r0/2;//z
						position_cpp[ 4 * i + 3 ] = p_type;
						velocity_cpp[ 4 * i + 0 ] = (1.f*(ix==0)-1*(ix==nx-1))/sqrt(3.f);//norm x
						velocity_cpp[ 4 * i + 1 ] = (1.f*(iy==0)-1*(iy==ny-1))/sqrt(3.f);//norm y
						velocity_cpp[ 4 * i + 2 ] =  1.f/sqrt(3.f);//norm z
						velocity_cpp[ 4 * i + 3 ] = p_type;
						i++;
						position_cpp[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position_cpp[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position_cpp[ 4 * i + 2 ] = (nz-1)*r0 + r0/2;//z
						position_cpp[ 4 * i + 3 ] = p_type;
						velocity_cpp[ 4 * i + 0 ] = (1*(ix==0)-1*(ix==nx-1))/sqrt(3.f);//norm x
						velocity_cpp[ 4 * i + 1 ] = (1*(iy==0)-1*(iy==ny-1))/sqrt(3.f);//norm y
						velocity_cpp[ 4 * i + 2 ] = -1.f/sqrt(3.f);//norm z
						velocity_cpp[ 4 * i + 3 ] = p_type;
						i++;
					}
					else //edges
					{
						position_cpp[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position_cpp[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position_cpp[ 4 * i + 2 ] =  0*r0 + r0/2;//z
						position_cpp[ 4 * i + 3 ] = p_type;
						velocity_cpp[ 4 * i + 0 ] =  1.f*((ix==0)-(ix==nx-1))/sqrt(2.f);//norm x
						velocity_cpp[ 4 * i + 1 ] =  1.f*((iy==0)-(iy==ny-1))/sqrt(2.f);//norm y
						velocity_cpp[ 4 * i + 2 ] =  1.f/sqrt(2.f);//norm z
						velocity_cpp[ 4 * i + 3 ] = p_type;
						i++;
						position_cpp[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position_cpp[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position_cpp[ 4 * i + 2 ] = (nz-1)*r0 + r0/2;//z
						position_cpp[ 4 * i + 3 ] = p_type;
						velocity_cpp[ 4 * i + 0 ] =  1.f*((ix==0)-(ix==nx-1))/sqrt(2.f);//norm x
						velocity_cpp[ 4 * i + 1 ] =  1.f*((iy==0)-(iy==ny-1))/sqrt(2.f);//norm y
						velocity_cpp[ 4 * i + 2 ] = -1.f/sqrt(2.f);//norm z
						velocity_cpp[ 4 * i + 3 ] = p_type;
						i++;
					}
				}
				else //planes
				{
						position_cpp[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position_cpp[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position_cpp[ 4 * i + 2 ] =  0*r0 + r0/2;//z
						position_cpp[ 4 * i + 3 ] = p_type;
						velocity_cpp[ 4 * i + 0 ] =  0;//norm x
						velocity_cpp[ 4 * i + 1 ] =  0;//norm y
						velocity_cpp[ 4 * i + 2 ] =  1;//norm z
						velocity_cpp[ 4 * i + 3 ] = p_type;
						i++;
						position_cpp[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position_cpp[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position_cpp[ 4 * i + 2 ] = (nz-1)*r0 + r0/2;//z
						position_cpp[ 4 * i + 3 ] = p_type;
						velocity_cpp[ 4 * i + 0 ] =  0;//norm x
						velocity_cpp[ 4 * i + 1 ] =  0;//norm y
						velocity_cpp[ 4 * i + 2 ] = -1;//norm z
						velocity_cpp[ 4 * i + 3 ] = p_type;
						i++;
				}
			}
		}

		// 2 - side walls OX-OZ and opposite
		for(ix=0;ix<nx;ix++)
		{
			for(iz=1;iz<nz-1;iz++)
			{
				//edges
				if((ix==0)||(ix==nx-1))
				{
					position_cpp[ 4 * i + 0 ] = ix*r0 + r0/2;//x
					position_cpp[ 4 * i + 1 ] =  0*r0 + r0/2;//y
					position_cpp[ 4 * i + 2 ] = iz*r0 + r0/2;//z
					position_cpp[ 4 * i + 3 ] = p_type;
					velocity_cpp[ 4 * i + 0 ] =  0;//norm x
					velocity_cpp[ 4 * i + 1 ] =  1.f/sqrt(2.f);//norm y
					velocity_cpp[ 4 * i + 2 ] =  1.f*((iz==0)-(iz==nz-1))/sqrt(2.f);//norm z
					velocity_cpp[ 4 * i + 3 ] = p_type;
					i++;
					position_cpp[ 4 * i + 0 ] = ix*r0 + r0/2;//x
					position_cpp[ 4 * i + 1 ] = (ny-1)*r0 + r0/2;//y
					position_cpp[ 4 * i + 2 ] = iz*r0 + r0/2;//z
					position_cpp[ 4 * i + 3 ] = p_type;
					velocity_cpp[ 4 * i + 0 ] =  0;//norm x
					velocity_cpp[ 4 * i + 1 ] = -1.f/sqrt(2.f);//norm y
					velocity_cpp[ 4 * i + 2 ] =  1.f*((iz==0)-(iz==nz-1))/sqrt(2.f);//norm z
					velocity_cpp[ 4 * i + 3 ] = p_type;
					i++;
				}
				else //planes
				{
					position_cpp[ 4 * i + 0 ] = ix*r0 + r0/2;//x
					position_cpp[ 4 * i + 1 ] =  0*r0 + r0/2;//y
					position_cpp[ 4 * i + 2 ] = iz*r0 + r0/2;//z
					position_cpp[ 4 * i + 3 ] = p_type;
					velocity_cpp[ 4 * i + 0 ] =  0;//norm x
					velocity_cpp[ 4 * i + 1 ] =  1;//norm y
					velocity_cpp[ 4 * i + 2 ] =  0;//norm z
					velocity_cpp[ 4 * i + 3 ] = p_type;
					i++;
					position_cpp[ 4 * i + 0 ] = ix*r0 + r0/2;//x
					position_cpp[ 4 * i + 1 ] = (ny-1)*r0 + r0/2;//y
					position_cpp[ 4 * i + 2 ] = iz*r0 + r0/2;//z
					position_cpp[ 4 * i + 3 ] = p_type;
					velocity_cpp[ 4 * i + 0 ] =  0;//norm x
					velocity_cpp[ 4 * i + 1 ] = -1;//norm y
					velocity_cpp[ 4 * i + 2 ] =  0;//norm z
					velocity_cpp[ 4 * i + 3 ] = p_type;
					i++;
				}
			}
		}

		// 3 - side walls OY-OZ and opposite
		for(iy=1;iy<ny-1;iy++)
		{
			for(iz=1;iz<nz-1;iz++)
			{
				position_cpp[ 4 * i + 0 ] =  0*r0 + r0/2;//x
				position_cpp[ 4 * i + 1 ] = iy*r0 + r0/2;//y
				position_cpp[ 4 * i + 2 ] = iz*r0 + r0/2;//z
				position_cpp[ 4 * i + 3 ] = p_type;
				velocity_cpp[ 4 * i + 0 ] =  1;//norm x
				velocity_cpp[ 4 * i + 1 ] =  0;//norm y
				velocity_cpp[ 4 * i + 2 ] =  0;//norm z
				velocity_cpp[ 4 * i + 3 ] = p_type;
				i++;
				position_cpp[ 4 * i + 0 ] = (nx-1)*r0 + r0/2;//x
				position_cpp[ 4 * i + 1 ] = iy*r0 + r0/2;//y
				position_cpp[ 4 * i + 2 ] = iz*r0 + r0/2;//z
				position_cpp[ 4 * i + 3 ] = p_type;
				velocity_cpp[ 4 * i + 0 ] = -1;//norm x
				velocity_cpp[ 4 * i + 1 ] =  0;//norm y
				velocity_cpp[ 4 * i + 2 ] =  0;//norm z
				velocity_cpp[ 4 * i + 3 ] = p_type;
				i++;
			}
		}
	}

	if(stage==0)
	{
		config->setParticleCount(numOfLiquidP + numOfBoundaryP + numOfElasticP);

		if(config->getParticleCount()<=0)
		{
			printf("\nWarning! Generated scene contains %d particles!\n",config->getParticleCount());
			exit(-2);
		}
	}
	else
	if(stage==1)
	{
		if(config->getParticleCount()!=i)
		{
			printf("\nWarning! Preliminary [%d] and final [%d] particle count are different\n",config->getParticleCount(),i);
			exit(-4);
		}

		
		for(int mli = 0/*membrane list index*/; mli<numOfElasticP; mli++)
		{
			for(int sli = 0 /*sublist index*/; sli<MAX_MEMBRANES_INCLUDING_SAME_PARTICLE; sli++)
			{
				particleMembranesList_cpp [mli*MAX_MEMBRANES_INCLUDING_SAME_PARTICLE + sli] = -1;// no membranes connected with current particle
			}
		}

		///////////////debug////////////
		int j;
		int ecc_total = 0;
		//float ix,iy,iz,jx,jy,jz;
		//int j_count=0;
		int muscleCounter = 0;
		int m_index[10640];
		float WXC = config->xmax*0.5f;
		float WYC = config->ymax*0.3f;
		float WZC = config->zmax*0.5f;
		//int sm_cnt = 0;
		//int array_k[MAX_NEIGHBOR_COUNT];
		for(i=numOfElasticP-numOfMembraneParticles;i<numOfElasticP;i++)
		{
			float dx2,dy2,dz2,r2_ij,r_ij;
			int q_i_start;
			int dq;//dorsal quadrant - "+1"=right, "-1"=left
			float muscle_color = 0.1f;
			ecc = 0;//!important!
			//        _____1_______      2       _____3________       
			for(j=0;j<numOfElasticP+numOfLiquidP+numOfBoundaryP;j++)
			{
				if(i!=j)
				{
					if(j==numOfElasticP) j+= numOfLiquidP;//skip liquid particles (they are located in the middle of memory) as candidates to spring connections


					dx2 = (position_cpp[ 4 * i + 0 ] - position_cpp[ 4 * j + 0 ]); dx2 *= dx2;
					dy2 = (position_cpp[ 4 * i + 1 ] - position_cpp[ 4 * j + 1 ]); dy2 *= dy2;
					dz2 = (position_cpp[ 4 * i + 2 ] - position_cpp[ 4 * j + 2 ]); dz2 *= dz2;
					r2_ij = dx2 + dy2 + dz2;
					r_ij = (float)sqrt(r2_ij);

					//if(r_ij<=r0*2*sqrt(/*3.2*/1.7))//grid = 1.5*r0
					if(r_ij<=r0*sqrt(/*2.3*/2.7/*2.7*/))//grid = 1.0*r0
					{
						elasticConnectionsData_cpp[ 4 * ( MAX_NEIGHBOR_COUNT * i + ecc) + 0 ] = ((float)j) + 0.1f;		// index of j-th particle in a pair connected with spring
						elasticConnectionsData_cpp[ 4 * ( MAX_NEIGHBOR_COUNT * i + ecc) + 1 ] = r_ij*simulationScale*0.95f;	// resting distance; that's why we use float type for elasticConnectionsData_cpp
						elasticConnectionsData_cpp[ 4 * ( MAX_NEIGHBOR_COUNT * i + ecc) + 2 ] = 0;						// type of connection; 0 - ordinary spring, 1 - muscle
						elasticConnectionsData_cpp[ 4 * ( MAX_NEIGHBOR_COUNT * i + ecc) + 3 ] = 0;						// not in use yet

						//define muscles
						if((position_cpp[ 4 * i + 2 ]<WZC+r0*95)&&(position_cpp[ 4 * j + 2 ]<WZC+r0*95))
						if((position_cpp[ 4 * i + 2 ]>WZC-r0*92)&&(position_cpp[ 4 * j + 2 ]>WZC-r0*92))
						if( (fabs(position_cpp[4*i+3]-2.2f)<=0.05) && (fabs(position_cpp[4*j+3]-2.2f)<=0.05) )//both points - i and j - are green
						if((dz2>4*dx2)&&(dz2>4*dy2)&&(dx2>4*dy2))
						{
							if(position_cpp[ 4 * i + 0 ]>WXC)
							{
								muscle_color = 1.1f;

								//DR and DL quadrant (dorsal) muscles mapping
								for(dq=-1;dq<=1;dq+=2)//dorsal quadrant - "-1"=right, "+1"=left
								{
									if(dq==1) q_i_start = 0; else q_i_start = 24*3; //muscle quadrant starting index
									
									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*97)&&(position_cpp[4*j+2]<WZC+r0*97)) //MDR01 || MDL01
									if((position_cpp[4*i+2]>WZC+r0*85.9)&&(position_cpp[4*j+2]>WZC+r0*85.9)) muscle_color = q_i_start + 1.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*95.0)&&(position_cpp[4*j+2]<WZC+r0*95.0)) //MDR02 || MDL02
									if((position_cpp[4*i+2]>WZC+r0*83.5)&&(position_cpp[4*j+2]>WZC+r0*83.5)) muscle_color = q_i_start + 2.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*86.5)&&(position_cpp[4*j+2]<WZC+r0*86.5)) //MDR03 || MDL03
									if((position_cpp[4*i+2]>WZC+r0*77.5)&&(position_cpp[4*j+2]>WZC+r0*77.5)) muscle_color = q_i_start + 3.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*84.5)&&(position_cpp[4*j+2]<WZC+r0*84.5)) //MDR04 || MDL04
									if((position_cpp[4*i+2]>WZC+r0*76.5)&&(position_cpp[4*j+2]>WZC+r0*76.5)) muscle_color = q_i_start + 4.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*82.5)&&(position_cpp[4*j+2]<WZC+r0*82.5)) 
									if((position_cpp[4*i+2]>WZC+r0*72.5)&&(position_cpp[4*j+2]>WZC+r0*72.5)) muscle_color = q_i_start + 4.5f;//x.5 = violet
								
									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*78.5)&&(position_cpp[4*j+2]<WZC+r0*78.5)) //MDR05 || MDL05
									if((position_cpp[4*i+2]>WZC+r0*66.9)&&(position_cpp[4*j+2]>WZC+r0*66.9)) muscle_color = q_i_start + 5.2f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*77.5)&&(position_cpp[4*j+2]<WZC+r0*77.5)) 
									if((position_cpp[4*i+2]>WZC+r0*65.9)&&(position_cpp[4*j+2]>WZC+r0*65.9)) muscle_color = q_i_start + 5.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*74.0)&&(position_cpp[4*j+2]<WZC+r0*74.0)) //MDR06 || MDL06
									if((position_cpp[4*i+2]>WZC+r0*55.0)&&(position_cpp[4*j+2]>WZC+r0*55.0)) muscle_color = q_i_start + 6.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*74.0)&&(position_cpp[4*j+2]<WZC+r0*74.0)) 
									if((position_cpp[4*i+2]>WZC+r0*54.5)&&(position_cpp[4*j+2]>WZC+r0*54.5)) muscle_color = q_i_start + 6.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*68.5)&&(position_cpp[4*j+2]<WZC+r0*68.5)) //MDR07 || MDL07
									if((position_cpp[4*i+2]>WZC+r0*51.0)&&(position_cpp[4*j+2]>WZC+r0*51.0)) muscle_color = q_i_start + 7.3f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*66.5)&&(position_cpp[4*j+2]<WZC+r0*66.5)) 
									if((position_cpp[4*i+2]>WZC+r0*49.5)&&(position_cpp[4*j+2]>WZC+r0*49.5)) muscle_color = q_i_start + 7.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*56.5)&&(position_cpp[4*j+2]<WZC+r0*56.5)) //MDR08 || MDL08
									if((position_cpp[4*i+2]>WZC+r0*40.0)&&(position_cpp[4*j+2]>WZC+r0*40.0)) muscle_color = q_i_start + 8.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*55.5)&&(position_cpp[4*j+2]<WZC+r0*55.5)) 
									if((position_cpp[4*i+2]>WZC+r0*38.5)&&(position_cpp[4*j+2]>WZC+r0*38.5)) muscle_color = q_i_start + 8.5f;//x.5 = violet

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*52.1)&&(position_cpp[4*j+2]<WZC+r0*52.1)) //MDR09 || MDL09
									if((position_cpp[4*i+2]>WZC+r0*33.5)&&(position_cpp[4*j+2]>WZC+r0*33.5)) muscle_color = q_i_start + 9.2f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*50.5)&&(position_cpp[4*j+2]<WZC+r0*50.5)) 
									if((position_cpp[4*i+2]>WZC+r0*32.5)&&(position_cpp[4*j+2]>WZC+r0*32.5)) muscle_color = q_i_start + 9.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*41.1)&&(position_cpp[4*j+2]<WZC+r0*41.1)) //MDR10 || MDL10
									if((position_cpp[4*i+2]>WZC+r0*22.5)&&(position_cpp[4*j+2]>WZC+r0*22.5)) muscle_color = q_i_start + 10.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*40.0)&&(position_cpp[4*j+2]<WZC+r0*40.0)) 
									if((position_cpp[4*i+2]>WZC+r0*21.5)&&(position_cpp[4*j+2]>WZC+r0*21.5)) muscle_color = q_i_start + 10.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-4*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-5*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*40.0)&&(position_cpp[4*j+2]<WZC+r0*40.0)) 
									if((position_cpp[4*i+2]>WZC+r0*20.5)&&(position_cpp[4*j+2]>WZC+r0*20.5)) muscle_color = q_i_start + 10.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)        &&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*34.5)&&(position_cpp[4*j+2]<WZC+r0*34.5)) //MDR11 || MDL11
									if((position_cpp[4*i+2]>WZC+r0*15.5)&&(position_cpp[4*j+2]>WZC+r0*15.5)) muscle_color = q_i_start + 11.3f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)   &&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*33.5)&&(position_cpp[4*j+2]<WZC+r0*33.5)) 
									if((position_cpp[4*i+2]>WZC+r0*14.5)&&(position_cpp[4*j+2]>WZC+r0*14.5)) muscle_color = q_i_start + 11.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*23.5)&&(position_cpp[4*j+2]<WZC+r0*23.5)) //MDR12 || MDL12
									if((position_cpp[4*i+2]>WZC+r0* 8.5)&&(position_cpp[4*j+2]>WZC+r0* 8.5)) muscle_color = q_i_start + 12.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*22.5)&&(position_cpp[4*j+2]<WZC+r0*22.5)) 
									if((position_cpp[4*i+2]>WZC+r0* 7.5)&&(position_cpp[4*j+2]>WZC+r0* 7.5)) muscle_color = q_i_start + 12.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-4*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-5*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*21.5)&&(position_cpp[4*j+2]<WZC+r0*21.5)) 
									if((position_cpp[4*i+2]>WZC+r0* 6.5)&&(position_cpp[4*j+2]>WZC+r0* 6.5)) muscle_color = q_i_start + 12.5f;//x.5 = violet

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*16.5)&&(position_cpp[4*j+2]<WZC+r0*16.5)) //MDR13 || MDL13
									if((position_cpp[4*i+2]>WZC+r0* 1.5)&&(position_cpp[4*j+2]>WZC+r0* 1.5)) muscle_color = q_i_start + 13.2f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*15.5)&&(position_cpp[4*j+2]<WZC+r0*15.5)) 
									if((position_cpp[4*i+2]>WZC+r0* 0.5)&&(position_cpp[4*j+2]>WZC+r0* 0.5)) muscle_color = q_i_start + 13.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 9.0)&&(position_cpp[4*j+2]<WZC+r0* 9.0)) //MDR14 || MDL14
									if((position_cpp[4*i+2]>WZC-r0* 2.5)&&(position_cpp[4*j+2]>WZC-r0* 2.5)) muscle_color = q_i_start + 14.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 8.5)&&(position_cpp[4*j+2]<WZC+r0* 8.5)) 
									if((position_cpp[4*i+2]>WZC-r0* 3.5)&&(position_cpp[4*j+2]>WZC-r0* 3.5)) muscle_color = q_i_start + 14.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-4*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-5*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 7.5)&&(position_cpp[4*j+2]<WZC+r0* 7.5)) 
									if((position_cpp[4*i+2]>WZC-r0* 4.5)&&(position_cpp[4*j+2]>WZC-r0* 4.5)) muscle_color = q_i_start + 14.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 2.0)&&(position_cpp[4*j+2]<WZC+r0* 2.0)) //MDR15 || MDL15
									if((position_cpp[4*i+2]>WZC-r0*14.5)&&(position_cpp[4*j+2]>WZC-r0*14.5)) muscle_color = q_i_start + 15.3f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)   &&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 1.5)&&(position_cpp[4*j+2]<WZC+r0* 1.5)) 
									if((position_cpp[4*i+2]>WZC-r0*15.5)&&(position_cpp[4*j+2]>WZC-r0*15.5)) muscle_color = q_i_start + 15.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0* 1.5)&&(position_cpp[4*j+2]<WZC-r0* 1.5)) //MDR16 || MDL16
									if((position_cpp[4*i+2]>WZC-r0*21.5)&&(position_cpp[4*j+2]>WZC-r0*21.5)) muscle_color = q_i_start + 16.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0* 2.5)&&(position_cpp[4*j+2]<WZC-r0* 2.5)) 
									if((position_cpp[4*i+2]>WZC-r0*22.5)&&(position_cpp[4*j+2]>WZC-r0*22.5)) muscle_color = q_i_start + 16.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-4*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-5*r0)) 
									if((position_cpp[4*i+2]<WZC-r0* 3.5)&&(position_cpp[4*j+2]<WZC-r0* 3.5)) 
									if((position_cpp[4*i+2]>WZC-r0*23.5)&&(position_cpp[4*j+2]>WZC-r0*23.5)) muscle_color = q_i_start + 16.5f;//x.5 = violet

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*14.0)&&(position_cpp[4*j+2]<WZC-r0*14.0)) //MDR17 || MDL17
									if((position_cpp[4*i+2]>WZC-r0*34.5)&&(position_cpp[4*j+2]>WZC-r0*34.5)) muscle_color = q_i_start + 17.2f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*14.7)&&(position_cpp[4*j+2]<WZC-r0*14.7)) 
									if((position_cpp[4*i+2]>WZC-r0*35.5)&&(position_cpp[4*j+2]>WZC-r0*35.5)) muscle_color = q_i_start + 17.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*20.0)&&(position_cpp[4*j+2]<WZC-r0*20.0)) //MDR18 || MDL18
									if((position_cpp[4*i+2]>WZC-r0*40.5)&&(position_cpp[4*j+2]>WZC-r0*40.5)) muscle_color = q_i_start + 18.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*21.5)&&(position_cpp[4*j+2]<WZC-r0*21.5)) 
									if((position_cpp[4*i+2]>WZC-r0*41.5)&&(position_cpp[4*j+2]>WZC-r0*41.5)) muscle_color = q_i_start + 18.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-4*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-5*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*22.5)&&(position_cpp[4*j+2]<WZC-r0*22.5)) 
									if((position_cpp[4*i+2]>WZC-r0*34.5)&&(position_cpp[4*j+2]>WZC-r0*34.5)) muscle_color = q_i_start + 18.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*34.0)&&(position_cpp[4*j+2]<WZC-r0*34.0)) //MDR19 || MDL19
									if((position_cpp[4*i+2]>WZC-r0*54.5)&&(position_cpp[4*j+2]>WZC-r0*54.5)) muscle_color = q_i_start + 19.3f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*34.5)&&(position_cpp[4*j+2]<WZC-r0*34.5)) 
									if((position_cpp[4*i+2]>WZC-r0*55.5)&&(position_cpp[4*j+2]>WZC-r0*55.5)) muscle_color = q_i_start + 19.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*39.5)&&(position_cpp[4*j+2]<WZC-r0*39.5)) //MDR20 || MDL20
									if((position_cpp[4*i+2]>WZC-r0*50.5)&&(position_cpp[4*j+2]>WZC-r0*50.5)) muscle_color = q_i_start + 20.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*40.5)&&(position_cpp[4*j+2]<WZC-r0*40.5)) 
									if((position_cpp[4*i+2]>WZC-r0*51.5)&&(position_cpp[4*j+2]>WZC-r0*51.5)) muscle_color = q_i_start + 20.5f;//x.5 = violet

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*53.0)&&(position_cpp[4*j+2]<WZC-r0*53.0)) //MDR21 || MDL21
									if((position_cpp[4*i+2]>WZC-r0*71.5)&&(position_cpp[4*j+2]>WZC-r0*71.5)) muscle_color = q_i_start + 21.2f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*54.0)&&(position_cpp[4*j+2]<WZC-r0*54.0)) 
									if((position_cpp[4*i+2]>WZC-r0*72.5)&&(position_cpp[4*j+2]>WZC-r0*72.5)) muscle_color = q_i_start + 21.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*50.0)&&(position_cpp[4*j+2]<WZC-r0*50.0)) //MDR22 || MDL22
									if((position_cpp[4*i+2]>WZC-r0*63.5)&&(position_cpp[4*j+2]>WZC-r0*63.5)) muscle_color = q_i_start + 22.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*50.5)&&(position_cpp[4*j+2]<WZC-r0*50.5)) 
									if((position_cpp[4*i+2]>WZC-r0*64.5)&&(position_cpp[4*j+2]>WZC-r0*64.5)) muscle_color = q_i_start + 22.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*70.0)&&(position_cpp[4*j+2]<WZC-r0*70.0)) //MDR23 || MDL23
									if((position_cpp[4*i+2]>WZC-r0*92.0)&&(position_cpp[4*j+2]>WZC-r0*92.0)) muscle_color = q_i_start + 23.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*71.5)&&(position_cpp[4*j+2]<WZC-r0*71.5)) //MDR24 || MDL24
									if((position_cpp[4*i+2]>WZC-r0*92.0)&&(position_cpp[4*j+2]>WZC-r0*92.0)) muscle_color = q_i_start + 24.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*62.5)&&(position_cpp[4*j+2]<WZC-r0*62.5)) 
									if((position_cpp[4*i+2]>WZC-r0*82.5)&&(position_cpp[4*j+2]>WZC-r0*82.5)) muscle_color = q_i_start + 24.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*63.5)&&(position_cpp[4*j+2]<WZC-r0*63.5)) 
									if((position_cpp[4*i+2]>WZC-r0*66.5)&&(position_cpp[4*j+2]>WZC-r0*66.5)) muscle_color = q_i_start + 24.5f;//x.5 = violet
								}

								elasticConnectionsData_cpp[ m_index[muscleCounter] = 4 * ( MAX_NEIGHBOR_COUNT * i + ecc) + 2 ] = muscle_color;// type of connection; 0 - ordinary spring, 1 or more - muscle
								muscleCounter++;
							}
							else
							{
								muscle_color = 1.1f;

								//VR and VL quadrant (ventral) muscles mapping
								for(dq=-1;dq<=1;dq+=2)//dorsal quadrant - "-1"=right, "+1"=left
								{
									if(dq==1) q_i_start = 24; else q_i_start = 24*2; //muscle quadrant starting index
									
									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*97)&&(position_cpp[4*j+2]<WZC+r0*97)) //MVR01 || MVL01
									if((position_cpp[4*i+2]>WZC+r0*85.9)&&(position_cpp[4*j+2]>WZC+r0*85.9)) muscle_color = q_i_start + 1.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*95.0)&&(position_cpp[4*j+2]<WZC+r0*95.0)) //MVR02 || MVL02
									if((position_cpp[4*i+2]>WZC+r0*83.5)&&(position_cpp[4*j+2]>WZC+r0*83.5)) muscle_color = q_i_start + 2.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*86.5)&&(position_cpp[4*j+2]<WZC+r0*86.5)) //MVR03 || MVL03
									if((position_cpp[4*i+2]>WZC+r0*77.5)&&(position_cpp[4*j+2]>WZC+r0*77.5)) muscle_color = q_i_start + 3.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*84.5)&&(position_cpp[4*j+2]<WZC+r0*84.5)) //MVR04 || MVL04
									if((position_cpp[4*i+2]>WZC+r0*76.5)&&(position_cpp[4*j+2]>WZC+r0*76.5)) muscle_color = q_i_start + 4.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*82.5)&&(position_cpp[4*j+2]<WZC+r0*82.5)) 
									if((position_cpp[4*i+2]>WZC+r0*72.5)&&(position_cpp[4*j+2]>WZC+r0*72.5)) muscle_color = q_i_start + 4.5f;//x.5 = violet
								
									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*78.5)&&(position_cpp[4*j+2]<WZC+r0*78.5)) //MVR05 || MVL05
									if((position_cpp[4*i+2]>WZC+r0*66.9)&&(position_cpp[4*j+2]>WZC+r0*66.9)) muscle_color = q_i_start + 5.2f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*77.5)&&(position_cpp[4*j+2]<WZC+r0*77.5)) 
									if((position_cpp[4*i+2]>WZC+r0*65.9)&&(position_cpp[4*j+2]>WZC+r0*65.9)) muscle_color = q_i_start + 5.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*74.0)&&(position_cpp[4*j+2]<WZC+r0*74.0)) //MVR06 || MVL06
									if((position_cpp[4*i+2]>WZC+r0*55.0)&&(position_cpp[4*j+2]>WZC+r0*55.0)) muscle_color = q_i_start + 6.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*74.0)&&(position_cpp[4*j+2]<WZC+r0*74.0)) 
									if((position_cpp[4*i+2]>WZC+r0*54.5)&&(position_cpp[4*j+2]>WZC+r0*54.5)) muscle_color = q_i_start + 6.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*68.5)&&(position_cpp[4*j+2]<WZC+r0*68.5)) //MVR07 || MVL07
									if((position_cpp[4*i+2]>WZC+r0*51.0)&&(position_cpp[4*j+2]>WZC+r0*51.0)) muscle_color = q_i_start + 7.3f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*66.5)&&(position_cpp[4*j+2]<WZC+r0*66.5)) 
									if((position_cpp[4*i+2]>WZC+r0*49.5)&&(position_cpp[4*j+2]>WZC+r0*49.5)) muscle_color = q_i_start + 7.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*56.5)&&(position_cpp[4*j+2]<WZC+r0*56.5)) //MVR08 || MVL08
									if((position_cpp[4*i+2]>WZC+r0*40.0)&&(position_cpp[4*j+2]>WZC+r0*40.0)) muscle_color = q_i_start + 8.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*55.5)&&(position_cpp[4*j+2]<WZC+r0*55.5)) 
									if((position_cpp[4*i+2]>WZC+r0*38.5)&&(position_cpp[4*j+2]>WZC+r0*38.5)) muscle_color = q_i_start + 8.5f;//x.5 = violet

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*52.1)&&(position_cpp[4*j+2]<WZC+r0*52.1)) //MVR09 || MVL09
									if((position_cpp[4*i+2]>WZC+r0*33.5)&&(position_cpp[4*j+2]>WZC+r0*33.5)) muscle_color = q_i_start + 9.2f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*50.5)&&(position_cpp[4*j+2]<WZC+r0*50.5)) 
									if((position_cpp[4*i+2]>WZC+r0*32.5)&&(position_cpp[4*j+2]>WZC+r0*32.5)) muscle_color = q_i_start + 9.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*41.1)&&(position_cpp[4*j+2]<WZC+r0*41.1)) //MVR10 || MVL10
									if((position_cpp[4*i+2]>WZC+r0*22.5)&&(position_cpp[4*j+2]>WZC+r0*22.5)) muscle_color = q_i_start + 10.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*40.0)&&(position_cpp[4*j+2]<WZC+r0*40.0)) 
									if((position_cpp[4*i+2]>WZC+r0*21.5)&&(position_cpp[4*j+2]>WZC+r0*21.5)) muscle_color = q_i_start + 10.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-4*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-5*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*40.0)&&(position_cpp[4*j+2]<WZC+r0*40.0)) 
									if((position_cpp[4*i+2]>WZC+r0*20.5)&&(position_cpp[4*j+2]>WZC+r0*20.5)) muscle_color = q_i_start + 10.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)        &&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*34.5)&&(position_cpp[4*j+2]<WZC+r0*34.5)) //MVR11 || MVL11
									if((position_cpp[4*i+2]>WZC+r0*15.5)&&(position_cpp[4*j+2]>WZC+r0*15.5)) muscle_color = q_i_start + 11.3f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)   &&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*33.5)&&(position_cpp[4*j+2]<WZC+r0*33.5)) 
									if((position_cpp[4*i+2]>WZC+r0*14.5)&&(position_cpp[4*j+2]>WZC+r0*14.5)) muscle_color = q_i_start + 11.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*23.5)&&(position_cpp[4*j+2]<WZC+r0*23.5)) //MVR12 || MVL12
									if((position_cpp[4*i+2]>WZC+r0* 8.5)&&(position_cpp[4*j+2]>WZC+r0* 8.5)) muscle_color = q_i_start + 12.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*22.5)&&(position_cpp[4*j+2]<WZC+r0*22.5)) 
									if((position_cpp[4*i+2]>WZC+r0* 7.5)&&(position_cpp[4*j+2]>WZC+r0* 7.5)) muscle_color = q_i_start + 12.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-4*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-5*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*21.5)&&(position_cpp[4*j+2]<WZC+r0*21.5)) 
									if((position_cpp[4*i+2]>WZC+r0* 6.5)&&(position_cpp[4*j+2]>WZC+r0* 6.5)) muscle_color = q_i_start + 12.5f;//x.5 = violet

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*16.5)&&(position_cpp[4*j+2]<WZC+r0*16.5)) //MVR13 || MVL13
									if((position_cpp[4*i+2]>WZC+r0* 1.5)&&(position_cpp[4*j+2]>WZC+r0* 1.5)) muscle_color = q_i_start + 13.2f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*15.5)&&(position_cpp[4*j+2]<WZC+r0*15.5)) 
									if((position_cpp[4*i+2]>WZC+r0* 0.5)&&(position_cpp[4*j+2]>WZC+r0* 0.5)) muscle_color = q_i_start + 13.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 9.0)&&(position_cpp[4*j+2]<WZC+r0* 9.0)) //MVR14 || MVL14
									if((position_cpp[4*i+2]>WZC-r0* 2.5)&&(position_cpp[4*j+2]>WZC-r0* 2.5)) muscle_color = q_i_start + 14.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 8.5)&&(position_cpp[4*j+2]<WZC+r0* 8.5)) 
									if((position_cpp[4*i+2]>WZC-r0* 3.5)&&(position_cpp[4*j+2]>WZC-r0* 3.5)) muscle_color = q_i_start + 14.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-4*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-5*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 7.5)&&(position_cpp[4*j+2]<WZC+r0* 7.5)) 
									if((position_cpp[4*i+2]>WZC-r0* 4.5)&&(position_cpp[4*j+2]>WZC-r0* 4.5)) muscle_color = q_i_start + 14.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 2.0)&&(position_cpp[4*j+2]<WZC+r0* 2.0)) //MVR15 || MVL15
									if((position_cpp[4*i+2]>WZC-r0*14.5)&&(position_cpp[4*j+2]>WZC-r0*14.5)) muscle_color = q_i_start + 15.3f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)   &&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 1.5)&&(position_cpp[4*j+2]<WZC+r0* 1.5)) 
									if((position_cpp[4*i+2]>WZC-r0*15.5)&&(position_cpp[4*j+2]>WZC-r0*15.5)) muscle_color = q_i_start + 15.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0* 1.5)&&(position_cpp[4*j+2]<WZC-r0* 1.5)) //MVR16 || MVL16
									if((position_cpp[4*i+2]>WZC-r0*21.5)&&(position_cpp[4*j+2]>WZC-r0*21.5)) muscle_color = q_i_start + 16.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0* 2.5)&&(position_cpp[4*j+2]<WZC-r0* 2.5)) 
									if((position_cpp[4*i+2]>WZC-r0*22.5)&&(position_cpp[4*j+2]>WZC-r0*22.5)) muscle_color = q_i_start + 16.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-4*r0) &&(position_cpp[4*i+1]*dq>WYC*dq-5*r0)) 
									if((position_cpp[4*i+2]<WZC-r0* 3.5)&&(position_cpp[4*j+2]<WZC-r0* 3.5)) 
									if((position_cpp[4*i+2]>WZC-r0*23.5)&&(position_cpp[4*j+2]>WZC-r0*23.5)) muscle_color = q_i_start + 16.5f;//x.5 = violet

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*14.0)&&(position_cpp[4*j+2]<WZC-r0*14.0)) //MVR17 || MVL17
									if((position_cpp[4*i+2]>WZC-r0*34.5)&&(position_cpp[4*j+2]>WZC-r0*34.5)) muscle_color = q_i_start + 17.2f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*14.7)&&(position_cpp[4*j+2]<WZC-r0*14.7)) 
									if((position_cpp[4*i+2]>WZC-r0*35.5)&&(position_cpp[4*j+2]>WZC-r0*35.5)) muscle_color = q_i_start + 17.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*20.0)&&(position_cpp[4*j+2]<WZC-r0*20.0)) //MVR18 || MVL18
									if((position_cpp[4*i+2]>WZC-r0*40.5)&&(position_cpp[4*j+2]>WZC-r0*40.5)) muscle_color = q_i_start + 18.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*21.5)&&(position_cpp[4*j+2]<WZC-r0*21.5)) 
									if((position_cpp[4*i+2]>WZC-r0*41.5)&&(position_cpp[4*j+2]>WZC-r0*41.5)) muscle_color = q_i_start + 18.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-4*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-5*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*22.5)&&(position_cpp[4*j+2]<WZC-r0*22.5)) 
									if((position_cpp[4*i+2]>WZC-r0*34.5)&&(position_cpp[4*j+2]>WZC-r0*34.5)) muscle_color = q_i_start + 18.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*34.0)&&(position_cpp[4*j+2]<WZC-r0*34.0)) //MVR19 || MVL19
									if((position_cpp[4*i+2]>WZC-r0*54.5)&&(position_cpp[4*j+2]>WZC-r0*54.5)) muscle_color = q_i_start + 19.3f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*34.5)&&(position_cpp[4*j+2]<WZC-r0*34.5)) 
									if((position_cpp[4*i+2]>WZC-r0*55.5)&&(position_cpp[4*j+2]>WZC-r0*55.5)) muscle_color = q_i_start + 19.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*39.5)&&(position_cpp[4*j+2]<WZC-r0*39.5)) //MVR20 || MVL20
									if((position_cpp[4*i+2]>WZC-r0*50.5)&&(position_cpp[4*j+2]>WZC-r0*50.5)) muscle_color = q_i_start + 20.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*40.5)&&(position_cpp[4*j+2]<WZC-r0*40.5)) 
									if((position_cpp[4*i+2]>WZC-r0*51.5)&&(position_cpp[4*j+2]>WZC-r0*51.5)) muscle_color = q_i_start + 20.5f;//x.5 = violet

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*53.0)&&(position_cpp[4*j+2]<WZC-r0*53.0)) //MVR21 || MVL21
									if((position_cpp[4*i+2]>WZC-r0*71.5)&&(position_cpp[4*j+2]>WZC-r0*71.5)) muscle_color = q_i_start + 21.2f;
									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*54.0)&&(position_cpp[4*j+2]<WZC-r0*54.0)) 
									if((position_cpp[4*i+2]>WZC-r0*72.5)&&(position_cpp[4*j+2]>WZC-r0*72.5)) muscle_color = q_i_start + 21.2f;//x.2 = red

									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*50.0)&&(position_cpp[4*j+2]<WZC-r0*50.0)) //MVR22 || MVL22
									if((position_cpp[4*i+2]>WZC-r0*63.5)&&(position_cpp[4*j+2]>WZC-r0*63.5)) muscle_color = q_i_start + 22.4f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*50.5)&&(position_cpp[4*j+2]<WZC-r0*50.5)) 
									if((position_cpp[4*i+2]>WZC-r0*64.5)&&(position_cpp[4*j+2]>WZC-r0*64.5)) muscle_color = q_i_start + 22.4f;//x.4 = magenta

									if((position_cpp[4*i+1]*dq<WYC*dq)&&(position_cpp[4*i+1]*dq>WYC*dq-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*70.0)&&(position_cpp[4*j+2]<WZC-r0*70.0)) //MVR23 || MVL23
									if((position_cpp[4*i+2]>WZC-r0*92.0)&&(position_cpp[4*j+2]>WZC-r0*92.0)) muscle_color = q_i_start + 23.3f;//x.3 = orange

									if((position_cpp[4*i+1]*dq<WYC*dq-1*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*71.5)&&(position_cpp[4*j+2]<WZC-r0*71.5)) //MVR24 || MVL24
									if((position_cpp[4*i+2]>WZC-r0*92.0)&&(position_cpp[4*j+2]>WZC-r0*92.0)) muscle_color = q_i_start + 24.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-2*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*62.5)&&(position_cpp[4*j+2]<WZC-r0*62.5)) 
									if((position_cpp[4*i+2]>WZC-r0*82.5)&&(position_cpp[4*j+2]>WZC-r0*82.5)) muscle_color = q_i_start + 24.5f;
									if((position_cpp[4*i+1]*dq<WYC*dq-3*r0)&&(position_cpp[4*i+1]*dq>WYC*dq-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*63.5)&&(position_cpp[4*j+2]<WZC-r0*63.5)) 
									if((position_cpp[4*i+2]>WZC-r0*66.5)&&(position_cpp[4*j+2]>WZC-r0*66.5)) muscle_color = q_i_start + 24.5f;//x.5 = violet								
								}

								elasticConnectionsData_cpp[ m_index[muscleCounter] = 4 * ( MAX_NEIGHBOR_COUNT * i + ecc) + 2 ] = muscle_color;// type of connection; 0 - ordinary spring, 1 or more - muscle
								muscleCounter++;
		
							}
						}
						ecc++;
						ecc_total++;
					}
				}
			}
		}

		
		/*
		if(muscleCounter==10408)
		{
			//for(i=0;i<muscleCounter;i++)

			//m_number[m_index[0]] = 1.2f;
			elasticConnectionsData_cpp[m_index[3]] = 1.2f;

		}*/


		//membraneData - the list containing triplets of indexes of particles forming triangular membranes; size: numOfMembranes*3
		//particleMembranesList - the list containing for each particle(involved in membranes) its list of these membranes indexes (which are stored in membraneData array)

		for(int _mc = 0; _mc < numOfMembranes*3; _mc++)
		{
			int particle_index;
			for(int sli=0/*sublist index*/;sli<MAX_MEMBRANES_INCLUDING_SAME_PARTICLE; sli++)
			{//search for the not filled yet cell (-1) in the list and fill it; 
				if(/**/particleMembranesList_cpp [particle_index = membraneData_cpp [_mc]*MAX_MEMBRANES_INCLUDING_SAME_PARTICLE+sli]/**/==-1)
				{
					particleMembranesList_cpp [particle_index] = _mc/3;
					break;//if there are no free cells, break, because we are limited with MAX_MEMBRANES_INCLUDING_SAME_PARTICLE per particle
				} 
				else
				{
					//_mc = _mc; //https://github.com/openworm/OpenWorm/issues/152
				}
			}		
		}
	}



	return;
}


//READ DEFAULT CONFIGURATATION FROM FILE IN CONFIGURATION FOLDER
int read_position = 0;
std::string owHelper::path = "./configuration/";
std::string owHelper::configFileName = "demo1";

/** TODO make description
 *
 */
void findValf(std::string & str, char delimiter, size_t & start, float & val){
	size_t end = str.find(delimiter, start + 1);
	if(end != std::string::npos){
		val =std::atof(str.substr(start, end - start).c_str());//TODO make check if this a number
		start = end;
	}else // it means usualy that we reached the end of string and there are no \t
		val = std::atof(str.substr(start).c_str());//TODO make check if this a number
}
/** TODO Documentation
 */
void findVali(std::string & str, char delimiter, size_t & start, int & val){
	size_t end = str.find(delimiter, start + 1);
	if(end != std::string::npos){
		val =std::atoi(str.substr(start, end - start).c_str());
		start = end;
	}else // it means usualy that we reached the end of string and there are no \t
		val = std::atoi(str.substr(start).c_str());
}
/** Preparing initial data before load full configuration
 *
 *  Before load configuration data from file (initial position and velocity,
 *  connection data if is and membrane data if is) Sibernetic
 *  should allocate a memory in RAM. Method starts with reading position file
 *  first 6 lines in file correspond to dimensions of boundary box, than it reads all file
 *  till the end an calculate numbers of elastic, liquid and boundary
 *  particles and total number too. Also this read membranes file
 *  for counting membranes numbers.
 *
 *  @param numOfMembranes
 *  reference to numOfMembrane variable
 *  @param config
 *  pointer to owConfigProrerty object it includes information about
 *  boundary box dimensions
 *  @param numOfLiquidP
 *  reference to numOfLiquidP variable
 *  @param numOfElasticP
 *  reference to numOfElasticP variable
 *  @param numOfBoundaryP
 *  reference to numOfBoundaryP variable
 */
void owHelper::preLoadConfiguration(int & numOfMembranes, owConfigProrerty * config, int & numOfLiquidP, int & numOfElasticP, int & numOfBoundaryP)
{
	try
	{
		int p_count = 0;
		std::string file_name = path + config->getCofigFileName();
		std::string inputStr;
		std::ifstream configFile (file_name.c_str(), std::ios_base::binary);
		float x, y, z, p_type;
		char delimiter = '\t';
		size_t pos;
		ELOADMODE mode = NOMODE;
		if( configFile.is_open() )
		{
			configFile >> config->xmin;
			configFile >> config->xmax;
			configFile >> config->ymin;
			configFile >> config->ymax;
			configFile >> config->zmin;
			configFile >> config->zmax;
			read_position = configFile.tellg();
			while( configFile.good() )
			{
				//configFile >> inputStr;
				std::getline(configFile,inputStr);
				if(inputStr == "[position]"){
					mode = POSITION;
					continue;
				}
				if(inputStr == "[velocity]"){
					mode = VELOCITY;
					continue;
				}
				if(inputStr == "[membranes]"){
					mode = MEMBRANE;
					continue;
				}
				if(inputStr == "[particleMemIndex]"){
					break;
				}
				switch(mode){
					case POSITION:{
						p_type = -1.1f;//reinitialize
						pos = 0;
						findValf(inputStr, delimiter,pos,x);
						findValf(inputStr, delimiter,pos,y);
						findValf(inputStr, delimiter,pos,z);
						findValf(inputStr, delimiter,pos,p_type);
						p_count++;
						switch((int)p_type){
							case LIQUID_PARTICLE:
								numOfLiquidP++;
								break;
							case ELASTIC_PARTICLE:
								numOfElasticP++;
								break;
							case BOUNDARY_PARTICLE:
								numOfBoundaryP++;
								break;
						}
						break;
					}
					case MEMBRANE:{
						numOfMembranes++;
						break;
					}
					default:
						continue;
				}
			}
		}else
			throw std::runtime_error("Could not open file configuration file");
		configFile.close();
		config->setParticleCount(p_count);
	}
	catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}
/** Load full configuration
 *
 *  Load configuration data from file (initial position and velocity,
 *  connection data if is and membrane data if is). Method starts with
 *  reading position file than velocity elastic connection and membranes data files
 *
 *  @param position_cpp
 *  pointer to position_cpp buffer
 *  @param velocity_cpp
 *  pointer to velocity_cpp buffer
 *  @param elasticConnections
 *  reference on pointer to elasticConnections buffer.
 *  In this function we allocate memory for elasticConnections.
 *  TODO: change it replace to owPhysicsFluidSimulator constructor.
 *  @param numOfLiquidP
 *  reference to numOfLiquidP variable
 *  @param numOfElasticP
 *  reference to numOfElasticP variable
 *  @param numOfBoundaryP
 *  reference to numOfBoundaryP variable
 *  @param numOfElasticConnections
 *  reference to numOfElasticConnections variable
 *  @param numOfMembranes
 *  reference to numOfMembranes variable
 *  @param membraneData_cpp
 *  pointer to membraneData_cpp buffer
 *  @param particleMembranesList_cpp
 *  pointer to particleMembranesList_cpp buffer
 *  @param config
 *  pointer to owConfigProrerty object it includes information about
 */
void owHelper::loadConfiguration(float *position_cpp, float *velocity_cpp, float *& elasticConnections,int & numOfLiquidP, int & numOfElasticP, int & numOfBoundaryP, int & numOfElasticConnections, int & numOfMembranes,int * membraneData_cpp, int *& particleMembranesList_cpp, owConfigProrerty * config)
{
	try
	{
		std::string file_name = path + config->getCofigFileName();
		std::string inputStr;
		std::ifstream configFile (file_name.c_str(), std::ios_base::binary);
		char delimiter = '\t';
		size_t pos;
		ELOADMODE mode = NOMODE;
		int i = 0;
		if( configFile.is_open() )
		{
			configFile.seekg(read_position);
			while( configFile.good() )
			{
				//configFile >> inputStr;
				std::getline(configFile,inputStr);
				if(inputStr == "[position]"){
					mode = POSITION;
					continue;
				}
				if(inputStr == "[velocity]"){
					i = 0;
					mode = VELOCITY;
					continue;
				}
				if(inputStr == "[connection]"){
					i = 0;
					mode = CONNECTION;
					continue;
				}
				if(inputStr == "[membranes]"){
					i = 0;
					mode = MEMBRANE;
					continue;
				}
				if(inputStr == "[particleMemIndex]"){
					i = 0;
					mode = PMEMINDEX;
					continue;
				}
				if(inputStr == "[end]")
					break;
				switch(mode){
					case POSITION:{
						float x, y, z, p_type;
						p_type = -1.1f;//reinitialize
						pos = 0;
						findValf(inputStr, delimiter,pos,x);
						findValf(inputStr, delimiter,pos,y);
						findValf(inputStr, delimiter,pos,z);
						findValf(inputStr, delimiter,pos,p_type);
						position_cpp[ 4 * i + 0 ] = x;
						position_cpp[ 4 * i + 1 ] = y;
						position_cpp[ 4 * i + 2 ] = z;
						position_cpp[ 4 * i + 3 ] = p_type;
						i++;
						break;
					}
					case VELOCITY:{
						float x, y, z, p_type;
						p_type = -1.1f;//reinitialize
						pos = 0;
						findValf(inputStr, delimiter,pos,x);
						findValf(inputStr, delimiter,pos,y);
						findValf(inputStr, delimiter,pos,z);
						findValf(inputStr, delimiter,pos,p_type);
						velocity_cpp[ 4 * i + 0 ] = x;
						velocity_cpp[ 4 * i + 1 ] = y;
						velocity_cpp[ 4 * i + 2 ] = z;
						velocity_cpp[ 4 * i + 3 ] = p_type;
						i++;
						break;
					}
					case CONNECTION:{
						pos = 0;
						float  jd, rij0, val1, val2;
						findValf(inputStr, delimiter,pos,jd);
						findValf(inputStr, delimiter,pos,rij0);
						findValf(inputStr, delimiter,pos,val1);
						findValf(inputStr, delimiter,pos,val2);
						elasticConnections[ 4 * i + 0 ] = jd;
						elasticConnections[ 4 * i + 1 ] = rij0 * simulationScale;
						elasticConnections[ 4 * i + 2 ] = val1;
						elasticConnections[ 4 * i + 3 ] = val2;
						i++;
						break;
					}
					case MEMBRANE:{
						pos = 0;
						int id, jd, kd;
						findVali(inputStr, delimiter,pos,id);
						findVali(inputStr, delimiter,pos,jd);
						findVali(inputStr, delimiter,pos,kd);
						membraneData_cpp[ 3 * i + 0 ] = id;
						membraneData_cpp[ 3 * i + 1 ] = jd;
						membraneData_cpp[ 3 * i + 2 ] = kd;
						i++;
						break;
					}
					case PMEMINDEX:{
						int id;
						pos = 0;
						findVali(inputStr, delimiter,pos,id);
						particleMembranesList_cpp[ i ] = id;
						i++;
						break;
					}
					default:
						continue;
				}
			}
		}else
			throw std::runtime_error("Could not open file configuration file");
		configFile.close();
		std::cout << "Configuration was loaded" << std::endl;
	}catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}
template<typename T>
std::ostream& binary_write(std::ostream& stream, const T& value){
    return stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
}
/** Load configuration from simulation to files
 *
 *  This method is required for work with "load config to file" mode.
 *  In this mode information about simulation's evolution is storing into a file
 *  on every step (every time it reads data block with size = c).
 *  If Sibernetic runs in this mode it means that
 *  no calculation on OpenCL device runs.
 *
 *  @param position
 *  pointer to position buffer
 *  @param config
 *  pointer to owConfigProrerty object it includes information about
 *  @param connections
 *  reference on pointer to elasticConnections buffer.
 *  @param membranes
 *  pointer to membranes buffer
 *  @param firstIteration
 *  if true it means that we first time record information
 *  to a file and on first iteration it put to
 *  the file info about dimensions of boundary box
 *  NOTE: next 2 parameters are an experimental
 *  @param filter_p
 *  pointer to filter particle buffer, if you need storing only
 *  a bunch of particles not all of them
 *  @param size
 *  size of filter_p array
 */
void owHelper::loadConfigurationToFile(float * position, owConfigProrerty * config, std::vector<int> & filter, float * connections, int * membranes, bool firstIteration){
	try{
		ofstream positionFile;
		if(firstIteration){
#if !EXPEREMENTAL_WRITE
			positionFile.open("./buffers/position_buffer.txt", ios::trunc);
			positionFile << config->xmin << "\n";
			positionFile << config->xmax << "\n";
			positionFile << config->ymin << "\n";
			positionFile << config->ymax << "\n";
			positionFile << config->zmin << "\n";
			positionFile << config->zmax << "\n";
#else
			positionFile.open("./buffers/position_buffer.txt", ios::trunc|ios::binary);
			binary_write(positionFile,config->xmin);
			binary_write(positionFile,config->xmax);
			binary_write(positionFile,config->ymin);
			binary_write(positionFile,config->ymax);
			binary_write(positionFile,config->zmin);
			binary_write(positionFile,config->zmax);
			binary_write(positionFile,40000.0f);
#endif
			if(!filter.empty()){
#if !EXPEREMENTAL_WRITE
				positionFile << filter.size() << "\n";
				positionFile << 0 << "\n";
#else
				binary_write(positionFile,(float)filter.size());
				binary_write(positionFile,0.0f);
#endif
			}
			else{
#if !EXPEREMENTAL_WRITE
				positionFile << numOfElasticP << "\n";
				positionFile << numOfLiquidP << "\n";
#else
				binary_write(positionFile,numOfElasticP);
				binary_write(positionFile,numOfLiquidP);
#endif
			}
		}else{
#if !EXPEREMENTAL_WRITE
			positionFile.open("./buffers/position_buffer.txt", ios::app);
#else
			positionFile.open("./buffers/position_buffer.txt", ios::app|ios::binary);
#endif
		}
		if(filter.empty()){
			for(int i=0;i < config->getParticleCount(); i++){
				if((int)position[ 4 * i + 3] != BOUNDARY_PARTICLE){
#if !EXPEREMENTAL_WRITE
				positionFile << position[i * 4 + 0] << "\t" << position[i * 4 + 1] << "\t" << position[i * 4 + 2] << "\t" << position[i * 4 + 3] << "\n";
#else
				binary_write(positionFile,position[i * 4 + 0]);
				binary_write(positionFile,position[i * 4 + 1]);
				binary_write(positionFile,position[i * 4 + 2]);
				binary_write(positionFile,position[i * 4 + 3]);
#endif

				}
			}
		}else{
			int i = 0;
			for(unsigned int index = 0; index<filter.size(); index++){
				i = filter[index];
#if !EXPEREMENTAL_WRITE
				positionFile << position[i * 4 + 0] << "\t" << position[i * 4 + 1] << "\t" << position[i * 4 + 2] << "\t" << position[i * 4 + 3] << "\n";
#else
				binary_write(positionFile,position[i * 4 + 0]);
				binary_write(positionFile,position[i * 4 + 1]);
				binary_write(positionFile,position[i * 4 + 2]);
				binary_write(positionFile,position[i * 4 + 3]);
#endif
			}
		}
		positionFile.close();
		if(firstIteration){
			ofstream connectionFile("./buffers/connection_buffer.txt", std::ofstream::trunc);
			int con_num = MAX_NEIGHBOR_COUNT * numOfElasticP;
			for(int i = 0; i < con_num; i++)
				connectionFile << connections[4 * i + 0] << "\t" << connections[4 * i + 1] << "\t" << connections[4 * i + 2] << "\t" << connections[4 * i + 3] << "\n";
			connectionFile.close();
			ofstream membranesFile("./buffers/membranes_buffer.txt", std::ofstream::trunc);
			membranesFile << numOfMembranes << "\n";
			for(int i = 0; i < numOfMembranes; i++)
				membranesFile << membranes[3 * i + 0] << "\t" << membranes[3 * i + 1] << "\t" << membranes[3 * i + 2] << "\n";
			membranesFile.close();
		}
	}catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}
/** Load configuration from simulation to files
 *
 *  TODO make description
 *
 *
 *
 *
 *
 *  @param position
 *  pointer to position buffer
 *  @param config
 *  pointer to owConfigProrerty object it includes information about
 *  @param connections
 *  reference on pointer to elasticConnections buffer.
 *  @param membranes
 *  pointer to membranes buffer
 *  @param firstIteration
 *  if true it means that we first time record information
 *  to a file and on first iteration it put to
 *  the file info about dimensions of boundary box
 *  NOTE: next 2 parameters are an experimental
 *  @param filter_p
 *  pointer to filter particle buffer, if you need storing only
 *  a bunch of particles not all of them
 *  @param size
 *  size of filter_p array
 */
void owHelper::loadConfigurationToFile(float * position, float * velocity, float * connections, int * membranes, int * particleMemIndex,const char * filename, owConfigProrerty * config){
	try{
		ofstream configFile;
		configFile.open(filename, std::ofstream::trunc);
		configFile << config->xmin << "\n";
		configFile << config->xmax << "\n";
		configFile << config->ymin << "\n";
		configFile << config->ymax << "\n";
		configFile << config->zmin << "\n";
		configFile << config->zmax << "\n";
		configFile << "[position]\n" ;
		for(int i=0;i < config->getParticleCount(); i++)
			configFile << position[i * 4 + 0] << "\t" << position[i * 4 + 1] << "\t" << position[i * 4 + 2] << "\t" << position[i * 4 + 3] << "\n";
		configFile << "[velocity]\n" ;
		for(int i=0;i < config->getParticleCount(); i++)
			configFile << velocity[i * 4 + 0] << "\t" << velocity[i * 4 + 1] << "\t" << velocity[i * 4 + 2] << "\t" << velocity[i * 4 + 3] << "\n";
		configFile << "[connection]\n" ;
		int con_num = MAX_NEIGHBOR_COUNT * numOfElasticP;
		for(int i = 0; i < con_num; i++)
			configFile << connections[4 * i + 0] << "\t" << connections[4 * i + 1] / simulationScale << "\t" << connections[4 * i + 2] << "\t" << connections[4 * i + 3] << "\n";
		configFile << "[membranes]\n";
		for(int i = 0; i < numOfMembranes; i++)
			configFile << membranes[3 * i + 0] << "\t" << membranes[3 * i + 1] << "\t" << membranes[3 * i + 2] << "\n";
		configFile << "[particleMemIndex]\n";
		int particleMemIndexCount = numOfElasticP*MAX_MEMBRANES_INCLUDING_SAME_PARTICLE;
		for(int i = 0; i < particleMemIndexCount; i++)
			configFile << particleMemIndex[i] << "\n";
		configFile << "[end]";
		configFile.close();
	}catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}
//This function needed for visualiazation buffered data
long position_index = 0;
ifstream positionFile;
/** Load configuration from file to simulation
 *
 *  This method is required for work with "load config from file" mode.
 *  In this mode information about simulation's evolution is taking from file
 *  on every step (every time it reads data block with size = PARTICLE_COUNT).
 *  If Sibernetic runs in this mode it means that
 *  no calculation on OpenCL device runs.
 *
 *  @param position
 *  pointer to position buffer
 *  @param connections
 *  reference on pointer to elasticConnections buffer.
 *  @param membranes
 *  pointer to membranes buffer
 *  @param config
 *  pointer to owConfigProrerty object it includes information about
 *  @param iteration
 *  if iteration==0 it means that we first time record information
 *  to a file and on first iteration it put to
 *  the file info about dimensions of boundary box
 */
void owHelper::loadConfigurationFromFile(float *& position, float *& connections, int *& membranes, owConfigProrerty * config, int iteration){
	try{
		if(iteration == 0)
			positionFile.open("./buffers/position_buffer.txt");
		int i = 0;
		float x, y, z, p_type;
		if( positionFile.is_open() )
		{
			if(iteration == 0){
				positionFile >> config->xmin;
				positionFile >> config->xmax;
				positionFile >> config->ymin;
				positionFile >> config->ymax;
				positionFile >> config->zmin;
				positionFile >> config->zmax;
				positionFile >> numOfElasticP;
				positionFile >> numOfLiquidP;
				config->setParticleCount(numOfElasticP + numOfLiquidP);
				position = new float[4 * config->getParticleCount()];
			}
			while( positionFile.good() &&  i < config->getParticleCount())
			{
				positionFile >> x >> y >> z >> p_type;
				position[i * 4 + 0] = x;
				position[i * 4 + 1] = y;
				position[i * 4 + 2] = z;
				position[i * 4 + 3] = p_type;
				i++;
			}
		}
		if(!positionFile.good()){
			positionFile.close();
			exit(0);
		}
		if(iteration == 0){

			ifstream connectionFile("./buffers/connection_buffer.txt");
			connections = new float[MAX_NEIGHBOR_COUNT * numOfElasticP * 4];
			if( connectionFile.is_open() )
			{
				int i = 0;
				float jd, rij0, val1, val2;
				while(connectionFile.good() && i < MAX_NEIGHBOR_COUNT * numOfElasticP){
					connectionFile >> jd >> rij0 >> val1 >> val2;
					connections[ 4 * i + 0 ] = jd;
					connections[ 4 * i + 1 ] = rij0;
					connections[ 4 * i + 2 ] = val1;
					connections[ 4 * i + 3 ] = val2;
					i++;
				}
			}
			connectionFile.close();
			ifstream membranesFile("./buffers/membranes_buffer.txt");
			if(membranesFile.is_open()){
				int m_count = 0;
				//membranesFile >> m_count;
				int i = 0;
				membranes = new int[4 * m_count];
				while(membranesFile.good() && i < m_count){
					membranesFile >> membranes[4 * i + 0] >> membranes[4 * i + 1] >> membranes[4 * i + 2] >> membranes[4 * i + 3];
					i++;
				}
			}
			membranesFile.close();
		}
	}catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}
/** Print value of elapsed time from last handling to watch_report method.
 *
 *  This function is required for logging time consumption info.
 *
 *  @param str
 *  represents output string format.
 */
void owHelper::watch_report( const char * str )
{
#if defined(_WIN32) || defined(_WIN64)
	QueryPerformanceCounter(&t2);
	printf(str,(t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart);
	t1 = t2;
	elapsedTime = (t2.QuadPart - t0.QuadPart) * 1000.0 / frequency.QuadPart;
#elif defined(__linux__)
	clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
	time_t sec = t2.tv_sec - t1.tv_sec;
	long nsec;
	if (t2.tv_nsec >= t1.tv_nsec) {
			nsec = t2.tv_nsec - t1.tv_nsec;
	} else {
			nsec = 1000000000 - (t1.tv_nsec - t2.tv_nsec);
			sec -= 1;
	}
	printf(str,(float)sec * 1000.f + (float)nsec/1000000.f);
	t1 = t2;
	elapsedTime =  (float)(t2.tv_sec - t0.tv_sec) * 1000.f + (float)(t2.tv_nsec - t0.tv_nsec)/1000000.f;
#elif defined(__APPLE__)
    uint64_t elapsedNano;
    static mach_timebase_info_data_t    sTimebaseInfo;

    if ( sTimebaseInfo.denom == 0 ) {
        (void) mach_timebase_info(&sTimebaseInfo);
    }

    t2 = mach_absolute_time();
    elapsedNano = (t2-t1) * sTimebaseInfo.numer / sTimebaseInfo.denom;
    printf(str, (float)elapsedNano/1000000.f );
    t1=t2;
    elapsedNano = (t2-t0) * sTimebaseInfo.numer / sTimebaseInfo.denom;
    elapsedTime = (float)elapsedNano/1000000.f;
#endif
}
