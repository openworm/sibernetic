#include "owHelper.h"
#include "owPhysicsConstant.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <sstream>
#include <string>
using namespace std;

extern int PARTICLE_COUNT;
extern int PARTICLE_COUNT_RoundedUp;
extern int local_NDRange_size;

owHelper::owHelper(void)
{
	refreshTime();
}

owHelper::~owHelper(void)
{
}
void owHelper::refreshTime()
{
#if defined(_WIN32) || defined (_WIN64)
    QueryPerformanceFrequency(&frequency);	// get ticks per second
    QueryPerformanceCounter(&t1);			// start timer
	t0 = t1;
#elif defined(__linux__)
	clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
	t0 = t1;
#endif
}
//For output float buffer
//Create File in which line element_size elements forn buffer
//global_size - size of buffer / element_size
void owHelper::log_bufferf(const float * buffer, const int element_size, const int global_size, const char * fileName)
{
	try{
		ofstream outFile (fileName);
		for(int i = 0; i < global_size; i++)
		{
			for(int j = 0; j < element_size; j++)
			{
				if(j < element_size - 1 )
					outFile << buffer[ i * element_size + j ] << "\t";
				else
					outFile << buffer[ i * element_size + j ] << "\n";
			}
		}
		outFile.close();
	}catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}
//For output int buffer
void owHelper::log_bufferi(const int * buffer, const int element_size, const int global_size, const char * fileName)
{
	try{
		ofstream outFile (fileName);
		for(int i = 0; i < global_size; i++)
		{
			for(int j = 0; j < element_size; j++)
			{
				if(j < element_size + 1 )
					outFile << buffer[ i * element_size + j ] << "\t";
				else
					outFile << buffer[ i * element_size + j ] << "\n";
			}
		}
		outFile.close();
	}catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}

int generateWormShell(int stage, int i_start,float *position_cpp, float *velocity_cpp, int &numOfMembranes, int *membraneData_cpp)
{
	//return 0;
	if( !((stage==0)||(stage==1)) ) return 0;

	float alpha;
	float wormBodyRadius;
	int pCount = 0;//total counter of particles being created within this function
	int i,j;
	float *positionVector;
	float *velocityVector;
	float xc = XMAX*0.5f;
	float yc = YMAX*0.3f;
	float zc = ZMAX*0.5f;
	int elasticLayers = 1;//starting value
	float PI = 3.1415926536f;
	int currSlice_pCount;
	int prevSlice_pCount;
	int currSlice_start;
	int prevSlice_start;
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

		if((wormBodyRadius>0.707*r0)&&
		   (wormBodyRadius<1.000*r0)) wormBodyRadius = 1.000*r0;

		if(wormBodyRadius<0.707*r0) { tip = 1; wormBodyRadius = 0.707f*r0; }//0.707 = sqrt(2)/2

		//alpha = 2*asin(0.5*r0/wormBodyRadius);//in radians
		//angle = alpha;

		elasticLayers = 1;

		if(stage==1)
		{
			positionVector = position_cpp + 4 * (pCount+i_start);
			positionVector[ 0 ] = xc + wormBodyRadius*cos(0.0);
			positionVector[ 1 ] = yc + wormBodyRadius*sin(0.0);
			positionVector[ 2 ] = zc + r0*j;
			positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow

			positionVector = position_cpp + 4 * (pCount+1+i_start);
			positionVector[ 0 ] = xc - wormBodyRadius*cos(0.0);
			positionVector[ 1 ] = yc - wormBodyRadius*sin(0.0);
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
				positionVector[ 0 ] = xc + wormBodyRadius*sin(0.0);
				positionVector[ 1 ] = yc + wormBodyRadius*cos(0.0);
				positionVector[ 2 ] = zc + r0*j;
				positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow

				positionVector = position_cpp + 4 * (pCount+1+i_start);
				positionVector[ 0 ] = xc - wormBodyRadius*sin(0.0);
				positionVector[ 1 ] = yc - wormBodyRadius*cos(0.0);
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
				if(wormBodyRadius>r0*(1.00))
				{
					if(stage==1)
					{
						positionVector = position_cpp + 4 * (pCount+i_start);
						positionVector[ 0 ] = xc + wormBodyRadius*cos(0.0);
						positionVector[ 1 ] = yc + wormBodyRadius*sin(0.0);
						positionVector[ 2 ] = zc + r0*j;
						positionVector[ 3 ] = 2.1f;// 2 = elastic matter, yellow

						positionVector = position_cpp + 4 * (pCount+1+i_start);
						positionVector[ 0 ] = xc - wormBodyRadius*cos(0.0);
						positionVector[ 1 ] = yc - wormBodyRadius*sin(0.0);
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

			if(wormBodyRadius<r0*0.707) break;
			alpha = 2*asin(0.5*r0/wormBodyRadius);//in radians//recalculate -- wormBodyRadius changed
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
				int n_non_muscle_particles = floor(non_muscle_angle / alpha)-1;// distance between each 2 radially adjacent particles will be r0 or more (not less); alpha corresponds to r0
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

int generateInnerWormLiquid(int stage, int i_start,float *position_cpp, float *velocity_cpp)
{
	//return 0;
	int segmentsCount;// = 10;
	float alpha;// = 2.f*3.14159f/segmentsCount;
	float coeff = 0.23f;
	float wormBodyRadius;// = h*coeff / sin(alpha/2);
	int pCount = 0;//particle counter
	int i;

	float *positionVector;
	float *velocityVector;
	float value;
	int elasticLayers;//starting from 2, because 1 is for outer shell and doesn't contain liquid particles
	float xc = XMAX*0.5f;
	float yc = YMAX*0.3f;
	float zc = ZMAX*0.5f;
	float PI = 3.1415926536f;
	float beta;
	float angle;

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
		wormBodyRadius = 6.0f*r0*sqrt(max(1.f-(1.0e-4f)*j*j,0.f)) - r0*(1+0.85);

		while(1)
		{
			if(wormBodyRadius>0.707*r0)
			{
				if(stage==1)
				{
					positionVector = position_cpp + 4 * (pCount+i_start);
					positionVector[ 0 ] = xc + wormBodyRadius*sin(0.0);
					positionVector[ 1 ] = yc + wormBodyRadius*cos(0.0);
					positionVector[ 2 ] = zc + r0*j;
					positionVector[ 3 ] = 1.1f;// liquid

					positionVector = position_cpp + 4 * (pCount+1+i_start);
					positionVector[ 0 ] = xc - wormBodyRadius*sin(0.0);
					positionVector[ 1 ] = yc - wormBodyRadius*cos(0.0);
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

			alpha = 2*asin(0.5*r0/wormBodyRadius);//in radians//recalculate -- wormBodyRadius changed
			
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
			int n_non_muscle_particles = floor(non_muscle_angle / (alpha*0.85) )-1;
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
			wormBodyRadius -= r0*0.85;
		}

	}

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


void owHelper::generateConfiguration(int stage, float *position_cpp, float *velocity_cpp, float *& elasticConnectionsData_cpp, int *membraneData_cpp, int & numOfLiquidP, int & numOfElasticP, int & numOfBoundaryP, int & numOfElasticConnections, int & numOfMembranes, int *particleMembranesList_cpp)
{
	float x,y,z;
	float p_type = LIQUID_PARTICLE;
	int i = 0;// particle counter
	int ix,iy,iz;
	int ecc = 0;//elastic connections counter

	int nx = (int)( ( XMAX - XMIN ) / r0 ); //X
	int ny = (int)( ( YMAX - YMIN ) / r0 ); //Y
	int nz = (int)( ( ZMAX - ZMIN ) / r0 ); //Z

	int nEx = 5*0;//7
	int nEy = 3*0;//4
	int nEz = 9*0;//25
	int nMuscles = 5;
	int nM,nMi,nMj;
	int wormIndex_start,wormIndex_end;
	int numOfMembraneParticles = generateWormShell(0,0,position_cpp,velocity_cpp, numOfMembranes, membraneData_cpp);

	if(stage==0)
	{
		numOfLiquidP = generateInnerWormLiquid(0,0,position_cpp,velocity_cpp);
		numOfElasticP = numOfMembraneParticles;
		numOfBoundaryP = 0;

		if(numOfElasticP<=0) elasticConnectionsData_cpp = NULL; else elasticConnectionsData_cpp = new float[ 4 * numOfElasticP * MAX_NEIGHBOR_COUNT ];
	}

	//=============== create worm body (elastic parts) ==================================================
	if(stage==1)
	{
		wormIndex_start = i;
		i += generateWormShell(1/*stage*/,i,position_cpp,velocity_cpp, numOfMembranes,membraneData_cpp);
		wormIndex_end = i;

		float r2ij;
		float dx2,dy2,dz2;

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
		i += generateInnerWormLiquid(stage,i,position_cpp,velocity_cpp);
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
		PARTICLE_COUNT = numOfLiquidP + numOfBoundaryP + numOfElasticP;
		PARTICLE_COUNT_RoundedUp = ((( PARTICLE_COUNT - 1 ) / local_NDRange_size ) + 1 ) * local_NDRange_size;

		if(PARTICLE_COUNT<=0) 
		{
			printf("\nWarning! Generated scene contains %d particles!\n",PARTICLE_COUNT);
			exit(-2);
		}
	}
	else
	if(stage==1)
	{
		if(PARTICLE_COUNT!=i) 
		{
			printf("\nWarning! Preliminary [%d] and final [%d] particle count are different\n",PARTICLE_COUNT,i);
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
		int array_j[MAX_NEIGHBOR_COUNT];
		//float ix,iy,iz,jx,jy,jz;
		//int j_count=0;
		int muscleCounter = 0;
		int m_index[10640];
		float m_number[10640];
		float WXC = XMAX*0.5f;
		float WYC = YMAX*0.3f;
		float WZC = ZMAX*0.5f;
		//int sm_cnt = 0;
		//int array_k[MAX_NEIGHBOR_COUNT];
		for(i=numOfElasticP-numOfMembraneParticles;i<numOfElasticP;i++)
		{
			float dx2,dy2,dz2,r2_ij,r_ij;
			int k;
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
						elasticConnectionsData_cpp[ 4 * ( MAX_NEIGHBOR_COUNT * i + ecc) + 1 ] = r_ij*simulationScale*0.95;	// resting distance; that's why we use float type for elasticConnectionsData_cpp
						elasticConnectionsData_cpp[ 4 * ( MAX_NEIGHBOR_COUNT * i + ecc) + 2 ] = 0;						// type of connection; 0 - ordinary spring, 1 - muscle
						elasticConnectionsData_cpp[ 4 * ( MAX_NEIGHBOR_COUNT * i + ecc) + 3 ] = 0;						// not in use yet

						//define muscles
						if((position_cpp[ 4 * i + 2 ]<WZC+r0*95)&&(position_cpp[ 4 * j + 2 ]<WZC+r0*95))
						if((position_cpp[ 4 * i + 2 ]>WZC-r0*92)&&(position_cpp[ 4 * j + 2 ]>WZC-r0*92))
						//if(position_cpp[ 4 * i + 0 ]>XMAX*0.5)
						if( (fabs(position_cpp[4*i+3]-2.2f)<=0.05) && (fabs(position_cpp[4*j+3]-2.2f)<=0.05) )//both green
						{
							if((dz2>4*dx2)&&(dz2>4*dy2)&&(dx2>4*dy2))
							{
								muscle_color = 1.1f;

								//DR quadrant muscles
								
								if(position_cpp[ 4 * i + 0 ]>WXC)
								{
									if((position_cpp[4*i+1]<WYC)      &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*97)&&(position_cpp[4*j+2]<WZC+r0*97)) //MDR01
									if((position_cpp[4*i+2]>WZC+r0*85.9)&&(position_cpp[4*j+2]>WZC+r0*85.9)) muscle_color = 1.2f;//x.2 = red

									if((position_cpp[4*i+1]<WYC-1*r0) &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*95.0)&&(position_cpp[4*j+2]<WZC+r0*95.0)) //MDR02
									if((position_cpp[4*i+2]>WZC+r0*83.5)&&(position_cpp[4*j+2]>WZC+r0*83.5)) muscle_color = 2.4f;//x.4 = magenta

									if((position_cpp[4*i+1]<WYC)      &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*86.5)&&(position_cpp[4*j+2]<WZC+r0*86.5)) //MDR03
									if((position_cpp[4*i+2]>WZC+r0*77.5)&&(position_cpp[4*j+2]>WZC+r0*77.5)) muscle_color = 3.3f;//x.3 = orange

									if((position_cpp[4*i+1]<WYC-1*r0) &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*84.5)&&(position_cpp[4*j+2]<WZC+r0*84.5)) //MDR04
									if((position_cpp[4*i+2]>WZC+r0*76.5)&&(position_cpp[4*j+2]>WZC+r0*76.5)) muscle_color = 4.5f;
									if((position_cpp[4*i+1]<WYC-2*r0) &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*82.5)&&(position_cpp[4*j+2]<WZC+r0*82.5)) 
									if((position_cpp[4*i+2]>WZC+r0*72.5)&&(position_cpp[4*j+2]>WZC+r0*72.5)) muscle_color = 4.5f;//x.5 = violet
								
									if((position_cpp[4*i+1]<WYC)      &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*78)&&(position_cpp[4*j+2]<WZC+r0*78)) //MDR05
									if((position_cpp[4*i+2]>WZC+r0*66.9)&&(position_cpp[4*j+2]>WZC+r0*66.9)) muscle_color = 5.2f;
									if((position_cpp[4*i+1]<WYC-1*r0) &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*77.5)&&(position_cpp[4*j+2]<WZC+r0*77.5)) 
									if((position_cpp[4*i+2]>WZC+r0*65.9)&&(position_cpp[4*j+2]>WZC+r0*65.9)) muscle_color = 5.2f;//x.2 = red

									if((position_cpp[4*i+1]<WYC-2*r0) &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*74.0)&&(position_cpp[4*j+2]<WZC+r0*74.0)) //MDR06
									if((position_cpp[4*i+2]>WZC+r0*55.0)&&(position_cpp[4*j+2]>WZC+r0*55.0)) muscle_color = 6.4f;
									if((position_cpp[4*i+1]<WYC-3*r0) &&(position_cpp[4*i+1]>WYC-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*74.0)&&(position_cpp[4*j+2]<WZC+r0*74.0)) 
									if((position_cpp[4*i+2]>WZC+r0*54.5)&&(position_cpp[4*j+2]>WZC+r0*54.5)) muscle_color = 6.4f;//x.4 = magenta

									if((position_cpp[4*i+1]<WYC)        &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*68.5)&&(position_cpp[4*j+2]<WZC+r0*68.5)) //MDR07
									if((position_cpp[4*i+2]>WZC+r0*51.0)&&(position_cpp[4*j+2]>WZC+r0*51.0)) muscle_color = 7.3f;
									if((position_cpp[4*i+1]<WYC-1*r0)   &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*66.5)&&(position_cpp[4*j+2]<WZC+r0*66.5)) 
									if((position_cpp[4*i+2]>WZC+r0*49.5)&&(position_cpp[4*j+2]>WZC+r0*49.5)) muscle_color = 7.3f;//x.3 = orange

									if((position_cpp[4*i+1]<WYC-2*r0) &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*56.5)&&(position_cpp[4*j+2]<WZC+r0*56.5)) //MDR08
									if((position_cpp[4*i+2]>WZC+r0*40.0)&&(position_cpp[4*j+2]>WZC+r0*40.0)) muscle_color = 8.5f;
									if((position_cpp[4*i+1]<WYC-3*r0) &&(position_cpp[4*i+1]>WYC-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*55.5)&&(position_cpp[4*j+2]<WZC+r0*55.5)) 
									if((position_cpp[4*i+2]>WZC+r0*38.5)&&(position_cpp[4*j+2]>WZC+r0*38.5)) muscle_color = 8.5f;//x.5 = violet

									if((position_cpp[4*i+1]<WYC)        &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*51.5)&&(position_cpp[4*j+2]<WZC+r0*51.5)) //MDR09
									if((position_cpp[4*i+2]>WZC+r0*33.5)&&(position_cpp[4*j+2]>WZC+r0*33.5)) muscle_color = 9.2f;
									if((position_cpp[4*i+1]<WYC-1*r0)   &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*50.0)&&(position_cpp[4*j+2]<WZC+r0*50.0)) 
									if((position_cpp[4*i+2]>WZC+r0*33.0)&&(position_cpp[4*j+2]>WZC+r0*33.0)) muscle_color = 9.2f;//x.2 = red

									if((position_cpp[4*i+1]<WYC-2*r0) &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*40.5)&&(position_cpp[4*j+2]<WZC+r0*40.5)) //MDR10
									if((position_cpp[4*i+2]>WZC+r0*22.5)&&(position_cpp[4*j+2]>WZC+r0*22.5)) muscle_color = 10.4f;
									if((position_cpp[4*i+1]<WYC-3*r0) &&(position_cpp[4*i+1]>WYC-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*40.0)&&(position_cpp[4*j+2]<WZC+r0*40.0)) 
									if((position_cpp[4*i+2]>WZC+r0*21.5)&&(position_cpp[4*j+2]>WZC+r0*21.5)) muscle_color = 10.4f;
									if((position_cpp[4*i+1]<WYC-4*r0) &&(position_cpp[4*i+1]>WYC-5*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*40.0)&&(position_cpp[4*j+2]<WZC+r0*40.0)) 
									if((position_cpp[4*i+2]>WZC+r0*20.5)&&(position_cpp[4*j+2]>WZC+r0*20.5)) muscle_color = 10.4f;//x.4 = magenta

									if((position_cpp[4*i+1]<WYC)        &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*34.5)&&(position_cpp[4*j+2]<WZC+r0*34.5)) //MDR11
									if((position_cpp[4*i+2]>WZC+r0*15.5)&&(position_cpp[4*j+2]>WZC+r0*15.5)) muscle_color = 11.3f;
									if((position_cpp[4*i+1]<WYC-1*r0)   &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*33.5)&&(position_cpp[4*j+2]<WZC+r0*33.5)) 
									if((position_cpp[4*i+2]>WZC+r0*14.5)&&(position_cpp[4*j+2]>WZC+r0*14.5)) muscle_color = 11.3f;//x.3 = orange

									if((position_cpp[4*i+1]<WYC-2*r0) &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*23.5)&&(position_cpp[4*j+2]<WZC+r0*23.5)) //MDR12
									if((position_cpp[4*i+2]>WZC+r0* 8.5)&&(position_cpp[4*j+2]>WZC+r0* 8.5)) muscle_color = 12.5f;
									if((position_cpp[4*i+1]<WYC-3*r0) &&(position_cpp[4*i+1]>WYC-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*22.5)&&(position_cpp[4*j+2]<WZC+r0*22.5)) 
									if((position_cpp[4*i+2]>WZC+r0* 7.5)&&(position_cpp[4*j+2]>WZC+r0* 7.5)) muscle_color = 12.5f;
									if((position_cpp[4*i+1]<WYC-4*r0) &&(position_cpp[4*i+1]>WYC-5*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*21.5)&&(position_cpp[4*j+2]<WZC+r0*21.5)) 
									if((position_cpp[4*i+2]>WZC+r0* 6.5)&&(position_cpp[4*j+2]>WZC+r0* 6.5)) muscle_color = 12.5f;//x.5 = violet

									if((position_cpp[4*i+1]<WYC)        &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*16.0)&&(position_cpp[4*j+2]<WZC+r0*16.0)) //MDR13
									if((position_cpp[4*i+2]>WZC+r0* 1.5)&&(position_cpp[4*j+2]>WZC+r0* 1.5)) muscle_color = 13.2f;
									if((position_cpp[4*i+1]<WYC-1*r0)   &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0*15.5)&&(position_cpp[4*j+2]<WZC+r0*15.5)) 
									if((position_cpp[4*i+2]>WZC+r0* 0.5)&&(position_cpp[4*j+2]>WZC+r0* 0.5)) muscle_color = 13.2f;//x.2 = red

									if((position_cpp[4*i+1]<WYC-2*r0) &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 9.0)&&(position_cpp[4*j+2]<WZC+r0* 9.0)) //MDR14
									if((position_cpp[4*i+2]>WZC-r0* 2.5)&&(position_cpp[4*j+2]>WZC-r0* 2.5)) muscle_color = 14.4f;
									if((position_cpp[4*i+1]<WYC-3*r0) &&(position_cpp[4*i+1]>WYC-4*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 8.5)&&(position_cpp[4*j+2]<WZC+r0* 8.5)) 
									if((position_cpp[4*i+2]>WZC-r0* 3.5)&&(position_cpp[4*j+2]>WZC-r0* 3.5)) muscle_color = 14.4f;
									if((position_cpp[4*i+1]<WYC-4*r0) &&(position_cpp[4*i+1]>WYC-5*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 7.5)&&(position_cpp[4*j+2]<WZC+r0* 7.5)) 
									if((position_cpp[4*i+2]>WZC-r0* 4.5)&&(position_cpp[4*j+2]>WZC-r0* 4.5)) muscle_color = 14.4f;//x.4 = magenta

									if((position_cpp[4*i+1]<WYC)        &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 2.0)&&(position_cpp[4*j+2]<WZC+r0* 2.0)) //MDR15
									if((position_cpp[4*i+2]>WZC-r0*14.5)&&(position_cpp[4*j+2]>WZC-r0*14.5)) muscle_color = 15.3f;
									if((position_cpp[4*i+1]<WYC-1*r0)   &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC+r0* 1.0)&&(position_cpp[4*j+2]<WZC+r0* 1.0)) 
									if((position_cpp[4*i+2]>WZC-r0*15.5)&&(position_cpp[4*j+2]>WZC-r0*15.5)) muscle_color = 15.3f;//x.3 = orange

									if((position_cpp[4*i+1]<WYC-2*r0) &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0* 1.5)&&(position_cpp[4*j+2]<WZC-r0* 1.5)) //MDR16
									if((position_cpp[4*i+2]>WZC-r0*21.5)&&(position_cpp[4*j+2]>WZC-r0*21.5)) muscle_color = 16.5f;
									if((position_cpp[4*i+1]<WYC-3*r0) &&(position_cpp[4*i+1]>WYC-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0* 2.5)&&(position_cpp[4*j+2]<WZC-r0* 2.5)) 
									if((position_cpp[4*i+2]>WZC-r0*22.5)&&(position_cpp[4*j+2]>WZC-r0*22.5)) muscle_color = 16.5f;
									if((position_cpp[4*i+1]<WYC-4*r0) &&(position_cpp[4*i+1]>WYC-5*r0)) 
									if((position_cpp[4*i+2]<WZC-r0* 3.5)&&(position_cpp[4*j+2]<WZC-r0* 3.5)) 
									if((position_cpp[4*i+2]>WZC-r0*23.5)&&(position_cpp[4*j+2]>WZC-r0*23.5)) muscle_color = 16.5f;//x.5 = violet

									if((position_cpp[4*i+1]<WYC)        &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*14.0)&&(position_cpp[4*j+2]<WZC-r0*14.0)) //MDR17
									if((position_cpp[4*i+2]>WZC-r0*34.5)&&(position_cpp[4*j+2]>WZC-r0*34.5)) muscle_color = 17.2f;
									if((position_cpp[4*i+1]<WYC-1*r0)   &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*15.0)&&(position_cpp[4*j+2]<WZC-r0*15.0)) 
									if((position_cpp[4*i+2]>WZC-r0*35.5)&&(position_cpp[4*j+2]>WZC-r0*35.5)) muscle_color = 17.2f;//x.2 = red

									if((position_cpp[4*i+1]<WYC-2*r0) &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*20.0)&&(position_cpp[4*j+2]<WZC-r0*20.0)) //MDR18
									if((position_cpp[4*i+2]>WZC-r0*40.5)&&(position_cpp[4*j+2]>WZC-r0*40.5)) muscle_color = 18.4f;
									if((position_cpp[4*i+1]<WYC-3*r0) &&(position_cpp[4*i+1]>WYC-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*21.5)&&(position_cpp[4*j+2]<WZC-r0*21.5)) 
									if((position_cpp[4*i+2]>WZC-r0*41.5)&&(position_cpp[4*j+2]>WZC-r0*41.5)) muscle_color = 18.4f;
									if((position_cpp[4*i+1]<WYC-4*r0) &&(position_cpp[4*i+1]>WYC-5*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*22.5)&&(position_cpp[4*j+2]<WZC-r0*22.5)) 
									if((position_cpp[4*i+2]>WZC-r0*34.5)&&(position_cpp[4*j+2]>WZC-r0*34.5)) muscle_color = 18.4f;//x.4 = magenta

									if((position_cpp[4*i+1]<WYC)        &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*34.0)&&(position_cpp[4*j+2]<WZC-r0*34.0)) //MDR19
									if((position_cpp[4*i+2]>WZC-r0*54.5)&&(position_cpp[4*j+2]>WZC-r0*54.5)) muscle_color = 19.3f;
									if((position_cpp[4*i+1]<WYC-1*r0)   &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*34.5)&&(position_cpp[4*j+2]<WZC-r0*34.5)) 
									if((position_cpp[4*i+2]>WZC-r0*55.5)&&(position_cpp[4*j+2]>WZC-r0*55.5)) muscle_color = 19.3f;//x.3 = orange

									if((position_cpp[4*i+1]<WYC-2*r0) &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*39.5)&&(position_cpp[4*j+2]<WZC-r0*39.5)) //MDR20
									if((position_cpp[4*i+2]>WZC-r0*50.5)&&(position_cpp[4*j+2]>WZC-r0*50.5)) muscle_color = 20.5f;
									if((position_cpp[4*i+1]<WYC-3*r0) &&(position_cpp[4*i+1]>WYC-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*40.5)&&(position_cpp[4*j+2]<WZC-r0*40.5)) 
									if((position_cpp[4*i+2]>WZC-r0*51.5)&&(position_cpp[4*j+2]>WZC-r0*51.5)) muscle_color = 20.5f;//x.5 = violet

									if((position_cpp[4*i+1]<WYC)        &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*53.0)&&(position_cpp[4*j+2]<WZC-r0*53.0)) //MDR21
									if((position_cpp[4*i+2]>WZC-r0*71.5)&&(position_cpp[4*j+2]>WZC-r0*71.5)) muscle_color = 21.2f;
									if((position_cpp[4*i+1]<WYC-1*r0)   &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*54.0)&&(position_cpp[4*j+2]<WZC-r0*54.0)) 
									if((position_cpp[4*i+2]>WZC-r0*72.5)&&(position_cpp[4*j+2]>WZC-r0*72.5)) muscle_color = 21.2f;//x.2 = red

									if((position_cpp[4*i+1]<WYC-2*r0)   &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*50.0)&&(position_cpp[4*j+2]<WZC-r0*50.0)) //MDR22
									if((position_cpp[4*i+2]>WZC-r0*63.5)&&(position_cpp[4*j+2]>WZC-r0*63.5)) muscle_color = 22.4f;
									if((position_cpp[4*i+1]<WYC-3*r0)   &&(position_cpp[4*i+1]>WYC-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*51.0)&&(position_cpp[4*j+2]<WZC-r0*51.0)) 
									if((position_cpp[4*i+2]>WZC-r0*64.5)&&(position_cpp[4*j+2]>WZC-r0*64.5)) muscle_color = 22.4f;//x.4 = magenta

									if((position_cpp[4*i+1]<WYC)        &&(position_cpp[4*i+1]>WYC-1*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*70.0)&&(position_cpp[4*j+2]<WZC-r0*70.0)) //MDR23
									if((position_cpp[4*i+2]>WZC-r0*91.5)&&(position_cpp[4*j+2]>WZC-r0*91.5)) muscle_color = 23.3f;//x.3 = orange

									if((position_cpp[4*i+1]<WYC-1*r0) &&(position_cpp[4*i+1]>WYC-2*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*71.5)&&(position_cpp[4*j+2]<WZC-r0*71.5)) //MDR24
									if((position_cpp[4*i+2]>WZC-r0*91.5)&&(position_cpp[4*j+2]>WZC-r0*91.5)) muscle_color = 24.5f;
									if((position_cpp[4*i+1]<WYC-2*r0) &&(position_cpp[4*i+1]>WYC-3*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*62.5)&&(position_cpp[4*j+2]<WZC-r0*62.5)) 
									if((position_cpp[4*i+2]>WZC-r0*82.5)&&(position_cpp[4*j+2]>WZC-r0*82.5)) muscle_color = 24.5f;
									if((position_cpp[4*i+1]<WYC-3*r0) &&(position_cpp[4*i+1]>WYC-4*r0)) 
									if((position_cpp[4*i+2]<WZC-r0*63.5)&&(position_cpp[4*j+2]<WZC-r0*63.5)) 
									if((position_cpp[4*i+2]>WZC-r0*66.0)&&(position_cpp[4*j+2]>WZC-r0*66.0)) muscle_color = 24.5f;//x.5 = violet

								}

								elasticConnectionsData_cpp[ m_index[muscleCounter] = 4 * ( MAX_NEIGHBOR_COUNT * i + ecc) + 2 ] = muscle_color;// type of connection; 0 - ordinary spring, 1 or more - muscle
								muscleCounter++;
							}
						}
						array_j[ecc] = j;
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
					_mc = _mc;
				}
			}		
		}
	}



	return;
}

void owHelper::preLoadConfiguration()
{
	try
	{
		PARTICLE_COUNT = 0;
		ifstream positionFile ("./configuration/position.txt");
		int i = 0;
		float x, y, z, p_type;
		if( positionFile.is_open() )
		{
			while( positionFile.good() )
			{
				p_type = -1.1f;//reinitialize 
				positionFile >> x >> y >> z >> p_type;
				if(p_type>=0) PARTICLE_COUNT++;//last line of a file can contain only "\n", then p_type thanks to reinitialization will indicate the problem via negative value
				else break;//end of file
			}
		}

		PARTICLE_COUNT_RoundedUp = ((( PARTICLE_COUNT - 1 ) / local_NDRange_size ) + 1 ) * local_NDRange_size;

		printf("\nConfiguration we are going to load contains %d particles. Now plan to allocate memory for them.\n",PARTICLE_COUNT);
	}
	catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}

void owHelper::loadConfiguration(float *position_cpp, float *velocity_cpp, float *& elasticConnections,int & numOfLiquidP, int & numOfElasticP, int & numOfBoundaryP, int & numOfElasticConnections)
{

	try
	{
		ifstream positionFile ("./configuration/position.txt");
		int i = 0;
		float x, y, z, p_type;
		if( positionFile.is_open() )
		{
			while( positionFile.good() && i < PARTICLE_COUNT )
			{
				positionFile >> x >> y >> z >> p_type;
				position_cpp[ 4 * i + 0 ] = x;
				position_cpp[ 4 * i + 1 ] = y;
				position_cpp[ 4 * i + 2 ] = z;
				position_cpp[ 4 * i + 3 ] = p_type;
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
				i++;
			}
			positionFile.close();
		}
		else 
			throw std::runtime_error("Could not open file position.txt");
		ifstream velocityFile ("./configuration/velocity.txt");
		i = 0;
		if( velocityFile.is_open() )
		{
			while( velocityFile.good() && i < PARTICLE_COUNT )
			{
				velocityFile >> x >> y >> z >> p_type;
				velocity_cpp[ 4 * i + 0 ] = x;
				velocity_cpp[ 4 * i + 1 ] = y;
				velocity_cpp[ 4 * i + 2 ] = z;
				velocity_cpp[ 4 * i + 3 ] = p_type;
				i++;
			}
			velocityFile.close();
		}
		else 
			throw std::runtime_error("Could not open file velocity.txt");
		//TODO NEXT BLOCK WILL BE new load of elastic connections
		if(numOfElasticP != 0){
			ifstream elasticConectionsFile ("./configuration/elasticconnections.txt");
			elasticConnections = new float[ 4 * numOfElasticP * MAX_NEIGHBOR_COUNT ];
			/*int numElasticConnections = 0;
			for(i=0;i<numOfElasticP * MAX_NEIGHBOR_COUNT;i++)
			{
				elasticConnections[ 4 * i + 0 ] = NO_PARTICLE_ID;
			}*/
			i = 0;
			float  jd, rij0, val1, val2;// Elastic connection particle jd - jparticle rij0 - distance between i and j, val1, val2 - doesn't have any useful information yet
			if( elasticConectionsFile.is_open() )
			{
				//elasticConectionsFile >> numElasticConnections;

				while( elasticConectionsFile.good())
				{
					jd = -10;
					elasticConectionsFile >> jd >> rij0 >> val1 >> val2;
					if(jd>=-1)
					{
						elasticConnections[ 4 * i + 0 ] = jd;
						elasticConnections[ 4 * i + 1 ] = rij0;
						elasticConnections[ 4 * i + 2 ] = val1;
						elasticConnections[ 4 * i + 3 ] = val2;
						i++;
					}
				}
			}
		}
	}catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}
//This Function is currently on testing stage
void owHelper::loadConfigurationFromOneFile(float *position_cpp, float *velocity_cpp, float *&elasticConnections, int &numOfLiquidP, int &numOfElasticP, int &numOfBoundaryP, int &numOfElasticConnections)
{
	try
	{
		ifstream configurationFile ("./configuration/configuration.txt");
		int i = 0;
		float x, y, z, p_type = -1.f;
		std::string line;
		const int isPositionBlock = 1;
		const int isVelocityBlock = 2;
		const int isElasticConnectionsBlock = 3;
		int block = 0;
		bool isNotString = true;
		if( configurationFile.is_open() )
		{
			bool firstString = true;
			while( configurationFile.good() && i < PARTICLE_COUNT )
			{
				std::getline(configurationFile, line);
				std::istringstream iss(line);
				isNotString = true;
				//iss >> numOfElasticConnections;
				if (!(iss >> x >> y >> z >> p_type)) { 
					if(line == "Position"){
						block = isPositionBlock;
						i = 0;
						isNotString = false;
					}
					if( line == "Velocity"){
						block = isVelocityBlock;
						i = 0;
						isNotString = false;
					}
					if( line =="ElasticConnection"){
						block = isElasticConnectionsBlock;
						i = 0;
						isNotString = false;
					}
				} if(isNotString){
					switch(block){
						case isPositionBlock: {
							position_cpp[ 4 * i + 0 ] = x;
							position_cpp[ 4 * i + 1 ] = y;
							position_cpp[ 4 * i + 2 ] = z;
							position_cpp[ 4 * i + 3 ] = p_type;
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
							i++;
							break;
						}
						case isVelocityBlock: {
							velocity_cpp[ 4 * i + 0 ] = x;
							velocity_cpp[ 4 * i + 1 ] = y;
							velocity_cpp[ 4 * i + 2 ] = z;
							velocity_cpp[ 4 * i + 3 ] = p_type;
							i++;
							break;
						}
						case isElasticConnectionsBlock: {
							if(firstString){
								numOfElasticConnections = (int)x;//TODO write Comments here
								elasticConnections = new float[ 4 * numOfElasticConnections ];
								firstString = false;//on fist string we save count of all elastic connection
							}else if (i < numOfElasticConnections){
								elasticConnections[ 4 * i + 0 ] = x;//id;
								elasticConnections[ 4 * i + 1 ] = y;//jd;
								elasticConnections[ 4 * i + 2 ] = z;//rij0;
								elasticConnections[ 4 * i + 3 ] = p_type;//val;
								i++;
							}
							break;
						}
					}
				}
			}
			configurationFile.close();
		}
		else 
			throw std::runtime_error("Could not open file configuration.txt");
	}catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}
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
#endif
}
