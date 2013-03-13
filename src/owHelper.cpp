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


void owHelper::generateConfiguration(int stage, float *position, float *velocity, float *& elasticConnections,int & numOfLiquidP, int & numOfElasticP, int & numOfBoundaryP, int & numOfElasticConnections)
{
	// we need to know at least 
	// 1) sizes of the box which contains the simulation within
	// [0..MAX_X],[0..MAX_Y],[0..MAX_Z]
	// 2) smoothing radius
	// local vectors system and template data for each structure we are going to generate here

	float x,y,z;
	float p_type = LIQUID_PARTICLE;
	int i = 0;// particle counter
	int ix,iy,iz;
	int ecc = 0;//elastic connections counter

	int nx = (int)( ( XMAX - XMIN ) / r0 ); //X
	int ny = (int)( ( YMAX - YMIN ) / r0 ); //Y
	int nz = (int)( ( ZMAX - ZMIN ) / r0 ); //Z

	int nEx = 7;
	int nEy = 3;
	int nEz = 23;

	if(stage==0)
	{
		numOfLiquidP = 0;
		numOfElasticP = nEx*nEy*nEz;
		numOfBoundaryP = 0;

		elasticConnections = new float[ 4 * numOfElasticP * NEIGHBOR_COUNT ];
	}

	//=============== create elastic particles ==================================================
	if(stage==1)
	{
		p_type = ELASTIC_PARTICLE;

		for(x=0;x<nEx;x+=1.f)
		for(y=0;y<nEy;y+=1.f)
		for(z=0;z<nEz;z+=1.f)
		{
			//write particle coordinates to corresponding arrays
			position[ 4 * i + 0 ] = XMAX/2+x*r0-nEx*r0/2;
			position[ 4 * i + 1 ] = YMAX/2+y*r0-nEy*r0/2 + YMAX*3/8;
			position[ 4 * i + 2 ] = ZMAX/2+z*r0-nEz*r0/2;
			position[ 4 * i + 3 ] = p_type;

			velocity[ 4 * i + 0 ] = 0;
			velocity[ 4 * i + 1 ] = 0;
			velocity[ 4 * i + 2 ] = 0;
			velocity[ 4 * i + 3 ] = p_type;

			i++;
		}


		for(int i_ec = 0; i_ec < numOfElasticP * NEIGHBOR_COUNT; i_ec++)
		{
			elasticConnections[ 4 * i_ec + 0 ] = NO_PARTICLE_ID;
			elasticConnections[ 4 * i_ec + 1 ] = 0;
			elasticConnections[ 4 * i_ec + 2 ] = 0;
			elasticConnections[ 4 * i_ec + 3 ] = 0;
		}

		float r2ij;
		float dx2,dy2,dz2;

		for(int i_ec = 0; i_ec < numOfElasticP; i_ec++)
		{
			ecc = 0;
			float test;

			for(int j_ec = 0; j_ec < numOfElasticP; j_ec++)
			{
				if(i_ec!=j_ec)
				{
					dx2 = (position[ 4 * i_ec + 0 ] - position[ 4 * j_ec + 0 ]);
					dy2 = (position[ 4 * i_ec + 1 ] - position[ 4 * j_ec + 1 ]);
					dz2 = (position[ 4 * i_ec + 2 ] - position[ 4 * j_ec + 2 ]);
					dx2 *= dx2;
					dy2 *= dy2;
					dz2 *= dz2;
					r2ij = dx2 + dy2 + dz2;

					if(r2ij<=r0*r0*3.05f)
					{
						elasticConnections[ 4 * ( NEIGHBOR_COUNT * i_ec + ecc) + 0 ] = ((float)j_ec) + 0.1f;//connect elastic particles 0 and 1
						elasticConnections[ 4 * ( NEIGHBOR_COUNT * i_ec + ecc) + 1 ] = (float)sqrt(r2ij)*simulationScale;
						elasticConnections[ 4 * ( NEIGHBOR_COUNT * i_ec + ecc) + 2 ] = test = 0 + 1.1*((dz2>100*dx2)&&(dz2>100*dy2));
						elasticConnections[ 4 * ( NEIGHBOR_COUNT * i_ec + ecc) + 3 ] = 0;
						ecc++;
					}

					if(ecc>=NEIGHBOR_COUNT) break;
				}
			}
		}

		//and connections between them
		/*
		elasticConnections[ 4 * 0 + 0 ] = 1.1f;//connect elastic particles 0 and 1
		elasticConnections[ 4 * 0 + 1 ] = r0*simulationScale;
		elasticConnections[ 4 * 0 + 2 ] = 0;
		elasticConnections[ 4 * 0 + 3 ] = 0;

		
		elasticConnections[ 4 * NEIGHBOR_COUNT + 0 ] = 0.1f;//connect elastic particles 0 and 1
		elasticConnections[ 4 * NEIGHBOR_COUNT + 1 ] = r0*simulationScale;
		elasticConnections[ 4 * NEIGHBOR_COUNT + 2 ] = 0;
		elasticConnections[ 4 * NEIGHBOR_COUNT + 3 ] = 0;
		*/
	}

	//============= create volume of liquid =========================================================================
	p_type = LIQUID_PARTICLE;

	for(x = 15*r0/2;x<(XMAX-XMIN)/5 +3*r0/2;x += r0)
	for(y =  3*r0/2;y<(YMAX-YMIN)   -3*r0/2;y += r0)
	for(z =  3*r0/2+(ZMAX-ZMIN)*1/5;z<(ZMAX-ZMIN)*4/5-3*r0/2;z += r0)
	{
						// stage==0 - preliminary run
		if(stage==1)	// stage==1 - final run
		{
			if(i>=numOfLiquidP+numOfElasticP) 
			{
				printf("\nWarning! Final particle count >= preliminary particle count!\n");
				exit(-3);
			}
			//write particle coordinates to corresponding arrays
			position[ 4 * i + 0 ] = x;
			position[ 4 * i + 1 ] = y;
			position[ 4 * i + 2 ] = z;
			position[ 4 * i + 3 ] = p_type;

			velocity[ 4 * i + 0 ] = 0;
			velocity[ 4 * i + 1 ] = 0;
			velocity[ 4 * i + 2 ] = 0;
			velocity[ 4 * i + 3 ] = p_type;//if particle type is already defined in 'position', we don't need its duplicate here, right?
		}

		i++; // necessary for both stages
	}
	// end

	if(stage==0) 
	{
		numOfLiquidP = i;// - numOfElasticP;
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
						position[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position[ 4 * i + 2 ] =  0*r0 + r0/2;//z
						position[ 4 * i + 3 ] = p_type;
						velocity[ 4 * i + 0 ] = (1.f*(ix==0)-1*(ix==nx-1))/sqrt(3.f);//norm x
						velocity[ 4 * i + 1 ] = (1.f*(iy==0)-1*(iy==ny-1))/sqrt(3.f);//norm y
						velocity[ 4 * i + 2 ] =  1.f/sqrt(3.f);//norm z
						velocity[ 4 * i + 3 ] = p_type;
						i++;
						position[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position[ 4 * i + 2 ] = (nz-1)*r0 + r0/2;//z
						position[ 4 * i + 3 ] = p_type;
						velocity[ 4 * i + 0 ] = (1*(ix==0)-1*(ix==nx-1))/sqrt(3.f);//norm x
						velocity[ 4 * i + 1 ] = (1*(iy==0)-1*(iy==ny-1))/sqrt(3.f);//norm y
						velocity[ 4 * i + 2 ] = -1.f/sqrt(3.f);//norm z
						velocity[ 4 * i + 3 ] = p_type;
						i++;
					}
					else //edges
					{
						position[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position[ 4 * i + 2 ] =  0*r0 + r0/2;//z
						position[ 4 * i + 3 ] = p_type;
						velocity[ 4 * i + 0 ] =  1.f*((ix==0)-(ix==nx-1))/sqrt(2.f);//norm x
						velocity[ 4 * i + 1 ] =  1.f*((iy==0)-(iy==ny-1))/sqrt(2.f);//norm y
						velocity[ 4 * i + 2 ] =  1.f/sqrt(2.f);//norm z
						velocity[ 4 * i + 3 ] = p_type;
						i++;
						position[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position[ 4 * i + 2 ] = (nz-1)*r0 + r0/2;//z
						position[ 4 * i + 3 ] = p_type;
						velocity[ 4 * i + 0 ] =  1.f*((ix==0)-(ix==nx-1))/sqrt(2.f);//norm x
						velocity[ 4 * i + 1 ] =  1.f*((iy==0)-(iy==ny-1))/sqrt(2.f);//norm y
						velocity[ 4 * i + 2 ] = -1.f/sqrt(2.f);//norm z
						velocity[ 4 * i + 3 ] = p_type;
						i++;
					}
				}
				else //planes
				{
						position[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position[ 4 * i + 2 ] =  0*r0 + r0/2;//z
						position[ 4 * i + 3 ] = p_type;
						velocity[ 4 * i + 0 ] =  0;//norm x
						velocity[ 4 * i + 1 ] =  0;//norm y
						velocity[ 4 * i + 2 ] =  1;//norm z
						velocity[ 4 * i + 3 ] = p_type;
						i++;
						position[ 4 * i + 0 ] = ix*r0 + r0/2;//x
						position[ 4 * i + 1 ] = iy*r0 + r0/2;//y
						position[ 4 * i + 2 ] = (nz-1)*r0 + r0/2;//z
						position[ 4 * i + 3 ] = p_type;
						velocity[ 4 * i + 0 ] =  0;//norm x
						velocity[ 4 * i + 1 ] =  0;//norm y
						velocity[ 4 * i + 2 ] = -1;//norm z
						velocity[ 4 * i + 3 ] = p_type;
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
					position[ 4 * i + 0 ] = ix*r0 + r0/2;//x
					position[ 4 * i + 1 ] =  0*r0 + r0/2;//y
					position[ 4 * i + 2 ] = iz*r0 + r0/2;//z
					position[ 4 * i + 3 ] = p_type;
					velocity[ 4 * i + 0 ] =  0;//norm x
					velocity[ 4 * i + 1 ] =  1.f/sqrt(2.f);//norm y
					velocity[ 4 * i + 2 ] =  1.f*((iz==0)-(iz==nz-1))/sqrt(2.f);//norm z
					velocity[ 4 * i + 3 ] = p_type;
					i++;
					position[ 4 * i + 0 ] = ix*r0 + r0/2;//x
					position[ 4 * i + 1 ] = (ny-1)*r0 + r0/2;//y
					position[ 4 * i + 2 ] = iz*r0 + r0/2;//z
					position[ 4 * i + 3 ] = p_type;
					velocity[ 4 * i + 0 ] =  0;//norm x
					velocity[ 4 * i + 1 ] = -1.f/sqrt(2.f);//norm y
					velocity[ 4 * i + 2 ] =  1.f*((iz==0)-(iz==nz-1))/sqrt(2.f);//norm z
					velocity[ 4 * i + 3 ] = p_type;
					i++;
				}
				else //planes
				{
					position[ 4 * i + 0 ] = ix*r0 + r0/2;//x
					position[ 4 * i + 1 ] =  0*r0 + r0/2;//y
					position[ 4 * i + 2 ] = iz*r0 + r0/2;//z
					position[ 4 * i + 3 ] = p_type;
					velocity[ 4 * i + 0 ] =  0;//norm x
					velocity[ 4 * i + 1 ] =  1;//norm y
					velocity[ 4 * i + 2 ] =  0;//norm z
					velocity[ 4 * i + 3 ] = p_type;
					i++;
					position[ 4 * i + 0 ] = ix*r0 + r0/2;//x
					position[ 4 * i + 1 ] = (ny-1)*r0 + r0/2;//y
					position[ 4 * i + 2 ] = iz*r0 + r0/2;//z
					position[ 4 * i + 3 ] = p_type;
					velocity[ 4 * i + 0 ] =  0;//norm x
					velocity[ 4 * i + 1 ] = -1;//norm y
					velocity[ 4 * i + 2 ] =  0;//norm z
					velocity[ 4 * i + 3 ] = p_type;
					i++;
				}
			}
		}

		// 3 - side walls OY-OZ and opposite
		for(iy=1;iy<ny-1;iy++)
		{
			for(iz=1;iz<nz-1;iz++)
			{
				position[ 4 * i + 0 ] =  0*r0 + r0/2;//x
				position[ 4 * i + 1 ] = iy*r0 + r0/2;//y
				position[ 4 * i + 2 ] = iz*r0 + r0/2;//z
				position[ 4 * i + 3 ] = p_type;
				velocity[ 4 * i + 0 ] =  1;//norm x
				velocity[ 4 * i + 1 ] =  0;//norm y
				velocity[ 4 * i + 2 ] =  0;//norm z
				velocity[ 4 * i + 3 ] = p_type;
				i++;
				position[ 4 * i + 0 ] = (nx-1)*r0 + r0/2;//x
				position[ 4 * i + 1 ] = iy*r0 + r0/2;//y
				position[ 4 * i + 2 ] = iz*r0 + r0/2;//z
				position[ 4 * i + 3 ] = p_type;
				velocity[ 4 * i + 0 ] = -1;//norm x
				velocity[ 4 * i + 1 ] =  0;//norm y
				velocity[ 4 * i + 2 ] =  0;//norm z
				velocity[ 4 * i + 3 ] = p_type;
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

void owHelper::loadConfiguration(float *position, float *velocity, float *& elasticConnections,int & numOfLiquidP, int & numOfElasticP, int & numOfBoundaryP, int & numOfElasticConnections)
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
				position[ 4 * i + 0 ] = x;
				position[ 4 * i + 1 ] = y;
				position[ 4 * i + 2 ] = z;
				position[ 4 * i + 3 ] = p_type;
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
				velocity[ 4 * i + 0 ] = x;
				velocity[ 4 * i + 1 ] = y;
				velocity[ 4 * i + 2 ] = z;
				velocity[ 4 * i + 3 ] = p_type;
				i++;
			}
			velocityFile.close();
		}
		else 
			throw std::runtime_error("Could not open file velocity.txt");
		//TODO NEXT BLOCK WILL BE new load of elastic connections
		if(numOfElasticP != 0){
			ifstream elasticConectionsFile ("./configuration/elasticconnections.txt");
			elasticConnections = new float[ 4 * numOfElasticP * NEIGHBOR_COUNT ];
			/*int numElasticConnections = 0;
			for(i=0;i<numOfElasticP * NEIGHBOR_COUNT;i++)
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
void owHelper::loadConfigurationFromOneFile(float *position, float *velocity, float *&elasticConnections, int &numOfLiquidP, int &numOfElasticP, int &numOfBoundaryP, int &numOfElasticConnections)
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
							position[ 4 * i + 0 ] = x;
							position[ 4 * i + 1 ] = y;
							position[ 4 * i + 2 ] = z;
							position[ 4 * i + 3 ] = p_type;
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
							velocity[ 4 * i + 0 ] = x;
							velocity[ 4 * i + 1 ] = y;
							velocity[ 4 * i + 2 ] = z;
							velocity[ 4 * i + 3 ] = p_type;
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
