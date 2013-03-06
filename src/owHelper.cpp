#include "owHelper.h"
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

		PARTICLE_COUNT_RoundedUp = ((( PARTICLE_COUNT - 1 ) / 256 ) + 1 ) * 256;

		printf("\nConfiguration we are going to load contains %d particles. Now plan to allocate memory for them.\n",PARTICLE_COUNT);
	}
	catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}

void owHelper::loadConfiguration(float *position, float *velocity, float *& elasticConnections,int & numOfLiquedP, int & numOfElasticP, int & numOfBoundaryP, int & numOfElasticConnections)
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
						numOfLiquedP++;
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
					if(jd>=-5)
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
//This Function is testing
void owHelper::loadConfigurationFromOneFile(float *position, float *velocity, float *&elasticConnections, int &numOfLiquedP, int &numOfElasticP, int &numOfBoundaryP, int &numOfElasticConnections)
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
									numOfLiquedP++;
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
