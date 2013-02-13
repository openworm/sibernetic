#include "owHelper.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
using namespace std;

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
void owHelper::log_bufferf(const float * buffer, const int ellement_size, const int global_size, const char * fileName)
{
	try{
		ofstream outFile (fileName);
		for(int i = 0; i < global_size; i++)
		{
			for(int j = 0; j < ellement_size; j++)
			{
				if(j < ellement_size - 1 )
					outFile << buffer[ i * ellement_size + j ] << "\t";
				else
					outFile << buffer[ i * ellement_size + j ] << "\n";
			}
		}
		outFile.close();
	}catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
}
//For output int buffer
void owHelper::log_bufferi(const int * buffer, const int ellement_size, const int global_size, const char * fileName)
{
	try{
		ofstream outFile (fileName);
		for(int i = 0; i < global_size; i++)
		{
			for(int j = 0; j < ellement_size; j++)
			{
				if(j < ellement_size + 1 )
					outFile << buffer[ i * ellement_size + j ] << "\t";
				else
					outFile << buffer[ i * ellement_size + j ] << "\n";
			}
		}
		outFile.close();
	}catch(std::exception &e){
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
		if(numOfElasticP != 0){
			ifstream elasticConectionsFile ("./configuration/elasticconnections.txt");
			i = 0;
			float id, jd, rij0, val;// Ellastic connection info id - i partical jd - jparticle rij0 - distance between i and j, val - doesn't have any useful information we use it only for vectorization
			if( elasticConectionsFile.is_open() )
			{
				bool firstString = true;
				while( elasticConectionsFile.good())
				{
					if(firstString){
						elasticConectionsFile >> numOfElasticConnections;
						elasticConnections = new float[ 4 * numOfElasticConnections ];
						firstString = false;//on fist string we save cout of all ellastic connection
					}else if (i < numOfElasticConnections){
						elasticConectionsFile >> id >> jd >> rij0 >> val;
						elasticConnections[ 4 * i + 0 ] = id;
						elasticConnections[ 4 * i + 1 ] = jd;
						elasticConnections[ 4 * i + 2 ] = rij0;
						elasticConnections[ 4 * i + 3 ] = val;
						i++;
					}else{
						break;
					}
				}
				elasticConectionsFile.close();
			}
		}
		else {
			numOfElasticConnections = 0;
			//throw std::runtime_error("Could not open file elasticConnections.txt");
		}
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
