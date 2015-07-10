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


//READ DEFAULT CONFIGURATATION FROM FILE IN CONFIGURATION FOLDER TODO move it into configuration struct
int read_position = 0;
std::string owHelper::path = "./configuration/";
std::string owHelper::suffix = "";
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
		std::string p_file_name = path + "position" + suffix + ".txt";
		std::ifstream positionFile (p_file_name.c_str(), std::ios_base::binary);
		float x, y, z, p_type;
		if( positionFile.is_open() )
		{
			positionFile >> config->xmin;
			positionFile >> config->xmax;
			positionFile >> config->ymin;
			positionFile >> config->ymax;
			positionFile >> config->zmin;
			positionFile >> config->zmax;
			read_position = positionFile.tellg();
			while( positionFile.good() )
			{
				p_type = -1.1f;//reinitialize
				positionFile >> x >> y >> z >> p_type;
				if(p_type>=0){
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
				}//last line of a file can contain only "\n", then p_type thanks to reinitialization will indicate the problem via negative value
				else break;//end of file
			}
		}
		positionFile.close();
		config->setParticleCount(p_count);

		printf("\nConfiguration we are going to load contains %d particles. Now plan to allocate memory for them.\n",config->getParticleCount());

		numOfMembranes = 0;
		std::string m_file_name = path + "membranes" + suffix + ".txt";
		std::ifstream membranesFile (m_file_name.c_str());
		int id, jd, kd;
		if( membranesFile.is_open() )
		{
			while( membranesFile.good() )
			{
				kd = -1;
				membranesFile >> id >> jd >> kd ;
				if(kd>=0)numOfMembranes++;//last line of a file can contain only "\n", then kd thanks to reinitialization will indicate the problem via negative value
			}
		}
		membranesFile.close();
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
		std::string p_file_name = path + "position" + suffix + ".txt";
		std::ifstream positionFile (p_file_name.c_str());
		int i = 0;
		float x, y, z, p_type;
		if( positionFile.is_open() )
		{
			positionFile.seekg(read_position);
			while( positionFile.good() && i < config->getParticleCount() )
			{
				positionFile >> x >> y >> z >> p_type;
				position_cpp[ 4 * i + 0 ] = x;
				position_cpp[ 4 * i + 1 ] = y;
				position_cpp[ 4 * i + 2 ] = z;
				position_cpp[ 4 * i + 3 ] = p_type;
				i++;
			}
			positionFile.close();
		}
		else
			throw std::runtime_error("Could not open file position.txt");
		std::cout << "Position is loaded" << std::endl;
		std::string v_file_name = path + "velocity" + suffix + ".txt";
		std::ifstream velocityFile (v_file_name.c_str());
		i = 0;
		if( velocityFile.is_open() )
		{
			while( velocityFile.good() && i < config->getParticleCount() )
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
		std::cout << "Velocity is loaded" << std::endl;
		if(numOfElasticP != 0){
			std::string c_file_name = path + "connection" + suffix + ".txt";
			std::ifstream elasticConectionsFile (c_file_name.c_str());
			elasticConnections = new float[ 4 * numOfElasticP * MAX_NEIGHBOR_COUNT ];
			i = 0;
			if( elasticConectionsFile.is_open() )
			{
				float  jd, rij0, val1, val2;// Elastic connection particle jd - jparticle rij0 - distance between i and j, val1, val2 - doesn't have any useful information yet
				//elasticConectionsFile >> numElasticConnections;

				while( elasticConectionsFile.good())
				{
					jd = -10;
					elasticConectionsFile >> jd >> rij0 >> val1 >> val2;
					if(jd>=-1)
					{
						elasticConnections[ 4 * i + 0 ] = jd;
						elasticConnections[ 4 * i + 1 ] = rij0 * simulationScale;
						elasticConnections[ 4 * i + 2 ] = val1;
						elasticConnections[ 4 * i + 3 ] = val2;
						i++;
					}
				}
				elasticConectionsFile.close();
			}
			else
				throw std::runtime_error("Could not open file connection.txt");
			std::cout << "Elastic Connection is loaded" << std::endl;
			if(numOfMembranes != 0){
				//Import Membranes
				//return;
				std::string m_file_name = path + "membranes" + suffix + ".txt";
				std::ifstream membranesFile (m_file_name.c_str());
				i = 0;
				if( membranesFile.is_open() )
				{
					int id, jd, kd;

					while( membranesFile.good() && i < numOfMembranes)
					{
						membranesFile >> id >> jd >> kd;
						membraneData_cpp[ 3 * i + 0 ] = id;
						membraneData_cpp[ 3 * i + 1 ] = jd;
						membraneData_cpp[ 3 * i + 2 ] = kd;
						i++;
					}
					membranesFile.close();
				}
				else
					throw std::runtime_error("Could not open file membranes.txt");
				std::cout << "Membranes is loaded" << std::endl;
				//Import Membranes
				std::string mi_file_name = path + "particleMembraneIndex" + suffix + ".txt";
				std::ifstream membranesIndexFile (mi_file_name.c_str());
				i = 0;
				if( membranesIndexFile.is_open())
				{
					int id;
					while( membranesIndexFile.good() && i < numOfElasticP*MAX_MEMBRANES_INCLUDING_SAME_PARTICLE)
					{
						membranesIndexFile >> id ;
						particleMembranesList_cpp[ i ] = id;
						i++;
					}
					membranesIndexFile.close();
				}
				else
					throw std::runtime_error("Could not open file particleMembraneIndex.txt");
				std::cout << "ParticleMembraneIndex is loaded" << std::endl;
			}
		}
	}catch(std::exception &e){
		std::cout << "ERROR: " << e.what() << std::endl;
		exit( -1 );
	}
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
void owHelper::loadConfigurationToFile(float * position, owConfigProrerty * config, float * connections, int * membranes, bool firstIteration, int * filter_p, int size ){
	try{
		ofstream positionFile;
		if(firstIteration){
			positionFile.open("./buffers/position_buffer.txt", std::ofstream::trunc);
			positionFile << config->xmin << "\n";
			positionFile << config->xmax << "\n";
			positionFile << config->ymin << "\n";
			positionFile << config->ymax << "\n";
			positionFile << config->zmin << "\n";
			positionFile << config->zmax << "\n";
			positionFile << numOfElasticP << "\n";
			positionFile << numOfLiquidP << "\n";
		}else{
			positionFile.open("./buffers/position_buffer.txt", std::ofstream::app);
		}
		if(size==0){
			for(int i=0;i < config->getParticleCount(); i++){
				if((int)position[ 4 * i + 3] != BOUNDARY_PARTICLE){
					positionFile << position[i * 4 + 0] << "\t" << position[i * 4 + 1] << "\t" << position[i * 4 + 2] << "\t" << position[i * 4 + 3] << "\n";
				}
			}
		}else{
			int i = 0;
			int index = 0;
			while(index!=size){
				i = filter_p[index];
				positionFile << position[i * 4 + 0] << "\t" << position[i * 4 + 1] << "\t" << position[i * 4 + 2] << "\t" << position[i * 4 + 3] << "\n";
				index++;
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
				membranesFile << membranes[4 * i + 0] << "\t" << membranes[4 * i + 1] << "\t" << membranes[4 * i + 2] << "\t" << membranes[4 * i + 3] << "\n";
			membranesFile.close();
		}
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
void owHelper::loadConfigurationFromFile_experemental(float *& position, float *& connections, int *& membranes, owConfigProrerty * config, int iteration){
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
				membranesFile >> m_count;
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
