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
#include <algorithm>


#if defined(__APPLE__) || defined (__MACOSX)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

#include "owHelper.h"
#include "owPhysicsConstant.h"


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

int read_position = 0;
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
 *  @param config
 *  pointer to owConfigProperty object it includes information about
 *  boundary box dimensions
 */
void owHelper::preLoadConfiguration(owConfigProperty * config)
{
	int p_count = 0;
	std::string file_name = config->getCofigPath() + config->getCofigFileName();
	std::string inputStr;
	std::ifstream configFile (file_name.c_str(), std::ios_base::binary);
	float x, y, z, p_type;
	char delimiter = '\t';
	size_t pos;
	LOADMODE mode = NOMODE;
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
			inputStr.erase( std::remove( inputStr.begin(), inputStr.end(), '\r' ), inputStr.end() );
			inputStr.erase( std::remove( inputStr.begin(), inputStr.end(), '\n' ), inputStr.end() );
			if( inputStr == "[position]"){
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
							config->numOfLiquidP++;
							break;
						case ELASTIC_PARTICLE:
							config->numOfElasticP++;
							break;
						case BOUNDARY_PARTICLE:
							config->numOfBoundaryP++;
							break;
					}
					break;
				}
				case MEMBRANE:{
					config->numOfMembranes++;
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
 *  @param membraneData_cpp
 *  pointer to membraneData_cpp buffer
 *  @param particleMembranesList_cpp
 *  pointer to particleMembranesList_cpp buffer
 *  @param config
 *  pointer to owConfigProperty object it includes information about
 */
void owHelper::loadConfiguration(float *position_cpp, float *velocity_cpp, float *& elasticConnections, int * membraneData_cpp, int *& particleMembranesList_cpp, owConfigProperty * config)
{
	std::string file_name = config->getCofigPath() + config->getCofigFileName();
	std::string inputStr;
	std::ifstream configFile (file_name.c_str(), std::ios_base::binary);
	char delimiter = '\t';
	size_t pos;
	LOADMODE mode = NOMODE;
	int i = 0;
	if( configFile.is_open() )
	{
		configFile.seekg(read_position);
		while( configFile.good() )
		{
			//configFile >> inputStr;
			std::getline(configFile,inputStr);
			inputStr.erase( std::remove( inputStr.begin(), inputStr.end(), '\r' ), inputStr.end() );
			inputStr.erase( std::remove( inputStr.begin(), inputStr.end(), '\n' ), inputStr.end() );
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
 *  pointer to owConfigProperty object it includes information about
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
void owHelper::loadConfigurationToFile(float * position, owConfigProperty * config, float * connections, int * membranes, bool firstIteration, int * filter_p, int size ){
	std::ofstream positionFile;
	std::string positionFileName = config->getLoadPath() + std::string("/position_buffer.txt");
	if(firstIteration){
		positionFile.open(positionFileName.c_str(), std::ofstream::trunc);
		if(!positionFile)
			throw std::runtime_error("There was a problem with creation of position file for logging check the path.");
		positionFile << config->xmin << "\n";
		positionFile << config->xmax << "\n";
		positionFile << config->ymin << "\n";
		positionFile << config->ymax << "\n";
		positionFile << config->zmin << "\n";
		positionFile << config->zmax << "\n";
		positionFile << config->numOfElasticP << "\n";
		positionFile << config->numOfLiquidP << "\n";
		positionFile << config->numOfBoundaryP << "\n";
		positionFile << config->getTimeStep() << "\n";
		positionFile << config->getLogStep() << "\n";
	}else{
		positionFile.open(positionFileName.c_str(), std::ofstream::app);
		if(!positionFile)
			throw std::runtime_error("There was a problem with creation of position file for logging Check the path.");
	}
	if(size==0){
		for(int i=0;i < config->getParticleCount(); i++){
			if((int)position[ 4 * i + 3] != BOUNDARY_PARTICLE || firstIteration ){
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
		std::string connectionFileName = config->getLoadPath() + std::string("/connection_buffer.txt");
		std::ofstream connectionFile(connectionFileName.c_str(), std::ofstream::trunc);
		if(!connectionFile)
			throw std::runtime_error("There was a problem with creation of connection data file for logging. Check the path.");
		int con_num = MAX_NEIGHBOR_COUNT * config->numOfElasticP;
		for(int i = 0; i < con_num; i++)
			connectionFile << connections[4 * i + 0] << "\t" << connections[4 * i + 1] << "\t" << connections[4 * i + 2] << "\t" << connections[4 * i + 3] << "\n";
		connectionFile.close();
		std::string membraneFileName = config->getLoadPath() + std::string("/membranes_buffer.txt");
		std::ofstream membranesFile(membraneFileName.c_str(), std::ofstream::trunc);
		if(!membranesFile)
			throw std::runtime_error("There was a problem with creation of membrane data file for logging. Check the path.");
		membranesFile << config->numOfMembranes << "\n";
		for(unsigned int i = 0; i < config->numOfMembranes; i++)
			membranesFile << membranes[3 * i + 0] << "\t" << membranes[3 * i + 1] << "\t" << membranes[3 * i + 2] << "\n";
		membranesFile.close();
	}
}
/** Load configuration from simulation to files
 *
 *  Make configuration file
 *  @param position
 *  pointer to position buffer
 *  @param config
 *  pointer to owConfigProperty object it includes information about
 *  @param connections
 *  reference on pointer to elasticConnections buffer.
 *  @param membranes
 *  pointer to membranes buffer
 */
void owHelper::loadConfigurationToFile(float * position, float * velocity, float * connections, int * membranes, int * particleMemIndex,const char * filename, owConfigProperty * config){
	std::ofstream configFile(filename, std::ofstream::trunc);
	if(!configFile)
		throw std::runtime_error("There was a problem with creation of file for saving configuration. Check the path.");
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
	int con_num = MAX_NEIGHBOR_COUNT * config->numOfElasticP;
	for(int i = 0; i < con_num; i++)
		configFile << connections[4 * i + 0] << "\t" << connections[4 * i + 1] / simulationScale << "\t" << connections[4 * i + 2] << "\t" << connections[4 * i + 3] << "\n";
	configFile << "[membranes]\n";
	for(unsigned int i = 0; i < config->numOfMembranes; ++i)
		configFile << membranes[3 * i + 0] << "\t" << membranes[3 * i + 1] << "\t" << membranes[3 * i + 2] << "\n";
	configFile << "[particleMemIndex]\n";
	int particleMemIndexCount = config->numOfElasticP*MAX_MEMBRANES_INCLUDING_SAME_PARTICLE;
	for(int i = 0; i < particleMemIndexCount; ++i)
		configFile << particleMemIndex[i] << "\n";
	configFile << "[end]";
	configFile.close();
}
//This function needed for visualiazation buffered data
long position_index = 0;
std::ifstream positionFile;
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
 *  pointer to owConfigProperty object it includes information about
 *  @param iteration
 *  if iteration==0 it means that we first time record information
 *  to a file and on first iteration it put to
 *  the file info about dimensions of boundary box
 */
bool owHelper::loadConfigurationFromFile(float *& position, float *& connections, int *& membranes, owConfigProperty * config, int iteration){
		if(iteration == 0){
			std::string positionFileName = config->getLoadPath() + std::string("/position_buffer.txt");
			positionFile.open(positionFileName.c_str());
		}
		unsigned int i = 0;
		float x, y, z, p_type;
		if( positionFile.is_open() )
		{
			if(iteration == 0){
				float valueF;
				int valueI;
				positionFile >> config->xmin;
				positionFile >> config->xmax;
				positionFile >> config->ymin;
				positionFile >> config->ymax;
				positionFile >> config->zmin;
				positionFile >> config->zmax;
				positionFile >> config->numOfElasticP;
				positionFile >> config->numOfLiquidP;
				positionFile >> config->numOfBoundaryP;
				positionFile >> valueF;
				positionFile >> valueI;
				config->setTimeStep(valueF);
				config->setLogStep(valueI);
				config->setParticleCount(config->numOfElasticP + config->numOfLiquidP + config->numOfBoundaryP);
				position = new float[4 * config->getParticleCount()];
			}
			while( positionFile.good() &&  i < config->getParticleCount())
			{
				if(static_cast<int>(p_type) != BOUNDARY_PARTICLE || iteration == 0){
					positionFile >> x >> y >> z >> p_type;
					position[i * 4 + 0] = x;
					position[i * 4 + 1] = y;
					position[i * 4 + 2] = z;
					position[i * 4 + 3] = p_type;
				}
				i++;
			}
			if(iteration == 0)
				config->setParticleCount(config->numOfElasticP + config->numOfLiquidP);
		}
		if(!positionFile.good()){
			positionFile.close();
			return false;
		}
		if(iteration == 0){
			std::string connectionFileName = config->getLoadPath() + std::string("/connection_buffer.txt");
			std::ifstream connectionFile(connectionFileName.c_str());
			connections = new float[MAX_NEIGHBOR_COUNT * config->numOfElasticP * 4];
			if( connectionFile.is_open() )
			{
				i = 0;
				float jd, rij0, val1, val2;
				while(connectionFile.good() && i < MAX_NEIGHBOR_COUNT * config->numOfElasticP){
					connectionFile >> jd >> rij0 >> val1 >> val2;
					connections[ 4 * i + 0 ] = jd;
					connections[ 4 * i + 1 ] = rij0;
					connections[ 4 * i + 2 ] = val1;
					connections[ 4 * i + 3 ] = val2;
					i++;
				}
			}
			connectionFile.close();
			std::string membraneFileName = config->getLoadPath() + std::string("/membranes_buffer.txt");
			std::ifstream membranesFile(membraneFileName.c_str());
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
		return true;
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
