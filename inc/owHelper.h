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

#ifndef OW_HELPER_H
#define OW_HELPER_H

#include <fstream>
#include <stdexcept>
#include <iostream>
#include <cstdlib>

#include "owConfigProperty.h"
#include "owOpenCLConstant.h"
#include "owPhysicsConstant.h"

#if defined(_WIN32) || defined (_WIN64)
	#include <windows.h>
#elif defined(__linux__)
	#include <time.h>
#elif defined(__APPLE__)
    #include <stddef.h>
#endif
/** Helper class contains a number of helper methods.
  */
class owHelper
{
public:
	owHelper(void);
	~owHelper(void);
	static void preLoadConfiguration( owConfigProperty * config );
	static void loadConfiguration( float *position_cpp, float *velocity_cpp, float *& elasticConnections, int * membraneData_cpp, int *& particleMembranesList_cpp, owConfigProperty * config );
	static void loadConfigurationToFile(float * position, owConfigProperty * config,float * connections=NULL, int * membranes=NULL, bool firstIteration = true, int * filter_p=NULL, int size=0);
	static void loadConfigurationFromFile(float *& position, float *& connections, int *& membranes, int iteration = 0);
	static bool loadConfigurationFromFile(float *& position, float *& connections, int *& membranes, owConfigProperty * config,int iteration = 0);
	static void loadConfigurationFromFile(float *& position, float *& velocity,float *& connections, int *& membranes, int *& particleMemIndex, owConfigProperty * config);
	static void loadConfigurationToFile(float * position, float * velocity, float * connections, int * membranes, int * particleMemIndex, const char * filename, owConfigProperty * config);
	void watch_report(const char *str);
	double getElapsedTime() { return elapsedTime; };
	void refreshTime();
	//For output buffer
	//Create File in which line element_size elements
	//global_size - size of buffer / element_size
	template<typename T> static void log_buffer(const T * buffer, const int element_size, const int global_size, const char * fileName)
	{
		std::ofstream outFile (fileName);
		if(!outFile){
			std::string error = "Couldn't create file check the file name: ";
			error += fileName;
			throw std::runtime_error(error);
		}
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
	}
private:
	enum LOADMODE { NOMODE=-1, POSITION, VELOCITY, CONNECTION, MEMBRANE, PMEMINDEX };
	double elapsedTime;
#if defined(_WIN32) || defined (_WIN64)
	LARGE_INTEGER frequency;				// ticks per second
	LARGE_INTEGER t0, t1, t2;				// ticks
	LARGE_INTEGER t3,t4;
#elif defined(__linux__)
	timespec t0, t1, t2;
	timespec t3,t4;
	double us;
#elif defined(__APPLE__)
    unsigned long t0, t1, t2;
#endif
};
#endif // #ifndef OW_HELPER_H
