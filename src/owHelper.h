#ifndef OW_HELPER_H
#define OW_HELPER_H

#include "owOpenCLConstant.h"
#if defined(_WIN32) || defined (_WIN64)
	#include <windows.h>
#elif defined(__linux__)
	#include <time.h>
#endif

#if defined(__APPLE__)
    #include <stddef.h>
#endif

class owHelper
{
public:
	owHelper(void);
	~owHelper(void);
	static void generateConfiguration(int stage, float *position, float *velocity, float *& elasticConnectionsData_cpp, int *membraneData_cpp, int & numOfLiquidP, int & numOfElasticP, int & numOfBoundaryP, int & numOfElasticConnections, int & numOfMembranes, int * particleMembranesList_cpp);
	static void preLoadConfiguration();
	static void loadConfiguration( float * position, float  * velocity, float *& elasticConnectionsData_cpp, int & numOfLiquidP, int & numOfElasticP, int & numOfBoundaryP, int & numOfElasticConnections );
	static void loadConfigurationFromOneFile(float * position, float  * velocity, float *& elasticConnectionsData_cpp, int & numOfLiquidP, int & numOfElasticP, int & numOfBoundaryP, int & numOfElasticConnections);
	static void loadConfigurationToFile(float * position, float * connections=NULL, int * membranes=NULL, bool firstIteration = true);
	static void loadConfigurationFromFile(float *& position, float *& connections, int *& membranes, int iteration = 0);
	static void loadConfigurationFromFile_experemental(float *& position, float *& connections, int *& membranes, int iteration = 0);
	static void log_bufferf( const float * buffer, const int element_size,const int global_size, const char * fileName);
	static void log_bufferi( const int * buffer, const int element_size,const int global_size, const char * fileName);
	void watch_report(const char *str);
	double get_elapsedTime() { return elapsedTime; };
	void refreshTime();	
private:
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
