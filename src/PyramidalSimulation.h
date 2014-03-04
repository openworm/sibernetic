#ifndef PYRAMIDALSIMULATION_H
#define PYRAMIDALSIMULATION_H
//#include "/usr/include/python2.7/Python.h"  //need to fix
//#define MS_NO_COREDLL
#if defined (__APPLE__)
  #include <Python.h>
#else
  #include "C:/Python27/include/Python.h"
#endif
#include <vector>

//#pragma comment( lib, "C:\\Python27\\libs\\python27.lib" )
//#pragma comment( lib, "C:/Python26/libs/python26.lib" )
//#pragma comment( lib, "C:/Python26/libs/python26.lib" )

using namespace std;



class PyramidalSimulation{

 private:
  PyObject *pName, *pModule, *pDict, *pFunc, *pValue, *pClass, *pInstance;
  vector<float> unpackPythonList(PyObject*);

 public:
  int setup();
  vector<float> run();

};

#endif
