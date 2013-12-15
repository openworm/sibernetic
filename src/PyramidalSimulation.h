#ifndef PYRAMIDALSIMULATION_H
#define PYRAMIDALSIMULATION_H
#include "/usr/include/python2.7/Python.h"  //need to fix
#include <vector>

//#pragma comment( lib, "C:\\Python27\\libs\\python27.lib" )

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
