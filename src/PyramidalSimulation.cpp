// to compile and run (temp notes with hardcoded paths)
// run the following commands from inside curdir:
//
// $ export PYTHONPATH="/home/mike/dev/cpp_pyramidal_integration/"
//OR
//export PYTHONPATH=$PYTHONPATH:/home/mike/dev/muscle_model/pyramidal_implementation/
// $ g++ main.cpp -l python2.7 -o sim -I /usr/include/python2.7/
// $ ./sim

#include <Python.h>
#include <iostream>
#include "PyramidalSimulation.h"
#include <algorithm>
#include <string>
#include <iterator>
#include <vector>


using namespace std;

int PyramidalSimulation::setup(){

  char python_module[] = "main_sim";
  char pyClass[] = "muscle_simulation";

  Py_Initialize();
  PyObject* pName;
  pName = PyString_FromString(python_module);
  pModule = PyImport_Import(pName);

  //printf("%s",PyGETENV("PYTHONPATH"));

  //    // Build the name of a callable class
  if (pModule != NULL){
    pClass = PyObject_GetAttrString(pModule,pyClass);
  }
  else {
    cout << "Module not loaded, have you set PYTHONPATH?" <<endl;
  }

  // Create an instance of the class
  if (PyCallable_Check(pClass))
    {
      pInstance = PyObject_CallObject(pClass, NULL);
      cout << "Pyramidal simulation class loaded!"<<endl;
    }
  else {
    cout << "Pyramidal simulation class not callable!"<<endl;
  }

  return 0;
};

vector<float> PyramidalSimulation::unpackPythonList(PyObject* pValue){

	Py_ssize_t size = PyList_Size(pValue);
	vector<float> test(96); //needs to change!
	printf("====\n");
	for (Py_ssize_t i = 0; i < size; i++) {
		float value;
		value = PyFloat_AsDouble(PyList_GetItem(pValue, i));
		test[i]= value;
	}

	return test;
};

vector<float> PyramidalSimulation::run(){
  // Call a method of the class
//	  pValue = PyObject_CallMethod(pInstance, "rrun
//	  un", NULL);
	printf("!!!checkpoint001!!!\n");
  pValue = PyObject_CallMethod(pInstance, "run", NULL);
printf("!!!checkpoint002!!!\n");
  if(PyList_Check(pValue)){
	   printf("!!!checkpoint003.1!!!\n");
	  vector<float> value_array;
	  value_array = PyramidalSimulation::unpackPythonList(pValue);
	  printf("!!!checkpoint003!!!\n");
	  return value_array;

  }


  else {
	  printf("!!!checkpoint004.1!!!\n");
	  vector<float> single_element_array(0);
	  single_element_array[0] = PyFloat_AsDouble(pValue);
	  printf("!!!checkpoint004!!!\n");
	  return single_element_array;
  }
   printf("!!!checkpoint005!!!\n");
};
