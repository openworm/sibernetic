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
#include <algorithm>
#include <string>
#include <iterator>
#include <vector>

#include "PyramidalSimulation.h"

using namespace std;

int PyramidalSimulation::setup(){

  char python_module[] = "main_sim";
  char pyClass[] = "muscle_simulation";

  // Initialize the Python interpreter
  Py_Initialize();
  PyObject* pName;
  // Convert the file name to a Python string.
  pName = PyString_FromString(python_module);
  const char* s = PyString_AsString(pName);
  printf("[debug] pName = \"%s\"\n",s);
  const char* s2 = Py_GetPath();
  printf("[debug] PyPath = \"%s\"\n",s2);
  
  // Import the file as a Python module.
  pModule = PyImport_Import(pName);
  if( PyErr_Occurred() ) PyErr_Print();  

  // Build the name of a callable class
  if (pModule != NULL){
    pClass = PyObject_GetAttrString(pModule,pyClass);
	if( PyErr_Occurred() ) PyErr_Print();  
  }
  else {
    cout << "Module not loaded, have you set PYTHONPATH?" <<endl;
  }

  // Create an instance of the class
  if (PyCallable_Check(pClass))
    {
      pInstance = PyObject_CallObject(pClass, NULL);
	  if( PyErr_Occurred() ) PyErr_Print(); 
      cout << "Pyramidal simulation class loaded!"<<endl;
    }
  else {
    cout << "Pyramidal simulation class not callable! Try: export PYTHONPATH=$PYTHONPATH:./src"<<endl;
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
// pValue = PyObject_CallMethod(pInstance, "rrun
// un", NULL);
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
