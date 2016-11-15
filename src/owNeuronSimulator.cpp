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
/*	owNeuronSimulator.cpp
 *	Implementation of owNeuronSimulator
 *  Created on: Jun 27, 2016
 *      Author: serg
 */
#include <stdexcept>
#include "owNeuronSimulator.h"

/** Constructor method for owNeuronSimulator.
 *
 *  @param muscleNumber
 *  Number of muscles in simulation
 *  @param time_step
 *  pass into NEURON for synchrony work
 *  @param simFileName
 *  File name with needed function default it's main.py
 *  from siberneti_NEURON and you don't need to change it actually
 *  You only need add path to sibenetic_NEURON folder into PYTHONPATH
 */
owNeuronSimulator::owNeuronSimulator(int muscleNumber,float timeStep, const std::string & modelFileName, const std::string & simFileName) {
	  // Initialize the Python interpreter
	  Py_Initialize();
	  PyObject* pName;
	  // Convert the file name to a Python string.
	  pName = PyString_FromString(simFileName.c_str());
	  const char* s = PyString_AsString(pName);
	  printf("[debug] pName = \"%s\"\n",s);
	  const char* s2 = Py_GetPath();
	  printf("[debug] PyPath = \"%s\"\n",s2);
	  // Import the file as a Python module.
	  pModule = PyImport_Import(pName);
	  if( PyErr_Occurred() ) PyErr_Print();
	  if (pModule == NULL){
	    throw std::runtime_error("Module not loaded, have you set PYTHONPATH? If yes just check have you added path to sibenetic_NEURON");
	  }
	  PyObject* sectionNames = PyList_New(muscleNumber); // size should be equal to number of muscles in model
	  for(int i=0; i<muscleNumber; ++i)
		  PyList_SetItem(sectionNames, i, PyString_FromString("SMDDR_mus"));    // there you actually need put section names (from you'r .hoc file)
	                                                                            // file with NEURON model
                                                                                // from which you want read info about signal (Voltage)
	                                                                            // Now it's hadrcoded for only one section SMDDR_mus it's a muscle section
	                                                                            // from model file in sibernetic_NEURON you can find it in folder
                                                                                // (path/to/sibewnretic_NEURON/model/c.elegans/ria_.hoc)
	  PyObject* dt = PyFloat_FromDouble(timeStep);                              // Create time step argument
	  PyObject* fileName = PyString_FromString(modelFileName.c_str());          // Create model file name argument
	  PyObject * args = PyTuple_Pack(3,dt,fileName,sectionNames);               // Create tuple of arguments for initialization
	  PyObject* myFunction = PyObject_GetAttrString(pModule,(char*)"run_init");
	  PyObject_CallObject(myFunction, args);                                    // Run initialization method
	  Py_DECREF(sectionNames);
	  Py_DECREF(dt);
	  Py_DECREF(fileName);
	  Py_DECREF(args);
}

/** Run one step of NEURON simulation and read the data
 *  @return vector with neuron activity
 */
std::vector<float> owNeuronSimulator::run(){
	// Call a method
	PyObject* myFunction = PyObject_GetAttrString(pModule,(char*)"run_sim_one_step");
	pValue = PyObject_CallObject(myFunction, NULL);
	//pValue = PyObject_CallMethod(nrn_sim, const_cast<char *>("one_step"), NULL);
	if(PyList_Check(pValue)){
	  std::vector<float> value_array;
	  value_array = unpackPythonList(pValue);
	  return value_array;

	}
	else {
	  std::vector<float> single_element_array(0);
	  single_element_array[0] = (float)PyFloat_AsDouble(pValue);
	  return single_element_array;
	}
}

owNeuronSimulator::~owNeuronSimulator() {
	// TODO Auto-generated destructor stub
}

