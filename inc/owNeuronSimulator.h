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
/* 	owNeuronSimulator.h
 * 	C++ part of NEURON <-> sibernetic bridge
 * 	It will work if yo have NEURON simulator installed
 * 	on your OS also you need install sibernetic_neuron bridge
 * 	you can download it from github https://github.com/openworm/sibernetic_NEURON
 *  Created on: Jun 27, 2016
 *      Author: serg
 */

#ifndef INC_OWNEURONSIMULATOR_H_
#define INC_OWNEURONSIMULATOR_H_

#if defined(_WIN32) || defined (_WIN64)
  #include "C:/Python27/include/Python.h" // TODO make it optional
#else
  #include <Python.h>
#endif
#include <vector>
#include <iostream>
#include "owINeuronSimulator.h"

class owNeuronSimulator: public owINeuronSimulator {
public:
	owNeuronSimulator(int muscleNumber, float timeStep, const std::string & modelFileName, const std::string & simFileName = "main");
	std::vector<float> run();
	virtual ~owNeuronSimulator();
};

#endif /* INC_OWNEURONSIMULATOR_H_ */
