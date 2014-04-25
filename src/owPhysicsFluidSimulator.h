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

#ifndef OW_PHYSICS_SIMULATOR_H
#define OW_PHYSICS_SIMULATOR_H

#include "owPhysicsConstant.h"
#include "owHelper.h"
#include "owOpenCLSolver.h"

class owPhysicsFluidSimulator
{
public:
	owPhysicsFluidSimulator(void);
	owPhysicsFluidSimulator(owHelper * helper);
	~owPhysicsFluidSimulator(void);
	float * getPosition_cpp() { return position_cpp; };
	float * getvelocity_cpp() { return velocity_cpp; };
	float * getDensity_cpp() { ocl_solver->read_density_buffer( density_cpp ); return density_cpp; };
	unsigned int * getParticleIndex_cpp() { ocl_solver->read_particleIndex_buffer( particleIndex_cpp ); return particleIndex_cpp; };
	//TODO helper functions delete after fix!!
	float * getElasticConnectionsData_cpp() { return elasticConnectionsData_cpp; };
	int   * getMembraneData_cpp() { return membraneData_cpp; };
	double  simulationStep(const bool load_to = false);

private:
	owOpenCLSolver * ocl_solver;
	float * position_cpp;				// everywhere in the code %variableName%_cpp means that we create 
	float * velocity_cpp;				// and initialize in 'ordinary' memory some data, which will be 
	float * elasticConnectionsData_cpp; // copied later to OpenCL buffer %variableName% 
	int	  * membraneData_cpp;
	int   * particleMembranesList_cpp;

	//Helper arrays
	float * density_cpp;
	unsigned int * particleIndex_cpp;
	float * acceleration_cpp;//TODO REMOVE after fixing
	owHelper * helper;
};

#endif //OW_PHYSICS_SIMULATOR_H
