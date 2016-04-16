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

/** owPhysicsFluidSimulator class contains
 *  realization of algorithms.
 */
class owPhysicsFluidSimulator
{
public:
	owPhysicsFluidSimulator(void);
	owPhysicsFluidSimulator(owHelper * helper, int argc, char ** argv);
	~owPhysicsFluidSimulator(void);
	/** Getter for position_cpp
	 *
	 *  Method doesn't need to request data from OpenCL device's memory ,
	 *  cause new information updating every time step in
	 *  owPhysicsFluidSimulator::simulationStep() method.
	 *
	 *  @return position_cpp
	 */
	float * getPosition_cpp() const { return position_cpp; };
	/** Getter for velocity_cpp buffer
	 *
	 *  When run this method information about new value of velocity
	 *  getting from OpenCL memory
	 *  it use owOpenCLSolver::read_velocity_buffer(...) method
	 *
	 *  @return velocity_cpp
	 */
	float * getvelocity_cpp() { ocl_solver->read_velocity_buffer(velocity_cpp,config); return velocity_cpp; };
	/** Getter for density_cpp buffer
	 *
	 *  When run this method information about new values of density
	 *  getting from OpenCL memory
	 *  it use owOpenCLSolver::read_density_buffer(...) method
	 *
	 *  @return density_cpp
	 */
	float * getDensity_cpp() { ocl_solver->read_density_buffer( density_cpp, config ); return density_cpp; };
	/** Getter for particleIndex_cpp buffer
	 *
	 *  When run this method information about new values of particleIndex
	 *  getting from OpenCL memory
	 *  it use owOpenCLSolver::read_particleIndex_buffer(...) method
	 *
	 *  @return particleIndex_cpp
	 */
	unsigned int * getParticleIndex_cpp() { ocl_solver->read_particleIndex_buffer( particleIndex_cpp, config ); return particleIndex_cpp; };
	/** Getter for elasticConnectionsData_cpp buffer
	 *
	 *  Method doesn't need to request data from OpenCL device's memory
	 *  information about elastic connection defines once when
	 *  initialization of configuration.
	 *
	 *  @return elasticConnectionsData_cpp
	 */
	float * getElasticConnectionsData_cpp() const { return elasticConnectionsData_cpp; };
	/** Getter for membraneData_cpp buffer
	 *
	 *  Method doesn't need to request data from OpenCL device's memory
	 *  information about elastic connection defines once when
	 *  initialization of configuration.
	 *
	 *  @return membraneData_cpp
	 */
	int   * getMembraneData_cpp() const { return membraneData_cpp; };
	/** Getter for muscle_activation_signal_cpp buffer
	 *
	 *
	 *  @return membraneData_cpp
	 */
	float * getMuscleAtcivationSignal() const { return muscle_activation_signal_cpp; }
	double  simulationStep(const bool load_to = false);
	/** Getter for config
	 *  @return config
	 */
	owConfigProperty * getConfig() const { return config; };
	/** Getter for iteration
	 *  @return iteration
	 */
	const int getIteration() const { return iterationCount; };
	void reset();
	void makeSnapshot();
private:
	owOpenCLSolver * ocl_solver;
	float * position_cpp;				// everywhere in the code %variableName%_cpp means that we create
	float * velocity_cpp;				// and initialize in 'ordinary' memory some data, which will be
	float * elasticConnectionsData_cpp; // copied later to OpenCL buffer %variableName%
	int	  * membraneData_cpp;
	int   * particleMembranesList_cpp;
	//Muscle contraction data buffer
	float * muscle_activation_signal_cpp;
	//Helper arrays for displaying information about density changes
	float * density_cpp;
	unsigned int * particleIndex_cpp;
	owConfigProperty * config;
	owHelper * helper;
	int iterationCount;
	void destroy();
};

#endif //OW_PHYSICS_SIMULATOR_H
