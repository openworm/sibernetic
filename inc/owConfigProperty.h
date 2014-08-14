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
/*
 * This file contains definition for struct configuration
 */
#ifndef OWCONFIGURATION_H_
#define OWCONFIGURATION_H_

#include "owOpenCLConstant.h"

struct owConfigProrerty{
	//This value defines boundary of box in which simulation is
	//Sizes of the box containing simulated 'world'
	//Sizes choice is realized this way because it should be proportional to smoothing radius h
public:
	float xmin;
	float xmax;
	float ymin;
	float ymax;
	float zmin;
	float zmax;
	int gridCellsX;
	int gridCellsY;
	int gridCellsZ;
	int gridCellCount;
	const int getParticleCount(){ return PARTICLE_COUNT; };
	void setParticleCount(int value){
		PARTICLE_COUNT = value;
		PARTICLE_COUNT_RoundedUp = ((( PARTICLE_COUNT - 1 ) / local_NDRange_size ) + 1 ) * local_NDRange_size;
	};
	const int getParticleCount_RoundUp(){ return PARTICLE_COUNT_RoundedUp; };
	void setDeviceType(int type){ preferable_device_type=type; };
	const int getDeviceType(){ return preferable_device_type; };
private:
	int PARTICLE_COUNT;
	int PARTICLE_COUNT_RoundedUp;
	int preferable_device_type;// 0-CPU, 1-GPU
};


#endif /* OWCONFIGURATION_H_ */
