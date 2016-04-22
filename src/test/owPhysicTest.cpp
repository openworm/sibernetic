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
 * owPhycisTest.cpp
 *
 *  Created on: Apr 29, 2014
 *      Author: Sergey Khayrulin
 *      email: s.khayrulin@gmail.com
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>        // std::abs

#include "owPhysicTest.h"


float calcPotentialEnergy( owConfigProperty *, float * );
float calcKineticEnergy( owConfigProperty *, float *, float * );
float get_len( float * );

float gravity = 9.81f;

/*******************************************************
 * CONSERVATION ENERGY TEST DESCRIPTION
 * TOTAL ENERGY OF SYSTEM SHOULD BE CONSTANT ALL TIME
 * E = mv^2/2 + mgh
 * *****************************************************/
void test_energy_conservation(int argc, char **argv){
	//owHelper::path = "./configuration/test/"; TODO FIX it
	owHelper * helper = new owHelper();
	owPhysicsFluidSimulator * fluid_simulation = new owPhysicsFluidSimulator(helper, argc, argv);
	float total_energy = 0.f;
	float kinetic_energy = 0.f;
	float potential_energy = 0.f;
	float * p_buffer;
	float * v_buffer;
	std::vector<float> energy_evolution_total;
	std::vector<float> energy_evolution_potential;
	std::vector<float> energy_evolution_kinetic;
	int counter = 0;
	std::cout << "===================" << "CONSERVATION ENERGY TEST START" << "========================" << std::endl;
	while(1){
		p_buffer = fluid_simulation->getPosition_cpp();
		v_buffer = fluid_simulation->getvelocity_cpp();
		potential_energy = calcPotentialEnergy(fluid_simulation->getConfig(),p_buffer);
		kinetic_energy = calcKineticEnergy(fluid_simulation->getConfig(),v_buffer,p_buffer);
		total_energy = kinetic_energy + potential_energy;
		energy_evolution_total.push_back(total_energy);
		energy_evolution_kinetic.push_back(kinetic_energy);
		energy_evolution_potential.push_back(potential_energy);
		fluid_simulation->simulationStep();
		if(counter == 5000)
			break;
		counter++;
	}
	owHelper::log_buffer(&energy_evolution_total[0], 1, energy_evolution_total.size(), "./logs/total_energy_distrib.txt");
	owHelper::log_buffer(&energy_evolution_kinetic[0], 1, energy_evolution_kinetic.size(), "./logs/kinetic_energy_distrib.txt");
	owHelper::log_buffer(&energy_evolution_potential[0], 1, energy_evolution_potential.size(), "./logs/potential_energy_distrib.txt");
	std::cout << "===================" << "CONSERVATION ENERGY TEST END  " << "========================" << std::endl;
	delete fluid_simulation;
	delete p_buffer;
	delete v_buffer;
}
float calcPotentialEnergy(owConfigProperty * config, float * p_buffer){
	float e = 0.f;
	float l = 0.f;
	for(int i=0;i<config->getParticleCount();i++){
		if((int)(p_buffer[4 * i + 3]) != BOUNDARY_PARTICLE){
			l = (p_buffer[4 * i + 1] <= r0) ? 0.f : p_buffer[4 * i + 1] * simulationScale;
			e += std::abs(l) * gravity * mass; //Y - coordinate is a h
		}
	}
	return e;
}
float calcKineticEnergy(owConfigProperty * config, float * v_buffer, float * p_buffer){
	float e = 0.f;
	for(int i=0;i < config->getParticleCount();i++){
		if((int)(p_buffer[4 * i + 3]) != BOUNDARY_PARTICLE){
			e += mass * pow(get_len(v_buffer + 4 * i + 0),2.0f)/2.0f; //
		}
	}
	return e;
}
float get_dist(float * from, float * to){
	return sqrt ( pow(to[0] - from[0], 2.f) + pow(to[1] - from[1], 2.f) + pow(to[2] - from[2], 2.f) + pow(to[3] - from[3], 2.f) );
}
float get_len(float * v){
	return sqrt( pow(v[0],2.0f) + pow(v[1],2.0f) + pow(v[2],2.0f));
}
