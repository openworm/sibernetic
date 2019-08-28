/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2017 OpenWorm.
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


#ifdef cl_amd_printf
	#pragma OPENCL EXTENSION cl_amd_printf : enable
	#define PRINTF_ON // this comment because using printf leads to very slow work on Radeon r9 290x on my machine
    	              // don't know why
#elif defined(cl_intel_printf)
	#pragma OPENCL EXTENSION cl_intel_printf : enable
	#define PRINTF_ON
#endif

#ifdef _DOUBLE_PRECISION
#ifdef cl_khr_fp64
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable
	
#elif defined(cl_amd_fp64)
	#pragma OPENCL EXTENSION cl_amd_fp64 : enable
	//#define _DOUBLE_PRECISION
#else
	#error "Double precision floating point not supported by OpenCL implementation."
#endif
#endif
#if defined(_WIN32) || defined(_WIN64)
#include "inc\\ocl_struct.h"
#else
#include "inc/ocl_struct.h"
#endif

#define SQRT( x ) native_sqrt( x )
#define DIVIDE( a, b ) native_divide( a, b )
#define DOT( a, b ) dot( a, b )

#define NO_PARTICLE_ID -1

#define LIQUID_PARTICLE   1
#define ELASTIC_PARTICLE  2
#define BOUNDARY_PARTICLE 3

typedef struct particle{
#ifdef _DOUBLE_PRECISION
	double4 pos;
	double4 pos_n_1;
	double4 vel;
	double4 acceleration;
	double4 acceleration_n_1;
	float4 acceleration_n_0_5;
#else
	float4 pos;
	float4 pos_n_1;
	float4 vel;
	float4 acceleration;
	float4 acceleration_n_1;
	float4 acceleration_n_0_5;
#endif
	int type_;
	int cell_id;
	int particle_id;
#ifdef _DOUBLE_PRECISION
	double density;
	double pressure;
	double viscosity;
	double mass;
#else
	float density;
	float pressure;
	float viscosity;
	float mass;
#endif
} particle;


/** Just for test
*/
__kernel void k_sort(
        __global struct particle* particles,
        __global int * index_array,
        __global int * result_index_array,
        int step,
        int total_size
){
	int id = get_global_id(0);
    if(id > total_size) {
        return;
    }


}

__kernel void k_fill_index_array(
        __global int * index_array,
        int total_size
){
    int id = get_global_id(0);
    if(id >= total_size){
        return;
	}
    index_array[id] = id;
}