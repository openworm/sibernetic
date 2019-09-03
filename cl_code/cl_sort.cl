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
        uint total_size,
        uint step,
        uint block_count
){
	int id = get_global_id(0);
    if(id >= block_count) {
        return;
    }
    int start = id * step;
    int end = start + step - 1;
    if(end >= total_size ) {
        end = total_size - 1;
    }
    if (end - start == 0){
        result_index_array[end] = index_array[end];
        return;
    }
    int len = end - start + 1;
    int mid = start + (len >> 1);

    int first_sub_array_start = start;
    int first_sub_array_end = mid - 1;

    int second_sub_array_start;
    int second_sub_array_end = end;

    if(len % 2 == 0) {
        second_sub_array_start = mid;
    } else {
        second_sub_array_start = mid + 1;
        result_index_array[mid] = index_array[mid];
    }

    int i=first_sub_array_start, j=second_sub_array_start;
    bool copy_from_left = false, copy_from_right = false;
    int idx = start;
    printf("\nSTEP - %d, LEN - %d, mid - %d, first start - %d, first end - %d, second start - %d, second end - %d\n", step, len, mid, first_sub_array_start,first_sub_array_end, second_sub_array_start, second_sub_array_end);
    while(1){
        if(particles[i].cell_id > particles[j].cell_id ){
            result_index_array[idx++] = index_array[j++];
            printf("\n====1 idx %d val %d\n", idx - 1, result_index_array[idx-1]);
        } else {
            result_index_array[idx++] = index_array[i++];
            printf("\n====2 idx %d val %d\n", idx - 1, result_index_array[idx-1]);
        }
        if(i > first_sub_array_end) {
            copy_from_right = true;
            break;
        }else if(j > second_sub_array_end){
            copy_from_left = true;
            break;
        }
    }
    if(copy_from_left) {
        while(i <= first_sub_array_end){
            result_index_array[idx++] = index_array[i++];
        }
    } else if(copy_from_right) {
        while(j <= second_sub_array_end){
            result_index_array[idx++] = index_array[j++];
        }
    }
}

__kernel void k_fill_index_array(__global int * index_array, uint total_size){
    int id = get_global_id(0);
    if(id >= total_size){
        return;
	}
    index_array[id] = id;
}