/*
 *                  RADIXSORT.CL
 *
 * "radixsort.cl" is a compendium of functions for the openCL
 * implementation of the Radix Sort algorithm.
 *
 * 2017 Project for the "Facultad de Ciencias Exactas, Ingenieria
 * y Agrimensura" (FCEIA), Rosario, Santa Fe, Argentina.
 *
 * Implementation by Paoloni Gianfranco and Soncini Nicolas.
 */
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
#endif

#define SQRT( x ) native_sqrt( x )
#define DIVIDE( a, b ) native_divide( a, b )
#define DOT( a, b ) dot( a, b )

#include "inc/radixsort.h"
#define NULL 0

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

__kernel void prepare(__global int* cell_id_list,
                      __global int* array,
                      __global struct particle* in_particles,
                      const uint real_size,
                      const uint size)
{
    uint g_id = (uint) get_global_id(0);
    if(g_id >= size){
        return;
    }
    if(g_id < real_size) {
        cell_id_list[g_id] = in_particles[g_id].cell_id;
        array[g_id] = g_id;
    } else {
        cell_id_list[g_id] = 1000000000;
        array[g_id] = g_id;
    }
}

/** COUNT KERNEL **/

__kernel void count(const __global int* input,
                    __global int* output,
                    __local int* local_histo,
                    const int pass,
                    const int nkeys,
                    const __global int * particles)
{
    uint g_id = (uint) get_global_id(0);
    uint l_id = (uint) get_local_id(0);
    uint l_size = (uint) get_local_size(0);

    uint group_id = (uint) get_group_id(0);
    uint n_groups = (uint) get_num_groups(0);

    //Set the buckets of each item to 0
    int i;
    for(i = 0; i < BUCK; i++) {
        local_histo[i * l_size + l_id] = 0;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    //Calculate elements to process per item
    int size = (nkeys / n_groups) / l_size;

    //Calculate where to start on the global array
    int start = g_id * size;
    for(i = 0; i < size; i++) {
        int key = particles[input[i + start]];
        //Extract the corresponding radix of the key
        key = ((key >> (pass * RADIX)) & (BUCK - 1));
        //Count the ocurrences in the corresponding bucket
        ++local_histo[key * l_size + l_id];
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    for(i = 0; i < BUCK; i++) {
        //"from" references the local buckets
        int from = i * l_size + l_id;
        //"to" maps to the global buckets
        int to = i * n_groups + group_id;
        //Map the local data to its global position
        output[l_size * to + l_id] = local_histo[from];
    }

    barrier(CLK_GLOBAL_MEM_FENCE);
}

/** SCAN KERNEL **/
__kernel void scan(__global int* input,
                   __global int* output,
                   __local int* local_scan,
                   __global int* block_sum)
{
    uint g_id = (uint) get_global_id(0);
    uint l_id = (uint) get_local_id(0);
    uint l_size = (uint) get_local_size(0);

    uint group_id = (uint) get_group_id(0);
    uint n_groups = (uint) get_num_groups(0);

    //Store data from global to local memory to operate
    local_scan[2 * l_id] = input[2 * g_id];
    local_scan[2 * l_id + 1] = input[2 * g_id + 1];

    //UP SWEEP
    int d, offset = 1;
    for(d = l_size; d > 0; d >>= 1){
        barrier(CLK_LOCAL_MEM_FENCE);
        if(l_id < d) {
            int a = offset * (2 * l_id + 1) - 1;
            int b = offset * (2 * l_id + 2) - 1;
            local_scan[b] += local_scan[a];
        }
        offset *= 2;
    }

    if (l_id == 0) {
        //Store the full sum on last item
        if(block_sum != NULL){
            block_sum[group_id] = local_scan[l_size * 2 - 1];
        }

        //Clear the last element
        local_scan[l_size * 2 - 1] = 0;
    }

    //DOWN SWEEP
    for(d = 1; d < (l_size*2); d *= 2) {
        offset >>= 1;
        barrier(CLK_LOCAL_MEM_FENCE);
        if(l_id < d) {
            int a = offset * (2 * l_id + 1) - 1;
            int b = offset * (2 * l_id + 2) - 1;
            int tmp = local_scan[a];
            local_scan[a] = local_scan[b];
            local_scan[b] += tmp;
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    //Write results from Local to Global memory
    output[2 * g_id]     = local_scan[2 * l_id];
    output[2 * g_id + 1] = local_scan[2 * l_id + 1];
}


/** COALESCE KERNEL **/
__kernel void coalesce(__global int* scan,
                       __global int* block_sums)
{

    uint g_id = (uint) get_global_id(0);
    uint group_id = (uint) get_group_id(0);

    int b = block_sums[group_id];

    //TODO: Probar pasar a memoria local
    scan[2 * g_id] += b;
    scan[2 * g_id + 1] += b;

    barrier(CLK_GLOBAL_MEM_FENCE);
}



/** REORDER KERNEL **/
__kernel void reorder(__global int* array,
                      __global int* histo,
                      __global int* output,
                      const int pass,
                      const int nkeys,
                      __local int* local_histo,
                      const __global int * particles)
{
    uint g_id = (uint) get_global_id(0);
    uint l_id = (uint) get_local_id(0);
    uint l_size = (uint) get_local_size(0);

    uint group_id = (uint) get_group_id(0);
    uint n_groups = (uint) get_num_groups(0);

    //Bring histo to local memory
    int i;
    for(i = 0; i < BUCK; i++){
        int to = i * n_groups + group_id;
        local_histo[i * l_size + l_id] =
                histo[l_size * to + l_id];
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    //Write to global memory in order
    int size = (nkeys / n_groups) / l_size;
    int start = g_id * size;
    for(i = 0; i < size; i++){
        int item = array[i + start];
        int key = particles[item];
        key = (key >> (pass * RADIX)) & (BUCK - 1);
        int pos = local_histo[key * l_size + l_id];
        ++local_histo[key * l_size + l_id];

        output[pos] = item;
    }

    barrier(CLK_GLOBAL_MEM_FENCE);
}

__kernel void shuffle_particles(__global struct particle* in_particles,
                                __global struct particle* out_particles,
                                const int size,
                                __global int* index_buffer)
{
    uint g_id = (uint) get_global_id(0);
    if(g_id >= size){
        return;
    }

    out_particles[g_id] = in_particles[index_buffer[g_id]];
}