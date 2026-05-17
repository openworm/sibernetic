// cuda_common.h — shared includes, CUDA error-check macro, and host I/O
// helpers for the differentiable Sibernetic substrate.
//
// Every .cu in this directory includes this. The .cu files compile
// independently with `nvcc -rdc=true` and device-link at the final
// invocation; kernels defined in shaders.cu are launched from drivers
// living in ops_*.cu via the prototypes in shaders.h.
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstdint>

#include <cuda_runtime.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define CUDA_CHECK(call)                                                   \
    do {                                                                   \
        cudaError_t err__ = (call);                                        \
        if (err__ != cudaSuccess) {                                        \
            fprintf(stderr, "CUDA error %s at %s:%d: %s\n",                \
                    #call, __FILE__, __LINE__, cudaGetErrorString(err__)); \
            std::exit(1);                                                  \
        }                                                                  \
    } while (0)

// Host I/O helpers — fp32 binary slurp/dump used by every CLI driver.
// Defined in cuda_common.cu.
float *read_floats_or_die(const char *path, size_t n_floats);
void   write_floats_or_die(const char *path, const float *buf, size_t n_floats);
