#ifndef __SPH__
#define __SPH__

#define PARTICLE_COUNT ( 32 * 1024 )
#define NEIGHBOR_COUNT 32

#ifndef M_PI
#define M_PI 3.1415927f
#endif

#define XMIN 0
#define XMAX 100
#define YMIN 0
#define YMAX 40
#define ZMIN 0
#define ZMAX 40


#endif

/*
OpenCL Hello World 

#define __CL_ENABLE_EXCEPTIONS
#define __NO_STD_VECTOR
#define __NO_STD_STRING

#include <CL/cl.hpp>
#include <cstdio>
#include <cstdlib>
#include <iostream>

const char * helloStr = "__kernel void hello(void) { }\n";

int main(void) {
   try {
      cl::Context context(CL_DEVICE_TYPE_GPU, 0, NULL, NULL, &err);
      cl::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
      cl::CommandQueue queue(context, devices[0], 0, &err);
      cl::Program::Sources source(1, std::make_pair(helloStr,strlen(helloStr)));
      cl::Program program_ = cl::Program(context, source);
      program_.build(devices);
      cl::Kernel kernel(program_, "hello", &err);
      cl::KernelFunctor func = kernel.bind(queue, cl::NDRange(4, 4), cl::NDRange(2, 2));
      func().wait();
   } catch (cl::Error err) {
      std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
   }
   return EXIT_SUCCESS;
}

*/