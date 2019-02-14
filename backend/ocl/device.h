#ifndef DEVICE_H
#define DEVICE_H
#include <iostream>
#include <string>
#if defined(__APPLE__) || defined(__MACOSX)
#include "OpenCL/cl.hpp"
#else
#include <CL/cl.hpp>
#endif
enum type { CPU, GPU };
struct device {
  device(cl::Device &d, int p_id, int id)
      : platform_id(p_id), dev_id(id), dev(d), context(dev) {
    init_params();
  }
  type t;
  int platform_id;
  int dev_id;
  int is_busy;
  std::string name;
  cl::Device dev;
  cl::Context context;
  size_t device_coumpute_unit_num; // criteria to sort devices
  bool operator<(const device &d1) {
    return device_coumpute_unit_num < d1.device_coumpute_unit_num;
  }
  void show_info() {
    char c_buffer[100];
    cl_int result;
    int work_group_size, comp_unints_count;
    // Print some information about chosen platform
    result = dev.getInfo(CL_DEVICE_NAME, &c_buffer); // CL_INVALID_VALUE = -30;
    if (result == CL_SUCCESS) {
      std::cout << "CL_CONTEXT_PLATFORM [" << platform_id
                << "]: CL_DEVICE_NAME [" << dev_id << "]:\t" << c_buffer << "\n"
                << std::endl;
    }
    result = dev.getInfo(CL_DEVICE_TYPE, &c_buffer);
    if (result == CL_SUCCESS) {
      std::cout << "CL_CONTEXT_PLATFORM [" << platform_id
                << "]: CL_DEVICE_TYPE [" << dev_id << "]:\t"
                << (((int)c_buffer[0] == CL_DEVICE_TYPE_CPU) ? "CPU" : "GPU")
                << std::endl;
    }
    result = dev.getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &work_group_size);
    if (result == CL_SUCCESS) {
      std::cout << "CL_CONTEXT_PLATFORM [" << platform_id
                << "]: CL_DEVICE_MAX_WORK_GROUP_SIZE [" << dev_id << "]: \t"
                << work_group_size << std::endl;
    }
    result = dev.getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &comp_unints_count);
    if (result == CL_SUCCESS) {
      std::cout << "CL_CONTEXT_PLATFORM [" << platform_id
                << "]: CL_DEVICE_MAX_COMPUTE_UNITS [" << dev_id << "]: \t"
                << comp_unints_count << std::endl;
    }
  }

private:
  void init_params() {
    char c_buffer[100];
    cl_int result;
    result = dev.getInfo(CL_DEVICE_NAME, &c_buffer);
    if (result != CL_SUCCESS) {
    }
    name = c_buffer;
    result = dev.getInfo(CL_DEVICE_TYPE, &c_buffer);
    t = ((int)c_buffer[0] == CL_DEVICE_TYPE_CPU) ? CPU : GPU;
    result =
        dev.getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &device_coumpute_unit_num);
  }
};

#endif