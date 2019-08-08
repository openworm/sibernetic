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
  size_t device_compute_unit_num; // criteria to sort devices
  size_t global_work_group_size; // criteria to sort devices
  size_t max_thread_count;
  float balance_coeff;
  bool operator<(const device &d1) {
    return  max_thread_count < d1.max_thread_count;
  }

  void show_info() {
    char c_buffer[100];
    cl_int result;
    size_t work_group_size, comp_unints_count;
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
    size_t work_group_size;
    result = dev.getInfo(CL_DEVICE_NAME, &c_buffer);
    if (result != CL_SUCCESS) {
    }
    name = c_buffer;
    dev.getInfo(CL_DEVICE_TYPE, &c_buffer);
    auto err = dev.getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &device_compute_unit_num);
    if(err != CL_SUCCESS || device_compute_unit_num > 1000){
        device_compute_unit_num = 1;
    }
    result = dev.getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &work_group_size);
    if(result != CL_SUCCESS){
        global_work_group_size = 256;
    } else {
        if(work_group_size > 1024) {
            global_work_group_size = 256;
        }
        else {
            global_work_group_size = work_group_size;
        }
    }
    max_thread_count = global_work_group_size * device_compute_unit_num;
  }
};

#endif