#ifndef OCLHELPER_H
#define OCLHELPER_H
#include "device.h"
#include <memory>
#include <queue>
#if defined(__APPLE__) || defined(__MACOSX)
#include "../inc/OpenCL/cl.hpp"
//	#include <OpenCL/cl_d3d10.h>
#else
#include <CL/cl.hpp>
#endif
std::priority_queue<std::shared_ptr<device>> get_dev_queue();
#endif // OCLHELPER_H