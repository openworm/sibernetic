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
struct device_ptr_less{
	bool operator()(const std::shared_ptr<device> &x, const std::shared_ptr<device> &y) const
	{ return *x < *y; }
};

typedef std::priority_queue<std::shared_ptr<device>, std::deque<std::shared_ptr<device>>, device_ptr_less> p_q;

p_q get_dev_queue();
#endif // OCLHELPER_H