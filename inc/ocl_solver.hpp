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
#ifndef OCLSOLVER_HPP
#define OCLSOLVER_HPP

#if defined(_WIN32) || defined(_WIN64)
#pragma comment(lib, "opencl.lib") // opencl.lib
#endif

#if defined(__APPLE__) || defined(__MACOSX)

#include "OpenCL/cl.hpp"

#else
#include <CL/cl.hpp>
#endif

#include "isolver.h"
#include "ocl_const.h"
#include "ocl_struct.h"
#include "particle.h"
#include "sph_model.hpp"
#include "util/error.h"
#include "device.h"
#include <fstream>
#include <iostream>

namespace sibernetic {
    namespace solver {
        using std::cout;
        using std::endl;
        using std::shared_ptr;
        using sibernetic::model::particle;
        using sibernetic::model::partition;
        using sibernetic::model::sph_model;
        using sibernetic::ocl_error;
// OCL constans block
#define QUEUE_EACH_KERNEL 1

        const int LOCAL_NDRANGE_SIZE = 256;

        template<class T = float>
        class ocl_solver : public i_solver {
            typedef shared_ptr<sph_model<T>> model_ptr;

        public:
            ocl_solver(model_ptr &m, shared_ptr<device> d) : model(m), dev(d) {
                try {
                    this->initialize_ocl();
                } catch (ocl_error &ex) {
                    throw;
                }
            }

            // TODO rename method!!!
            virtual void init_model(const partition &p) {
                this->p = p;
                init_buffers();
                init_kernels();
                copy_buffer_to_device((void *) &(model->get_particles()[p.start]),
                                      b_particles, p.size() * sizeof(particle<T>));
            }

            ~ocl_solver() {}

            virtual void run_neighbour_search() {
              this->run_init_ext_particles();
            }

            virtual void run_physic() {}
            virtual void run() {
                run_neighbour_search();
                run_physic();
            }
        private:
            model_ptr model;
            partition p;
            shared_ptr<device> dev;
            std::string msg = dev->name + '\n';
            const std::string cl_program_file = "cl_code//sph_cl_code.cl";
            cl::Kernel k_init_ext_particles;
            cl::Kernel k_hash_particles;
            cl::Buffer b_particles;
            cl::Buffer b_ext_particles;
            cl::CommandQueue queue;
            cl::Program program;

            void init_buffers() {
                create_ocl_buffer("particles", b_particles, CL_MEM_READ_WRITE,
                                  p.size() * sizeof(particle<T>));
                create_ocl_buffer("ext_particles", b_ext_particles, CL_MEM_READ_WRITE,
                                  p.size() * sizeof(extend_particle));
            }

            void init_kernels() {
                create_ocl_kernel("k_init_ext_particles", k_init_ext_particles);
                create_ocl_kernel("k_hash_particles", k_hash_particles);
            }

            virtual void init_ext_particles() {}

            void initialize_ocl() {
                int err;
                queue = cl::CommandQueue(dev->context, dev->dev, 0, &err);
                if (err != CL_SUCCESS) {
                    throw ocl_error(msg + "Failed to create command queue");
                }
                std::ifstream file(cl_program_file);
                if (!file.is_open()) {
                    throw ocl_error(msg + "Could not open file with OpenCL program check "
                                          "input arguments oclsourcepath: " +
                                    cl_program_file);
                }
                std::string programSource(std::istreambuf_iterator<char>(file),
                                          (std::istreambuf_iterator<char>()));
                // TODO fix this to param
                if (0) {
                    programSource =
                            "#define _DOUBLE_PRECISSION\n" +
                            programSource; // not now it needs double extension check on device
                }
                cl::Program::Sources source(
                        1, std::make_pair(programSource.c_str(), programSource.length() + 1));
                program = cl::Program(dev->context, source);
#if defined(__APPLE__)
                err = program.build("-g -cl-opt-disable -I .");
#else
#if INTEL_OPENCL_DEBUG
                err = program.build(OPENCL_DEBUG_PROGRAM_PATH + "-g -cl-opt-disable -I .");
#else
                err = program.build("-I .");
#endif
#endif
                if (err != CL_SUCCESS) {
                    std::string compilationErrors;
                    compilationErrors = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev->dev);
                    msg += make_msg(msg, "Compilation failed: ", compilationErrors);
                    throw ocl_error(msg);
                }
                std::cout
                        << msg
                        << "OPENCL program was successfully build. Program file oclsourcepath: "
                        << cl_program_file << std::endl;
                return;
            }

            void create_ocl_buffer(const char *name, cl::Buffer &b,
                                   const cl_mem_flags flags, const int size) {
                int err;
                b = cl::Buffer(dev->context, flags, size, NULL, &err);
                if (err != CL_SUCCESS) {
                    std::string error_m =
                            make_msg("Buffer creation failed: ", name, " Error code is ", err);
                    throw ocl_error(error_m);
                }
            }

            void create_ocl_kernel(const char *name, cl::Kernel &k) {
                int err;
                k = cl::Kernel(program, name, &err);
                if (err != CL_SUCCESS) {
                    std::string error_m =
                            make_msg("Kernel creation failed: ", name, " Error code is ", err);
                    throw ocl_error(error_m);
                }
            }

            void copy_buffer_to_device(const void *host_b, cl::Buffer &ocl_b,
                                       const int size) {
                // Actualy we should check  size and type
                int err = queue.enqueueWriteBuffer(ocl_b, CL_TRUE, 0, size, host_b);
                if (err != CL_SUCCESS) {
                    std::string error_m =
                            make_msg("Copy buffer to device is failed error code is ", err);
                    throw ocl_error(error_m);
                }
                queue.finish();
            }

            void copy_buffer_from_device(void *host_b, const cl::Buffer &ocl_b,
                                         const int size) {
                // Actualy we should check  size and type
                int err = queue.enqueueReadBuffer(ocl_b, CL_TRUE, 0, size, host_b);
                if (err != CL_SUCCESS) {
                    std::string error_m =
                            make_msg("Copy buffer from device is failed error code is ", err);
                    throw ocl_error(error_m);
                }
                queue.finish();
            }

            int run_init_ext_particles() {
                this->kernel_runner(
                        this->k_init_ext_particles,
                        model->size(),
                        0,
                        this->b_ext_particles
                );
            }

            int run_hash_particles() {
                this->kernel_runner(
                        this->k_hash_particles,
                        model->size(),
                        0,
                        this->b_ext_particles
                );
            }

            template<typename U, typename... Args>
            int kernel_runner(
                    cl::Kernel &ker,
                    int dim,
                    int arg_pos,
                    U& arg,
                    Args... args) {
                ker.setArg(arg_pos, arg);
                return kernel_runner(ker, dim, ++arg_pos, args...);
            }

            template<typename U>
            int kernel_runner(cl::Kernel &ker, int dim, int arg_pos, U& arg) {
                ker.setArg(arg_pos, arg);
                int err = queue.enqueueNDRangeKernel(
                        k_init_ext_particles, cl::NullRange,
                        cl::NDRange(dim),
#if defined(__APPLE__)
                        cl::NullRange, nullptr, nullptr);
#else
                cl::NDRange(LOCAL_NDRANGE_SIZE), nullptr, nullptr);
#endif
#if QUEUE_EACH_KERNEL
                queue.finish();
#endif
                if (err != CL_SUCCESS) {
                    throw ocl_error(make_msg("An ERROR is appearing during work of kernel _runClearBuffers ", err));
                }
                return err;
            }

        };
    } // namespace solver
} // namespace sibernetic
#endif // OCLSOLVER_HPP