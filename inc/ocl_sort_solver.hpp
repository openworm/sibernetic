#include <utility>

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
#ifndef OCL_SORT_SOLVER_HPP
#define OCL_SORT_SOLVER_HPP

#if defined(_WIN32) || defined(_WIN64)
#pragma comment(lib, "opencl.lib") // opencl.lib
#endif

#if defined(__APPLE__) || defined(__MACOSX)

#include "OpenCL/cl.hpp"

#else

#include <CL/cl.hpp>

#endif

#include "isort_solver.h"
#include "ocl_const.h"
#include "ocl_struct.h"
#include "particle.h"
#include "sph_model.hpp"
#include "util/error.h"
#include "util/json.hpp"
#include "device.h"
#include <fstream>
#include <iostream>
#include <unordered_map>

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

		template<class T = float>
		class ocl_sort_solver : public i_sort_solver {
			typedef shared_ptr<sph_model<T>> model_ptr;

		public:
			ocl_sort_solver(
					model_ptr &m,
					shared_ptr<device> d,
					size_t idx,
					LOGGING_MODE log_mode = LOGGING_MODE::NO):
				model(m),
				dev(std::move(d)),
				device_index(idx),
				log_mode(log_mode),
                GREATEST_NEAR_POW_OF_2(2)
			{
				try {
					this->initialize_ocl();
					result_index=std::vector<int>(model->size());
				} catch (ocl_error &ex) {
					throw;
				}
			}
            std::shared_ptr<device> get_device() override {
                return this->dev;
			}
			// TODO rename method!!!
			void init_model() {
                init_buffers();
                init_kernels();
                unsigned int n = model->size();
                if (n && !(n & (n - 1)))
                    GREATEST_NEAR_POW_OF_2 = n;

                while (GREATEST_NEAR_POW_OF_2 < n)
                    GREATEST_NEAR_POW_OF_2 <<= 1;
			}

			~ocl_sort_solver() override = default;
			void _debug_(){
				copy_buffer_from_device(&(model->get_particles()[0]), b_particles, model->size() * sizeof(particle<T>), 0);
				std::string big_s = "[";
				for(auto p: model->get_particles()){

					big_s += "{\"particle\": ";
					big_s += p.jsonify();
					big_s += "}";
					//break;
				}
				big_s += "]";
				std::ofstream debug_file("debug");
				debug_file << big_s;
				debug_file.close();
			}

			void sort() override {
                copy_buffer_to_device((void *) &(model->get_particles()[0]),
                                      b_particles, 0, model->size() * sizeof(particle<T>));
                run_init_index_array();
				run_sort();
			}

		private:
			model_ptr model;
			size_t device_index;
            shared_ptr<device> dev;
			std::string msg = dev->name + '\n';
			const std::string cl_program_file = "cl_code//cl_sort.cl";
			LOGGING_MODE log_mode;
            unsigned int GREATEST_NEAR_POW_OF_2;
			cl::Kernel k_sort;
			cl::Kernel k_fill_index_array;
			
			cl::Buffer b_particles;
            cl::Buffer b_index_array;
            cl::Buffer b_swap_index_array;
			cl::CommandQueue queue;
			cl::Program program;

            std::unordered_map<int, int> _result_map;
            std::vector<int> result_index;

            void init_buffers() {
				create_ocl_buffer("particles", b_particles, CL_MEM_READ_WRITE,
                                  model->size() * sizeof(particle<T>));
				create_ocl_buffer("index_array", b_index_array, CL_MEM_READ_WRITE,
				                  model->size() * sizeof(int));
                create_ocl_buffer("swap_index_array", b_swap_index_array, CL_MEM_READ_WRITE,
                                  model->size() * sizeof(int));

                copy_buffer_to_device((void *) &(model->get_particles()[0]),
                                      b_particles, 0, model->size() * sizeof(particle<T>));
			}

			void init_kernels() {
				create_ocl_kernel("k_sort", k_sort);
				create_ocl_kernel("k_fill_index_array", k_fill_index_array);
			}

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
			}

			void create_ocl_buffer(const char *name, cl::Buffer &b,
			                       const cl_mem_flags flags, const int size) {
				int err;
				b = cl::Buffer(dev->context, flags, size, nullptr, &err);
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

			void copy_buffer_to_device(
					const void *host_b,
					cl::Buffer &ocl_b,
			        const size_t offset,
			        const size_t size) {
				// Actually we should check  size and type
				int err = queue.enqueueWriteBuffer(ocl_b, CL_TRUE, offset, size, host_b);
				if (err != CL_SUCCESS) {
					std::string error_m =
							make_msg("Copy buffer to device is failed error code is ", err);
					throw ocl_error(error_m);
				}
				queue.finish();
			}

			void copy_buffer_from_device(void *host_b, const cl::Buffer &ocl_b,
			                             const size_t size, size_t offset) {
				// Actualy we should check  size and type
				int err = queue.enqueueReadBuffer(ocl_b, CL_TRUE, offset, size, host_b);
				if (err != CL_SUCCESS) {
					std::string error_m =
							make_msg("Copy buffer from device is failed error code is ", err);
					throw ocl_error(error_m);
				}
				queue.finish();
			}

			int run_init_index_array() {
				if(log_mode == LOGGING_MODE::FULL)
					std::cout << "run init index buffer --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_fill_index_array,
                        model->size(),
						0,
						this->b_index_array,
                        static_cast<unsigned int>(model->size())
				);
			}

			int run_sort() {
                if(log_mode == LOGGING_MODE::FULL)
                    std::cout << "run run sort --> " << dev->name << std::endl;
                bool need_swap = false;
//                for(auto p: model->get_particles()){
//                    std::cout << p.particle_id << "\t" << p.cell_id << std::endl;
//                }
				for(unsigned int step = 2; step <= GREATEST_NEAR_POW_OF_2; step<<=1){
                    auto _count = (model->size() % step == 0) ? model->size() / step : model->size() / step + 1;
				    if(need_swap) {
                        this->kernel_runner(
                                this->k_sort,
                                _count,
                                0,
                                this->b_particles,
                                this->b_swap_index_array,
                                this->b_index_array,
                                static_cast<unsigned int>(model->size()),
                                step,
                                static_cast<unsigned int>(_count)
                        );
                        need_swap = false;
                    } else {
                        this->kernel_runner(
                                this->k_sort,
                                _count,
                                0,
                                this->b_particles,
                                this->b_index_array,
                                this->b_swap_index_array,
                                static_cast<unsigned int>(model->size()),
                                step,
                                static_cast<unsigned int>(_count)
                        );
                        need_swap = true;
				    }
//				    if(step == 128) {
//                        std::cout << "================" << std::endl;
//                        shuffle_particles(need_swap);
//                        for (auto p: model->get_particles()) {
//
//                            std::cout << p.particle_id << "\t" << p.cell_id << std::endl;
//                        }
//                        std::cout << "================" << std::endl;
//                    }
				}
				shuffle_particles(need_swap);
			}
			void shuffle_particles(bool need_swap){
                if(need_swap){
                    copy_buffer_from_device(&(result_index[0]), b_swap_index_array, model->size() * sizeof(int), 0);
                } else {
                    copy_buffer_from_device(&(result_index[0]), b_index_array, model->size() * sizeof(int), 0);
                }
                for(int i =0;i<result_index.size();++i){
                    _result_map[result_index[i]] = i;
                }
                //return 1;
                particle<T> tmp;
                for(int i=0; i < result_index.size();++i){
                    if(result_index[i] != i && _result_map[i] != -1) {
                        auto from = i;
                        auto to = _result_map[i];
                        tmp = model->get_particles()[from];
                        while(to != -1){
                            std::swap<particle<T>>(model->get_particles()[to], tmp);
                            _result_map[from] = -1;
                            from = to;
                            to = _result_map[to];
                        }
                    }
                }
            }
			template<typename U, typename... Args>
			int kernel_runner(
					cl::Kernel &ker,
					unsigned int dim,
					unsigned int arg_pos,
					U &arg,
					Args... args) {
				ker.setArg(arg_pos, arg);
				return kernel_runner(ker, dim, ++arg_pos, args...);
			}

			template<typename U>
			int kernel_runner(cl::Kernel &ker, unsigned int dim, unsigned int arg_pos, U &arg) {
			    auto lk = this->dev->global_work_group_size;
				ker.setArg(arg_pos, arg);
				auto dim_round_up = (((dim - 1) / lk) + 1) * lk;
				int err = queue.enqueueNDRangeKernel(
						ker, cl::NullRange, cl::NDRange(dim_round_up),
#if defined(__APPLE__)
						cl::NullRange, nullptr, nullptr);
#else
						cl::NDRange(lk), nullptr, nullptr);
#endif
#if QUEUE_EACH_KERNEL
				queue.finish();
#endif
				if (err != CL_SUCCESS) {
					throw ocl_error(make_msg("An ERROR is appearing during work of kernel", err));
				}
				return err;
			}

		};
	} // namespace solver
} // namespace sibernetic
#endif // OCL_SORT_SOLVER_HPP