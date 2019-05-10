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
			ocl_solver(model_ptr &m, shared_ptr<device> d, size_t idx):
				model(m),
				dev(std::move(d)),
				device_index(idx)
			{
				try {
					this->initialize_ocl();
				} catch (ocl_error &ex) {
					throw;
				}
			}

			// TODO rename method!!!
			void init_model(const partition &p) override {
				this->p = p;
				prev_part_size = p.total_size();
				init_buffers();
				init_kernels();
				copy_buffer_to_device((void *) &(model->get_particles()[p.ghost_start]),
				                      b_particles, 0, p.total_size() * sizeof(particle<T>));
			}

			~ocl_solver() override = default;

			void neighbour_search() override {
				run_init_ext_particles();
				run_hash_particles();
				run_clear_grid_hash();
				run_fill_particle_cell_hash();
				sync();
				run_neighbour_search();
			}

			void physic() override {
				int iter = 0;
				while(iter < sibernetic::model::PCI_ITER_COUNT) {
					run_compute_density();
					run_compute_forces_init_pressure();
					run_predict_positions();
					run_predict_density();
					run_correct_pressure();
					run_compute_pressure_force_acceleration();
					++iter;
				}
				run_integrate();
			}

			void sync() override {
				is_synchronizing = true;

				std::vector<particle<T>> tmp(p.size());
//				copy_buffer_from_device(
//						&(tmp[0]),
//						b_particles,
//						p.size() * sizeof(particle<T>),
//						p.offset() * sizeof(particle<T>));
//
//				for(int i=0;i<tmp.size();++i){
//					if (tmp[i].pos[0] != model->get_particles()[p.start + i].pos[0] ||
//					    tmp[i].pos[1] != model->get_particles()[p.start + i].pos[1] ||
//					    tmp[i].pos[2] != model->get_particles()[p.start + i].pos[2]){
//						std::cout << "Difference " << i << std::endl;
//						std::cout << "Difference x " << tmp[i].pos[0] - model->get_particles()[p.start + i].pos[0] << std::endl;
//						std::cout << "Difference y " << tmp[i].pos[1] - model->get_particles()[p.start + i].pos[1] << std::endl;
//						std::cout << "Difference z " << tmp[i].pos[2] - model->get_particles()[p.start + i].pos[2] << std::endl;
//						break;
//					}
//				}
				copy_buffer_from_device(
						&(model->get_particles()[p.start]),
						b_particles,
						p.size() * sizeof(particle<T>),
						p.offset() * sizeof(particle<T>));

				if(model->set_ready()){
					model->sync();
				} else {
					while(is_synchronizing);
				}

				copy_buffer_to_device(
						(void *) &(model->get_particles()[p.ghost_start]),
						b_particles,
						0,
						p.total_size() * sizeof(particle<T>));
				prev_part_size = p.total_size();


				/*BOTTOM IS IMPROVED VERSION*/
//				if(p.ghost_start == 0){
//					copy_buffer_to_device(
//							(void *) &(model->get_particles()[p.end]),
//				            b_particles,
//							p.end * sizeof(particle<T>),
//							(p.ghost_end - p.end) * sizeof(particle<T>));
//				} else if(p.ghost_end == model->size() - 1){
//					copy_buffer_to_device(
//							(void *) &(model->get_particles()[p.ghost_start]),
//							b_particles,
//							0,
//							(p.start - p.ghost_start) * sizeof(particle<T>));
//				} else {
//					copy_buffer_to_device(
//							(void *) &(model->get_particles()[p.end]),
//							b_particles,
//							p.end * sizeof(particle<T>),
//							(p.ghost_end - p.end) * sizeof(particle<T>));
//					copy_buffer_to_device(
//							(void *) &(model->get_particles()[p.ghost_start]),
//							b_particles,
//							0,
//							(p.start - p.ghost_start) * sizeof(particle<T>));
//				}
			}

		void run() override {
			int i = 0;
			while(i++ < 10000) {
				neighbour_search();
				physic();
			}
			std::cout << "?????????????? END OF WORK ????????????????" << dev->name << std::endl;
		}
		void unfreeze() override {
			is_synchronizing = false;
		}
		private:
			bool is_synchronizing;
			model_ptr model;
			size_t prev_part_size;
			size_t device_index;
			partition p;
			shared_ptr<device> dev;
			std::string msg = dev->name + '\n';

			const std::string cl_program_file = "cl_code//sph_cl_code.cl";
			cl::Kernel k_init_ext_particles;
			cl::Kernel k_hash_particles;
			cl::Kernel k_clear_grid_hash;
			cl::Kernel k_fill_particle_cell_hash;
			cl::Kernel k_neighbour_search;

			cl::Kernel k_compute_density;
			cl::Kernel k_compute_forces_init_pressure;
			cl::Kernel k_predict_positions;
			cl::Kernel ker_predict_density;
			cl::Kernel ker_correct_pressure;
			cl::Kernel ker_compute_pressure_force_acceleration;
			cl::Kernel k_integrate;

			cl::Buffer b_particles;
			cl::Buffer b_ext_particles;
			cl::Buffer b_grid_cell_id_list;
			cl::CommandQueue queue;
			cl::Program program;


			void init_buffers() {
				create_ocl_buffer("particles", b_particles, CL_MEM_READ_WRITE,
				                  p.total_size() * sizeof(particle<T>));
				create_ocl_buffer("ext_particles", b_ext_particles, CL_MEM_READ_WRITE,
				                  p.size() * sizeof(extend_particle));
				create_ocl_buffer("b_grid_cell_id_list", b_grid_cell_id_list, CL_MEM_READ_WRITE,
				                  p.total_cell_count() * sizeof(int));
			}

			void init_kernels() {
				create_ocl_kernel("k_init_ext_particles", k_init_ext_particles);
				create_ocl_kernel("k_hash_particles", k_hash_particles);
				create_ocl_kernel("k_clear_grid_hash", k_clear_grid_hash);
				create_ocl_kernel("k_fill_particle_cell_hash", k_fill_particle_cell_hash);
				create_ocl_kernel("k_neighbour_search", k_neighbour_search);

				create_ocl_kernel("k_compute_density", k_compute_density);
				create_ocl_kernel("k_compute_forces_init_pressure", k_compute_forces_init_pressure);
				create_ocl_kernel("k_predict_positions", k_predict_positions);
				create_ocl_kernel("ker_predict_density", ker_predict_density);
				create_ocl_kernel("ker_correct_pressure", ker_correct_pressure);
				create_ocl_kernel("ker_compute_pressure_force_acceleration", ker_compute_pressure_force_acceleration);
				create_ocl_kernel("k_integrate", k_integrate);
			}

			void init_ext_particles() override {}

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
							make_msg("Copy buffer to device is faqiled error code is ", err);
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

			int run_init_ext_particles() {
				std::cout << "run init_ext_particles --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_init_ext_particles,
						p.size(),
						0,
						this->b_ext_particles,
						p.size()
				);
			}

			int run_hash_particles() {
				std::cout << "run hash_particles --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_hash_particles,
						p.size(),
						0,
						this->b_particles,
						model->get_cell_num_x(),
						model->get_cell_num_y(),
						model->get_cell_num_z(),
						sibernetic::model::GRID_CELL_SIZE_INV,
						p.size(),
						p.offset(),
						p.limit()
				);
			}

			int run_clear_grid_hash() {
				std::cout << "run clear_grid_hash --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_clear_grid_hash,
						p.total_cell_count(),
						0,
						this->b_grid_cell_id_list,
						p.total_cell_count()
				);
			}

			int run_fill_particle_cell_hash() {
				std::cout << "run fill_particle_cell_hash --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_fill_particle_cell_hash,
						p.total_size(),
						0,
						this->b_grid_cell_id_list,
						this->b_particles,
						p.start_ghost_cell_id,
						p.total_size()
				);
			}

			int run_neighbour_search() {
				std::cout << "run neighbour_search --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_neighbour_search,
						p.size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						this->b_grid_cell_id_list,
						model->get_cell_num_x(),
						model->get_cell_num_y(),
						model->get_cell_num_z(),
						model->get_total_cell_num(),
						p.start_cell_id,
						sibernetic::model::H,
						sibernetic::model::GRID_CELL_SIZE,
						sibernetic::model::GRID_CELL_SIZE_INV,
						model->get_config()["simulation_scale"],
						model->get_config()["x_min"],
						model->get_config()["y_min"],
						model->get_config()["z_min"],
						p.size(),
						p.offset(),
						p.limit()
				);
			}

			int run_compute_density() {
				std::cout << "run compute_density --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_compute_density,
						p.size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["mass_mult_wpoly6_coefficient"],
						model->get_config()["h_scaled_2"],
						p.size(),
						p.offset(),
						p.limit()
				);
			}

			int run_compute_forces_init_pressure() {
				std::cout << "run compute_forces_init_pressure --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_compute_forces_init_pressure,
						p.size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["surf_tens_coeff"],
						model->get_config()["mass_mult_divgrad_viscosity_coefficient"],
						model->get_config()["h_scaled"],
						model->get_config()["gravity_x"],
						model->get_config()["gravity_y"],
						model->get_config()["gravity_z"],
						p.size(),
						p.offset(),
						p.limit()
				);
			}

			void run_predict_positions(){
				std::cout << "run predict_positions --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_predict_positions,
						p.size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["simulation_scale_inv"],
						model->get_config()["time_step"],
						sibernetic::model::R_0,
						p.size(),
						p.offset(),
						p.limit()
				);
			}

			void run_predict_density(){
				std::cout << "run predict_density --> " << dev->name << std::endl;
				this->kernel_runner(
						this->ker_predict_density,
						p.size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["mass_mult_wpoly6_coefficient"],
						sibernetic::model::H,
						model->get_config()["simulation_scale"],
						p.size(),
						p.offset(),
						p.limit()
				);
			}
			void run_correct_pressure(){
				std::cout << "run run_correct_pressure --> " << dev->name << std::endl;
				this->kernel_runner(
						this->ker_correct_pressure,
						p.size(),
						0,
						this->b_particles,
						sibernetic::model::DENSITY_WATER,
						model->get_config()["delta"],
						p.size(),
						p.offset(),
						p.limit()
				);
			}
			void run_compute_pressure_force_acceleration(){
				std::cout << "run run_compute_pressure_force_acceleration --> " << dev->name << std::endl;
				this->kernel_runner(
						this->ker_compute_pressure_force_acceleration,
						p.size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["mass_mult_grad_wspiky_coefficient"],
						model->get_config()["h_scaled"],
						model->get_config()["simulation_scale"],
						model->get_config()["delta"],
						sibernetic::model::DENSITY_WATER,
						p.size(),
						p.offset(),
						p.limit()
				);
			}
			void run_integrate(){
				std::cout << "run run_integrate --> " << dev->name << std::endl;
				this->kernel_runner(
						this->k_integrate,
						p.size(),
						0,
						this->b_ext_particles,
						this->b_particles,
						model->get_config()["simulation_scale_inv"],
						model->get_config()["time_step"],
						sibernetic::model::R_0,
						1,
						2,
						p.size(),
						p.offset(),
						p.limit()
				);
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
				ker.setArg(arg_pos, arg);
				auto dim_round_up = (((dim - 1) / LOCAL_NDRANGE_SIZE) + 1) * LOCAL_NDRANGE_SIZE;
				int err = queue.enqueueNDRangeKernel(
						ker, cl::NullRange, cl::NDRange(dim_round_up),
#if defined(__APPLE__)
						cl::NullRange, nullptr, nullptr);
#else
						cl::NDRange(LOCAL_NDRANGE_SIZE), nullptr, nullptr);
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
#endif // OCLSOLVER_HPP