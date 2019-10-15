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
#ifndef SPHMODEL_HPP
#define SPHMODEL_HPP

#include "particle.h"
#include "partition.h"
#include "util/error.h"
#include "util/abstract_reader.h"
#include "util/json_reader.hpp"
#include "abstract_model.hpp"
#include "isolver.h"
#include "isort_solver.h"
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <bitset>
#include <time.h>


#include <cstdlib>

namespace sibernetic {
	namespace model {
		/* const block end */
		template<class T = float, class container = std::vector<particle<T>>>
		class sph_model : public particle_model<T> {
			typedef std::map<std::string, T> sph_config;

		public:
			sph_model(const std::string &config_file, abstract_reader<T> *serializer = new json_reader<T>())
					: serializer(serializer), ready_to_sync(0), iteration(0) {
				config = {
					{"particles",                   T()},
					{"x_max",                       T()},
					{"x_min",                       T()},
					{"y_max",                       T()},
					{"y_min",                       T()},
					{"z_max",                       T()},
					{"z_min",                       T()},
					{"mass",                        T()},
					{"time_step",                   T()},
					{"simulation_scale",            T()},
					{"rho0",                        T()},
					{"time_step",                   T()},
					{"h_scaled_2",                  T()},
				};
				this->serializer->serialize(config_file, this);
				init_vars();
				for (particle<T> &p: particles) {
					this->calc_grid_id(p);
				}
				//arrange_particles();
				iteration = 0;
				clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
				t0 = t1;
				std::cout << "Model was loaded: " << particles.size() << " particles."
				          << std::endl;
			}

			const sph_config &get_config() const { return config; }

			const container &get_particles() const { return particles; }

			container &get_particles() { return particles; }

			const particle<T> &get_particle(const int index) const {
				return this->particles.at(index);
			}

			particle<T> &get_particle(const int index) {
				return this->particles.at(index);
			}

			void set_particle(int index, const particle<T> &p) {
				this->particles.at(index) = p;
			}

			void push_back(const particle<T> &p) {
				this->particles.push_back(p);
			}

			size_t size() const { return particles.size(); }

			/** Make partition for device
			 */
			void make_partition(size_t dev_count) {
				next_partition = 0;
				if (dev_count == 1) {
					//partitions.push_back(partition{0, static_cast<size_t>(size() - 1)});
					push_partition(0, size() - 1);
					return;
				}
				auto part_size = static_cast<size_t>(size() / dev_count);
				size_t start = 0;
				size_t end = 1;
				for (size_t i = 0; i < dev_count; ++i) {
					part_size = static_cast<size_t>(particles.size() * balance_coeff[i]);
					if (i == dev_count - 1)
						push_partition(start, static_cast<size_t>(size() - 1));
					else {
						if (particles[start + part_size].cell_id != particles[start + part_size + 1].cell_id) {
							push_partition(start, start + part_size);
							start += part_size + 1;
						} else {
							end = start + part_size;
							while (particles[end].cell_id == particles[end + 1].cell_id) {
								if (end + 1 > size()) {
									break;
								}
								++end;
							}
							if ((particles[end].cell_id + 1) % (cell_num_y * cell_num_z) != 0) {
								auto rest = (particles[end].cell_id + 1) % (cell_num_y * cell_num_z);
								auto last_cell = particles[end].cell_id + (cell_num_y * cell_num_z - rest) + 1;
								while (particles[end].cell_id != last_cell) {
									++end;
								}
								push_partition(start, end - 1);
								start = end;
							} else {
								push_partition(start, end);
								start = end + 1;
							}
						}
					}
				}
			}

			void set_balance_vector(const std::vector<float> &balance_coeff){
			    this->balance_coeff = balance_coeff;
			}

			partition &get_next_partition() {
				++next_partition;
				return partitions[next_partition - 1];
			}

			sph_config &get_config() {
				return this->config;
			}

			void sync() {
				arrange_particles();
//				for(size_t p_index=0; p_index< partitions.size(); ++p_index){
//					if(p_index == 0){
//						sync_right_segment(p_index, 0);
//					} else {
//						partitions[p_index].start = partitions[p_index - 1].end + 1;
//						sync_right_segment(p_index, partitions[p_index].start);
//						if(p_index == partitions.size() - 1)
//							sync_left_segment(p_index, partitions[p_index].end);
//					}
//				}

                update_partition_distrib();

				clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
				time_t sec = t2.tv_sec - t1.tv_sec;
				long nsec;
				if (t2.tv_nsec >= t1.tv_nsec) {
					nsec = t2.tv_nsec - t1.tv_nsec;
				} else {
					nsec = 1000000000 - (t1.tv_nsec - t2.tv_nsec);
					sec -= 1;
				}

				t1 = t2;
				elapsedTime = (float)(t2.tv_sec - t0.tv_sec) * 1000.f +
				              (float)(t2.tv_nsec - t0.tv_nsec) / 1000000.f;

				std::cout << "========New Iteration====== " << iteration << "====time=== " << (float)sec * 1000.f + (float)nsec / 1000000.f << std::endl;
				++iteration;
				ready_to_sync = 0;
				for(auto it: *solvers){
					it->unfreeze();
				}
			}

			void sync_right_segment(size_t p_index, size_t particle_start){
				for(size_t i = particle_start; i < particles.size(); ++i){
					if(particles[i].cell_id == partitions[p_index].end_cell_id
					   && particles[i + 1].cell_id != partitions[p_index].end_cell_id){
						partitions[p_index].end = i;
					} else if(particles[i].cell_id == partitions[p_index].end_ghost_cell_id
					          && particles[i + 1].cell_id != partitions[p_index].end_ghost_cell_id){
						partitions[p_index].ghost_end = i;
						break;
					}
				}
			}

			void sync_left_segment(size_t p_index, size_t particle_start){
				for(size_t i = particle_start; i > 0; --i){
					if(particles[i].cell_id == partitions[p_index].start_ghost_cell_id
					   && particles[i - 1].cell_id != partitions[p_index].start_ghost_cell_id){
						partitions[p_index].ghost_start = i;
						break;
					}
				}
			}

            bool set_ready() {
				std::lock_guard<std::mutex> guard(sync_mutex);
				++ready_to_sync;
				if (ready_to_sync == solvers->size()) {
					return true;
				}
				return false;
			}

			~sph_model() {
				delete this->serializer;
			}

			const int get_cell_num_x() {
				return cell_num_x;
			}

			const int get_cell_num_y() {
				return cell_num_y;
			}

			const int get_cell_num_z() {
				return cell_num_z;
			}

			const int get_total_cell_num() {
				return total_cell_num;
			}

            std::shared_ptr<sibernetic::solver::i_sort_solver> get_sort_solver() {
                return sort_solver;
            }

            void set_sort_solver(std::shared_ptr<sibernetic::solver::i_sort_solver> new_solver) {
                sort_solver = new_solver;
            }

			const std::vector<partition>& get_partition() {
				return partitions;
			}

			void set_solver_container(std::vector<std::shared_ptr<sibernetic::solver::i_solver>> *s){
				solvers = s;
			}
		private:
		    std::vector<float> balance_coeff;
			abstract_reader<T> *serializer;
			size_t next_partition;
			std::shared_ptr<sibernetic::solver::i_sort_solver> sort_solver;
			// vars block end
			int cell_num_x;
			int cell_num_y;
			int cell_num_z;
			int total_cell_num;

			timespec t0, t1, t2;
			timespec t3, t4;
			double us;
			double elapsedTime;

			std::vector<std::shared_ptr<sibernetic::solver::i_solver>> *solvers;
			std::atomic<int> iteration;
			std::atomic<int> ready_to_sync;

			std::mutex sync_mutex;
			container particles;
			sph_config config;
			std::map<std::string, T> phys_consts;
			std::vector<partition> partitions;

			void push_partition(size_t start, size_t end) {
				auto p = prepare_partition(start, end);
				partitions.push_back(p);
			}

            void update_partition_distrib() {
                size_t dev_count = balance_coeff.size();
                next_partition = 0;
                if (dev_count == 1) {
                    update_partition(0, size() - 1, 0);
                    return;
                }
                auto part_size = static_cast<size_t>(size() / dev_count);
                size_t start = 0;
                size_t end = 1;
                for (size_t i = 0; i < dev_count; ++i) {
                    part_size = static_cast<size_t>((float)particles.size() * balance_coeff[i]);
                    if (i == dev_count - 1)
                        update_partition(start, static_cast<size_t>(size() - 1), i);
                    else {
                        if (particles[start + part_size].cell_id != particles[start + part_size + 1].cell_id) {
                            update_partition(start, start + part_size, i);
                            start += part_size + 1;
                        } else {
                            end = start + part_size;
                            while (particles[end].cell_id == particles[end + 1].cell_id) {
                                if (end + 1 > size()) {
                                    break;
                                }
                                ++end;
                            }
                            if ((particles[end].cell_id + 1) % (cell_num_y * cell_num_z) != 0) {
                                auto rest = (particles[end].cell_id + 1) % (cell_num_y * cell_num_z);
                                auto last_cell = particles[end].cell_id + (cell_num_y * cell_num_z - rest) + 1;
                                while (particles[end].cell_id != last_cell) {
                                    ++end;
                                }
                                update_partition(start, end - 1, i);
                                start = end;
                            } else {
                                update_partition(start, end, i);
                                start = end + 1;
                            }
                        }
                    }
                }
            }

            void update_partition(size_t start, size_t end, size_t p_id) {
                auto p = prepare_partition(start, end);
                partitions[p_id].start = p.start;
                partitions[p_id].end = p.end;
                partitions[p_id].ghost_start = p.ghost_start;
                partitions[p_id].ghost_end = p.ghost_end;
                partitions[p_id].start_cell_id = p.start_cell_id;
                partitions[p_id].end_cell_id = p.end_cell_id;
                partitions[p_id].start_ghost_cell_id = p.start_ghost_cell_id;
                partitions[p_id].end_ghost_cell_id = p.end_ghost_cell_id;
            }

            partition prepare_partition(size_t start, size_t end) {
                auto start_cell_id = particles[start].cell_id;
                auto end_cell_id = particles[end].cell_id;
                auto ghost_end = end;
                auto ghost_start = start;
                size_t start_ghost_cell_id = 0, end_ghost_cell_id = end_cell_id;
                if (start_cell_id != 0) {
                    start_ghost_cell_id = start_cell_id - cell_num_y * cell_num_z;
                    if (start_ghost_cell_id < 1) {
                        ghost_start = 0;
                    } else {
                        while ( particles[ghost_start].cell_id > start_ghost_cell_id - 1) {
                            if(ghost_start > 0) {
                                --ghost_start;
                            } else {
                                break;
                            }
                        }
                        if(ghost_start != 0) {
                            ghost_start++;
                        }
                    }
                }

                if (end_cell_id != total_cell_num - 1) {
                    end_ghost_cell_id = end_cell_id + cell_num_y * cell_num_z;
                    if (end_ghost_cell_id + 1 >= total_cell_num) {
                        ghost_end = particles.size() - 1;
                    } else {
                        while (particles[ghost_end].cell_id < end_ghost_cell_id + 1) {
                            if(ghost_end < size() - 1) {
                                ++ghost_end;
                            } else {
                                break;
                            }
                        }
                        if(ghost_end != size() - 1) {
                            --ghost_end;
                        }
                    }
                }

                return partition{
                    start,
                    end,
                    ghost_start,
                    ghost_end,
                    start_cell_id,
                    end_cell_id,
                    start_ghost_cell_id,
                    end_ghost_cell_id
                };
			}

			/** Init variables for simulation
			 */
			void init_vars() {
				cell_num_x =
						static_cast<int>((config["x_max"] - config["x_min"]) * GRID_CELL_SIZE_INV);
				cell_num_y =
						static_cast<int>((config["y_max"] - config["y_min"]) * GRID_CELL_SIZE_INV);
				cell_num_z =
						static_cast<int>((config["z_max"] - config["z_min"]) * GRID_CELL_SIZE_INV);
				total_cell_num = cell_num_x * cell_num_y * cell_num_z;

				config["h_scaled"] = H * config["simulation_scale"];
				config["simulation_scale_inv"] = 1.0f / config["simulation_scale"];
				config["h_scaled_2"] = config["h_scaled"] * config["h_scaled"];
			}

			/**Arrange particles according its cell id
			 * it will need for future clustering
			 * particles array on several devices.
			 * TODO make sort parallel
			 */
			void arrange_particles() {
			    if(this->sort_solver == nullptr) {
                    std::sort(particles.begin(), particles.end(),
                              [](const particle<T> &p1, const particle<T> &p2) {
                                  return p1.cell_id < p2.cell_id;
                              });
                } else {
			        this->sort_solver->sort();
			    }
			}

			// Addition methods
			/** TODO Description here
			 */
			void calc_grid_id(particle<T> &p) {
				int A, B, C;
				A = static_cast<int>(p.pos[0] * GRID_CELL_SIZE_INV);
				B = static_cast<int>(p.pos[1] * GRID_CELL_SIZE_INV);
				C = static_cast<int>(p.pos[2] * GRID_CELL_SIZE_INV);
				//p.cell_id = A + B * cell_num_x + cell_num_x * cell_num_y * C; // this stats indexing from x component
				//p.cell_id = (B + C * cell_num_y + cell_num_y * cell_num_z * A); // now will indexing from y
                p.cell_id = ( std::rand() % ( 10 + 1 ) );
			}
		};
	} // namespace model
} // namespace sibernetic
#endif // SPHMODEL_HPP
