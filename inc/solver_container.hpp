#ifndef SOLVER_CONTAINER_HPP
#define SOLVER_CONTAINER_HPP

#include "isolver.h"
#include "ocl_const.h"
#include "ocl_solver.hpp"
#include "isort_solver.h"
#include "sph_model.hpp"
#include "util/ocl_helper.h"
#include "util/error.h"
#include <string>
#include <vector>
#include <thread>
#include <cmath>

namespace sibernetic {
	namespace solver {
		using model::sph_model;
		using std::shared_ptr;
		using sibernetic::solver::ocl_solver;
        using sibernetic::solver::i_sort_solver;
		enum EXECUTION_MODE {
			ONE,
			ALL
		};

		template<class T = float>
		class solver_container {
			typedef shared_ptr<sph_model<T>> model_ptr;

		public:
			solver_container(const solver_container &) = delete;

			solver_container &operator=(const solver_container &) = delete;

			/** Maer's singleton
			 */
			static solver_container &instance(model_ptr &model,
					EXECUTION_MODE mode = EXECUTION_MODE::ONE, size_t dev_count = 0,
			        SOLVER_TYPE s_t = OCL) {
				static solver_container s(model, mode, s_t, dev_count);
				model->set_solver_container(s.solvers());
				return s;
			}

			void run() {
				std::vector<std::thread> t_pool;
				std::for_each(
						_solvers.begin(), _solvers.end(), [&, this](std::shared_ptr<i_solver> &s) {
							t_pool.emplace_back(std::thread(solver_container::run_solver, std::ref(s)));
						});
				std::for_each(t_pool.begin(), t_pool.end(), [](std::thread &t) { t.join(); });
			}

			static void run_solver(std::shared_ptr<i_solver> &s) {
				s->run();
			}
			std::vector<std::shared_ptr<i_solver>>* solvers(){
				return &_solvers;
			}
		private:
			explicit solver_container(model_ptr &model, EXECUTION_MODE mode = EXECUTION_MODE::ONE,
			                          SOLVER_TYPE s_type = OCL, size_t dev_count = 0) {
				try {
					p_q dev_q = get_dev_queue();
					size_t device_index = 0;
					while (!dev_q.empty()) {
						try {
							std::shared_ptr<ocl_solver<T>> solver(
									new ocl_solver<T>(model, dev_q.top(), device_index));
							_solvers.push_back(solver);
							std::cout << "************* DEVICE *************" << std::endl;
							dev_q.top()->show_info();
							std::cout << "**********************************" << std::endl;
							++device_index;
							if(mode == EXECUTION_MODE::ONE){
								break;
							} else if(mode == EXECUTION_MODE::ALL){
								if(dev_count > 0 && dev_count == device_index){
									break;
								}
							}
						} catch (ocl_error &ex) {
							std::cout << ex.what() << std::endl;
						}
						dev_q.pop();
					}

					if (_solvers.size()) {
                        init_weights();
                        model->set_balance_vector(this->weights);
						model->make_partition(_solvers.size()); // TODO to think about is in future we
                        //model->make_partition(_solvers.size(), std::vector<float>{0.5, 0.4, 0.1});
						// can't init one or more
						// devices
						// obvious we should reinit partitions case ...
						int i=0;
						for (auto s : _solvers) {
						    s->get_device()->balance_coeff = weights[i++];
							s->init_model(&(model->get_next_partition()));
						}
					} else {
						throw ocl_error("No OpenCL devices were initialized.");
					}
				} catch (sibernetic::ocl_error &err) {
					throw;
				}
			}

            void init_weights(){
			    size_t total_compute_power = 0;
			    std::for_each(_solvers.begin(), _solvers.end(),
			            [&total_compute_power](std::shared_ptr<i_solver> s) -> void {
			                total_compute_power += s->get_device()->max_thread_count;
			        }
			    );
                float prev = 0.f;
                float rest = 1.f;
			    for(int i=0; i < _solvers.size(); ++i){
                    float cur;
			        if(i == _solvers.size() - 1) {
			            cur = ceilf(rest * 100) / 100;
			        }else {
			            auto s = _solvers[i];
			            cur = static_cast<float>(s->get_device()->max_thread_count) / static_cast<float>(total_compute_power);
			            cur = (cur * 100 + 0.5f) / 100.f;
			            rest -= cur;
			        }
			        weights.push_back(cur);
			    }
			}

			~solver_container() = default;
            std::vector<float> weights;
			std::vector<std::shared_ptr<i_solver>> _solvers;
            std::shared_ptr<i_sort_solver> sort_solver;
		};
	} // namespace solver
} // namespace sibernetic

#endif // SOLVER_CONTAINER_HPP