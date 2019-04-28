#ifndef SOLVER_CONTAINER_HPP
#define SOLVER_CONTAINER_HPP

#include "isolver.h"
#include "ocl_const.h"
#include "ocl_solver.hpp"
#include "sph_model.hpp"
#include "util/ocl_helper.h"
#include "util/error.h"
#include <string>
#include <vector>
#include <thread>

namespace sibernetic {
	namespace solver {
		using model::sph_model;
		using std::shared_ptr;
		using sibernetic::solver::ocl_solver;

		template<class T = float>
		class solver_container {
			typedef shared_ptr<sph_model<T>> model_ptr;

		public:
			solver_container(const solver_container &) = delete;

			solver_container &operator=(const solver_container &) = delete;

			/** Maer's singleton
			 */
			static solver_container &instance(model_ptr &model, size_t devices_number = 1,
			                                  SOLVER_TYPE s_t = OCL) {
				static solver_container s(model, devices_number, s_t);
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
			explicit solver_container(model_ptr &model, size_t devices_number = 1,
			                          SOLVER_TYPE s_type = OCL) {
				try {
					std::priority_queue<std::shared_ptr<device>> dev_q = get_dev_queue();
					size_t device_index = 0;
					while (!dev_q.empty()) {
						try {
							std::shared_ptr<ocl_solver<T>> solver(
									new ocl_solver<T>(model, dev_q.top(), device_index));
							_solvers.push_back(solver);
							++device_index;
						} catch (ocl_error &ex) {
							std::cout << ex.what() << std::endl;
						}
						dev_q.pop();
					}
					if (_solvers.size()) {
						model->make_partition(_solvers.size()); // TODO to think about is in future we
						// can't init one or more
						// devices
						// obvious we should reinit partitions case ...
						for (auto s : _solvers) {
							s->init_model(model->get_next_partition());
						}
					} else {
						throw ocl_error("No OpenCL devices were initialized.");
					}
				} catch (sibernetic::ocl_error &err) {
					throw;
				}
			}

			~solver_container() = default;

			std::vector<std::shared_ptr<i_solver>> _solvers;
		};
	} // namespace solver
} // namespace sibernetic

#endif // SOLVER_CONTAINER_HPP