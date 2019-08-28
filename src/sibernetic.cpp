/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
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
#include "solver_container.hpp"
#include "util/arg_parser.h"
#include "graph.h"
#include <iostream>
#include <thread>

using sibernetic::model::sph_model;
using sibernetic::solver::solver_container;
using sibernetic::solver::EXECUTION_MODE;
using sibernetic::graphics::graph;
using sibernetic::graphics::UI_MODE;
using sibernetic::graphics::g_config;


int main(int argc, char **argv) {
  arg_parser prsr(argc, argv);
  if (prsr.check_arg("-h") || prsr.check_arg("--help") ||
      prsr.check_arg("-?") || prsr.check_arg("-help")) {
    return arg_parser::show_usage();
  }
  std::string model_name;
  auto mode = EXECUTION_MODE::ONE;
  auto ui_mode = UI_MODE::CLI;
  size_t dev_count = 0;
  if (prsr.check_arg("-f")) {
    model_name = prsr.get_arg("-f");
  } else {
    model_name = "config/tmp";
  }
  if (prsr.check_arg("--multi_dev")) {
  	auto v = prsr.get_arg_value("--multi_dev");
  	if(v != "ALL"){
		dev_count = std::stol(v);
  	}
  	mode = EXECUTION_MODE::ALL;
  }
  if (prsr.check_arg("--with_g")) {
    ui_mode = UI_MODE::OGL;
  }
  try {
    std::shared_ptr<sph_model<float>> model(new sph_model<float>(model_name));
    auto config = new g_config{
      model->get_config()["x_min"],
      model->get_config()["y_min"],
      model->get_config()["z_min"],
      model->get_config()["x_max"],
      model->get_config()["y_max"],
      model->get_config()["z_max"],
    };
    graph::config = config;
    graph::model = model;
    solver_container<float> &s_con = solver_container<float>::instance(model, mode, dev_count);
    model->get_sort_solver()->sort();
//    if(ui_mode == UI_MODE::OGL) {
//      std::thread graph_thread(graph::run, argc, argv);
//      s_con.run();
//      graph_thread.join();
//    } else {
//      s_con.run();
//    }
  } catch (sibernetic::parser_error &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (sibernetic::ocl_error &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
