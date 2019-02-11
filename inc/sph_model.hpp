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
#include "util/error.h"
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <string>
namespace sibernetic {
namespace model {
enum LOADMODE { NOMODE = -1, PARAMS, MODEL, POS, VEL };
const float H = 3.34f;
const float H_INV = 1.f / H;
const float GRID_CELL_SIZE = 2.0f * H;
const float GRID_CELL_SIZE_INV = 1 / GRID_CELL_SIZE;
const float R_0 = 0.5f * H;
/* const block end */
struct partition {
  /**each device has its own partition
   * in which we define where starts
   * and end particles for this device.
   */
  size_t start;
  size_t end;
  size_t size() const { return end - start; }
};
template <class T = float, class container = std::vector<particle<T>>>
class sph_model {
  typedef std::map<std::string, T> sph_config;

public:
  sph_model(const std::string &config_file) {
    config = {{"particles", T()}, {"x_max", T()}, {"x_min", T()},
              {"y_max", T()},     {"y_min", T()}, {"z_max", T()},
              {"z_min", T()},     {"mass", T()},  {"time_step", T()},
              {"rho0", T()}};
    read_model(config_file);
    arrange_particles();
    std::cout << "Model was loaded: " << particles.size() << " particles."
              << std::endl;
  }
  const sph_config &get_config() const { return config; }
  const container &get_particles() const { return particles; }
  container &set_particles() { return particles; }
  int size() const { return particles.size(); }
  /** Make partition for device
   */
  void make_partition(size_t dev_count) {
    next_partition = 0;
    if (dev_count == 1) {
      partitions.push_back(partition{0, static_cast<size_t>(size() - 1)});
      return;
    }
    size_t part_size = static_cast<size_t>(size() / dev_count);
    size_t start = 0 * part_size;
    size_t end = 1 * part_size;
    for (size_t i = 0; i < dev_count; ++i) {
      if (i == dev_count - 1)
        partitions.push_back(partition{start, static_cast<size_t>(size() - 1)});
      else {
        if (particles[end].cell_id != particles[end + 1].cell_id) {
          partitions.push_back(partition{start, end});
          start = end;
        } else {
          for (; end < particles.size() - 1; ++end) {
            if (particles[end].cell_id != particles[end + 1].cell_id) {
              ++end;
              break;
            }
          }
          partitions.push_back(partition{start, end});
          start = end;
        }
      }
    }
  }
  const partition &get_next_partition() {
    ++next_partition;
    return partitions[next_partition - 1];
  }

private:
  size_t next_partition;
  // vars block end
  int cell_num_x;
  int cell_num_y;
  int cell_num_z;
  long total_cell_num;

  container particles;
  sph_config config;
  std::map<std::string, T> phys_consts;
  std::vector<partition> partitions;
  /** Init variables for simulation
   */
  void init_vars() {
    cell_num_x =
        static_cast<int>((config["x_max"] - config["x_min"]) / GRID_CELL_SIZE);
    cell_num_y =
        static_cast<int>((config["y_max"] - config["y_min"]) / GRID_CELL_SIZE);
    cell_num_z =
        static_cast<int>((config["z_max"] - config["z_min"]) / GRID_CELL_SIZE);
    total_cell_num = cell_num_x * cell_num_y * cell_num_z;
  }
  std::shared_ptr<std::array<T, 4>> get_vector(const std::string &line) {
    std::shared_ptr<std::array<T, 4>> v(new std::array<T, 4>());
    std::stringstream ss(line);
    ss >> (*v)[0] >> (*v)[1] >> (*v)[2] >> (*v)[3]; // TODO check here!!!
    return v;
  }
  /**Model reader
   * Read the model from file and load into memory
   */
  void read_model(const std::string &model_file) {
    std::ifstream file(model_file.c_str(), std::ios_base::binary);
    LOADMODE mode = NOMODE;
    bool is_model_mode = false;
    int index = 0;
    if (file.is_open()) {
      while (file.good()) {
        std::string cur_line;
        std::getline(file, cur_line);
        cur_line.erase(std::remove(cur_line.begin(), cur_line.end(), '\r'),
                       cur_line.end()); // crlf win fix
        auto i_space = cur_line.find_first_not_of(" ");
        auto i_tab = cur_line.find_first_not_of("\t");
        if (i_space) {
          cur_line.erase(cur_line.begin(), cur_line.begin() + i_space);
        }
        if (i_tab) {
          cur_line.erase(cur_line.begin(), cur_line.begin() + i_tab);
        }
        if (cur_line.compare("parametrs[") == 0) {
          mode = PARAMS;
          continue;
        } else if (cur_line.compare("model[") == 0) {
          mode = MODEL;
          is_model_mode = true;
          init_vars();
          continue;
        } else if (cur_line.compare("position[") == 0) {
          mode = POS;
          continue;
        } else if (cur_line.compare("velocity[") == 0) {
          mode = VEL;
          continue;
        } else if (cur_line.compare("]") == 0) {
          mode = NOMODE;
          continue;
        }
        if (mode == PARAMS) {
          std::regex rgx("^\\s*(\\w+)\\s*:\\s*(\\d+(\\.\\d*([eE]?[+-]?\\d+)?)?)"
                         "\\s*(//.*)?$");
          std::smatch matches;
          if (std::regex_search(cur_line, matches, rgx)) {
            if (matches.size() > 2) {
              if (config.find(matches[1]) != config.end()) {
                config[matches[1]] = static_cast<T>(stod(matches[2].str()));
                continue;
              }
            } else {
              std::string msg = sibernetic::make_msg(
                  "Problem with parsing parametrs:", matches[0].str(),
                  "Please check parametrs.");
              throw parser_error(msg);
            }
          } else {
            throw parser_error(
                "Please check parameters section there are no parametrs.");
          }
        }
        if (is_model_mode) {
          switch (mode) {
          case POS: {
            particle<T> p;
            p.pos = *get_vector(cur_line);
            calc_grid_id(p);
            particles.push_back(p);
            break;
          }
          case VEL: {
            if (index >= particles.size())
              throw parser_error(
                  "Config file problem. Velocities more than partiles is.");
            particles[index].vel = *get_vector(cur_line);
            ++index;
            break;
          }
          default: { break; }
          }
        }
      }
    } else {
      throw parser_error(
          "Check your file name or path there is no file with name " +
          model_file);
    }
    file.close();
  }
  /**Arrange particles according its cell id
   * it will need for future clustering
   * particles array on several devices.
   */
  void arrange_particles() {
    std::sort(particles.begin(), particles.end(),
              [](const particle<T> &p1, const particle<T> &p2) {
                return p1.cell_id < p2.cell_id;
              });
  }

  // Addition methods
  /** TODO Description here
   */
  void calc_grid_id(particle<T> &p) {
    int A, B, C;
    A = static_cast<int>(p.pos[0] * GRID_CELL_SIZE_INV);
    B = static_cast<int>(p.pos[1] * GRID_CELL_SIZE_INV);
    C = static_cast<int>(p.pos[2] * GRID_CELL_SIZE_INV);
    p.cell_id = A + B * cell_num_x + cell_num_x * cell_num_y * C; //
  }
};
} // namespace model
} // namespace sibernetic
#endif // SPHMODEL_HPP
