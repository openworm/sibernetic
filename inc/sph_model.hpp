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
#include "util/custom_reader.hpp"
#include "abstract_model.hpp"
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace sibernetic {
namespace model {
const float H = 3.34f;
const float H_INV = 1.f / H;
const float GRID_CELL_SIZE = 2.0f * H;
const float GRID_CELL_SIZE_INV = 1 / GRID_CELL_SIZE;
const float R_0 = 0.5f * H;
/* const block end */
template <class T = float, class container = std::vector<particle<T>>>
class sph_model: public particle_model<T> {
  typedef std::map<std::string, T> sph_config;

public:
  sph_model(const std::string &config_file, abstract_reader<T> * serializer = new custom_reader<T>()):serializer(serializer) {
    config = {{"particles", T()}, {"x_max", T()}, {"x_min", T()},
              {"y_max", T()},     {"y_min", T()}, {"z_max", T()},
              {"z_min", T()},     {"mass", T()},  {"time_step", T()},
              {"rho0", T()}};
    this->serializer->serialize(config_file, this);
    arrange_particles();
    init_vars();
    std::cout << "Model was loaded: " << particles.size() << " particles."
              << std::endl;
  }
  const sph_config &get_config() const { return config; }
  const container &get_particles() const { return particles; }
  container &get_particles() { return particles; }
  const particle<T> & get_particle(const int index) const {
    return this->particles.at(index);
  }
  particle<T> & get_particle(const int index){
    return this->particles.at(index);
  }
  void set_particle(int index,const particle<T> & p) {
    if(p.cell_id == 0) {
      this->particles.at(index) = p;
      calc_grid_id(this->particles.at(index));
    } else {
      this->particles.at(index) = p;
    }
  }
  void push_back(const particle<T>& p){
    if(p.cell_id == 0) {
      this->particles.push_back(p);
      calc_grid_id(this->particles.back());
    } else {
      this->particles.push_back(p);
    }
  }
  size_t size() const { return particles.size(); }
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
  sph_config & get_config(){
    return this->config;
  }
  ~sph_model(){
    delete this->serializer;
  }
private:
  abstract_reader<T> * serializer;
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
  /**Arrange particles according its cell id
   * it will need for future clustering
   * particles array on several devices.
   * TODO make sort parralel
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
