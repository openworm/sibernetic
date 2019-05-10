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
#ifndef PARTICLE_H
#define PARTICLE_H

#include <algorithm>
#include <array>
#include <sstream>

namespace sibernetic {
namespace model {
// TODO write the docs
// Write why alligment on 16 bytes is important!!
// Also in this structure we store all enought information for particle which we
// want to load from device it should be as small as possible for optimal
template <class T, size_t dim = 4> struct alignas(16) particle {
  typedef std::array<T, dim> container;
  container pos;
  container pos_n_1;
  container vel;
  container accelation;
  container accelation_n_1;
  container accelation_n_0_5;
  int type;
  int cell_id;
  int particle_id;
  T density;
  T pressure;
  T viscosity;
  T mass;
  size_t get_dim() const { return dim; }
  std::string pos_str() {
    std::stringstream s;
    std::for_each(pos.begin(), pos.end(), [&s](T c) { s << c << ' '; });
    s << '\n';
    return s.str();
  }
  bool operator<(const particle<T> &p) { return cell_id < p.cell_id; }
};
} // namespace model
} // namespace sibernetic

#endif