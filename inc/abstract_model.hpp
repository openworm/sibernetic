//
// Created by sergey on 12.02.19.
//

#ifndef SIBERNETIC_ABSTRACT_MODEL_HPP
#define SIBERNETIC_ABSTRACT_MODEL_HPP

#include "particle.h"
#include <map>

namespace sibernetic {
    namespace model {
	    const float H = 3.34f;
	    const float H_INV = 1.f / H;
	    const float GRID_CELL_SIZE = 2.0f * H;
	    const float GRID_CELL_SIZE_INV = 1 / GRID_CELL_SIZE;
	    const float R_0 = 0.5f * H;
	    const float DEFAULT_MASS = 20.00e-13f;
		const float DEFAULT_DENSITY = 1000.f;
		const float DEFAULT_MU = 0.00005f;
		const float DENSITY_WATER = 1000.0f;
		const int PCI_ITER_COUNT = 3;
	    template<class T> class particle_model {
        public:
          virtual std::map<std::string, T> & get_config() = 0;

          virtual const std::map<std::string, T> & get_config() const = 0;

          virtual size_t size() const = 0;

          virtual const particle<T> &get_particle(int) const = 0;

          virtual  particle<T> &get_particle(int) = 0;

          virtual void set_particle(int index, const particle<T> &) = 0;

          virtual void push_back(const particle<T> &) = 0;

          virtual void sync() = 0;

          virtual bool set_ready() = 0;
    };
}
}
#endif //SIBERNETIC_ABSTRACT_MODEL_HPP
