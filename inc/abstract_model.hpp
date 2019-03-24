//
// Created by sergey on 12.02.19.
//

#ifndef SIBERNETIC_ABSTRACT_MODEL_HPP
#define SIBERNETIC_ABSTRACT_MODEL_HPP

#include "particle.h"
#include <map>

namespace sibernetic {
    namespace model {
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
