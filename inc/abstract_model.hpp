//
// Created by sergey on 12.02.19.
//

#ifndef SIBERNETIC_ABSTRACT_MODEL_HPP
#define SIBERNETIC_ABSTRACT_MODEL_HPP

#include "particle.h"

namespace sibernetic {
    namespace model {
        template<class T> IParticleModel {
        public:
          const particle &get_particle(const int) const = 0;

          particle &get_particle(const int) = 0;

          void set_particle(const int index, const partile &) = 0;

          void push_back(const particle &) = 0;
    };
}
}
#endif //SIBERNETIC_ABSTRACT_MODEL_HPP
