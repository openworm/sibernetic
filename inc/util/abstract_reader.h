//
// Created by sergey on 12.02.19.
//

#ifndef SIBERNETIC_ABSTRARACT_READER_H
#define SIBERNETIC_ABSTRARACT_READER_H

#include "../abstract_model.hpp"

template<class T> class AbstractReader{
public:
    void serialize(const std::string &, IParticleModel<T> * ) = 0;
};

#endif //SIBERNETIC_ABSTARACT_READER_H
