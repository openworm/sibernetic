//
// Created by sergey on 12.02.19.
//

#ifndef SIBERNETIC_ABSTRARACT_READER_H
#define SIBERNETIC_ABSTRARACT_READER_H

#include "../abstract_model.hpp"

template<class T> class abstract_reader{
public:
    virtual void serialize(const std::string &, sibernetic::model::particle_model<T> * ) = 0;
};

#endif //SIBERNETIC_ABSTARACT_READER_H
