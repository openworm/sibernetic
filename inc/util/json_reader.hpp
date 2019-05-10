//
// Created by serg on 31.03.19.
//

#ifndef SIBERNETIC_JSON_READER_HPP
#define SIBERNETIC_JSON_READER_HPP

#include <fstream>
#include "json.hpp"
#include "abstract_reader.h"

using json = nlohmann::json;
using out_of_range = nlohmann::detail::out_of_range;

template<class T> class json_reader: public abstract_reader<T>{
public:
	void serialize(const std::string & file_name, sibernetic::model::particle_model<T> *  model) override{
		std::ifstream file(file_name.c_str(), std::ios_base::binary);
		json j;
		file >> j;
		auto o = j.at("parameters");
		for (json::iterator it = o.begin(); it != o.end(); ++it) {
			model->get_config()[it.key()] =  it.value();
		}
		int particle_id = 0;
		// range-based for
		for (auto& element : j.at("model")) {
			sibernetic::model::particle<T> p;
			p.pos = element.at("position");
			p.vel = element.at("velocity");
			p.type = element.at("type");
			p.particle_id = particle_id++;
			try {
				p.mass = element.at("mass");
			}
			catch(out_of_range & e){
				p.mass = sibernetic::model::DEFAULT_MASS;
			}
			try {
				p.density = element.at("density");
			}
			catch(out_of_range & e){
				p.density = sibernetic::model::DEFAULT_DENSITY;
			}

			model->push_back(p);
		}
	}
};


#endif //SIBERNETIC_JSON_READER_HPP
