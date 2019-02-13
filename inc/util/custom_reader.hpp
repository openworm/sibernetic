//
// Created by sergey on 12.02.19.
//

#ifndef SIBERNETIC_JSON_READER_HPP
#define SIBERNETIC_JSON_READER_HPP

#include "abstract_reader.h"
#include "error.h"

namespace sibernetic {
    namespace model {

        template<class T>
        class custom_reader : abstract_reader<T> {
            enum LOADMODE {
                NOMODE = -1, PARAMS, MODEL, POS, VEL
            };
        public:
            serialize(const string &file_name, IParticleModel <T> *model) {
                std::ifstream file(file_name.c_str(), std::ios_base::binary);
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
                        if (cur_line.compare("parameters[") == 0) {
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
                                    model->push_back(p)
                                    break;
                                }
                                case VEL: {
                                    if (index >= particles.size())
                                        throw parser_error(
                                                "Config file problem. Velocities more than partiles is.");
                                    model->get_particle(index).vel = *get_vector(cur_line);
                                    ++index;
                                    break;
                                }
                                default: {
                                    break;
                                }
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

        private:
            std::shared_ptr <std::array<T, 4>> get_vector(const std::string &line) {
                std::shared_ptr <std::array<T, 4>> v(new std::array<T, 4>());
                std::stringstream ss(line);
                ss >> (*v)[0] >> (*v)[1] >> (*v)[2] >> (*v)[3]; // TODO check here!!!
                return v;
            }
        };
    }
}
#endif //SIBERNETIC_JSON_READER_HPP
