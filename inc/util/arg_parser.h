/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
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
#ifndef ARG_PARSER
#define ARG_PARSER
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
// Class for parsing command line arguments
class arg_parser {
public:
  arg_parser(int argc, char **argv);
  bool check_arg(const std::string &arg) const;
  const std::string &get_arg(const std::string &arg) const;
  const std::string get_arg_value(const std::string &arg) const;
  static int show_usage() {
    std::string version = "0.0.1";
    std::cout
        << "\nsibernetic v" << version << "\n  This is a C++/OpenCL "
        << "implementation of the SPH algorithm supplemented with"
        << "with many posibilities"
        << "a set of biomechanics related features"
        << "Usage: ./bin/sibernetic [OPTION]\n\n"
        << "    --multi_dev                Run without on all available "
           "devices\n"
        << "                               but default it will run only "
           "one.\n\n"
        << "    -f <config_file>           Path to configuration file.\n\n"
        << "    -help, -h, -?, --help      Print this information\n\n"
        << "Full documentation at: <https://github.com/openworm/sibernetic>\n"
        << "Please report any bugs/issues "
        << "to: <https://github.com/openworm/sibernetic/issues>\n";
    return EXIT_SUCCESS;
  }
private:
  std::unordered_map<std::string, std::string> arguments;
};
#endif // ARG_PARSER