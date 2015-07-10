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

#include <stdio.h>
#include <iostream>
#include "owWorldSimulation.h"
#include "owPhysicTest.h"

#include "owWorldSimulation.h"

bool load_from_file = false;

int main(int argc, char **argv) {

    if (argc == 1) {
        std::cout << "Sibernetic: no arguments specified, run method executing\n";
        run(argc, argv);
    } else {
        bool graph = true;
        bool load_to = false;
        bool run_tests = false;

        char help_flag[] = "-help";
        char no_graphics_flag[] = "-no_g";
        char save_flag[] = "-l_to";
        char load_flag[] = "-l_from";
        char test_flag[] = "-test";

        for (int i = 1; i < argc; i++) {
            if (strncmp(argv[i], help_flag, 6) == 0) { // print usage information
                std::cout << "\nSibernetic\n  This is a C++ implementation of the Contractile SPH (Electrofluid) algorithm applied to C. elegans locomotion\n\n";
                std::cout << "  Usage:  ./Release/Sibernetic [OPTION]\n";
                std::cout << "    " << no_graphics_flag << "       Run without graphics\n";
                std::cout << "    " << save_flag << "       Save simulation results to disk\n";
                std::cout << "    " << load_flag << "     Load simulation results from disk\n";
                std::cout << "    " << test_flag << "       Run some tests\n";
                std::cout << "    " << help_flag << "       Print this information\n\n";
                std::cout << "  Please report any bugs/issues on: https://github.com/openworm/Smoothed-Particle-Hydrodynamics/issues\n\n";
                exit(0);
            }
            if (strncmp(argv[i], no_graphics_flag, 5) == 0) // run without graphics
                graph = false;
            if (strncmp(argv[i], save_flag, 5) == 0) { // run load config to file mode
                std::cout << save_flag << " flag: Sibernetic will save simulation results to disk\n";
                graph = false;
                load_to = true;
            }
            if (strncmp(argv[i], load_flag, 7) == 0) { // run load config from file mode
                graph = true;
                load_from_file = true;
            }
            if (strncmp(argv[i], test_flag, 5) == 0) { // run tests
                run_tests = true;
            }
        }
        if (run_tests) {
            test_energy_conservation();
        } else
            run(argc, argv, graph, load_to);
    }
    return 0;
}
