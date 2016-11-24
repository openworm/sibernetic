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
#include "owVtkExport.h"

bool load_from_file = false;
bool load_to = false;

int usage(){
	std::string version = "0.0.5b";
	std::cout << "\nSibernetic v" << version << "\n  This is a C++ "
	          << "implementation of the Contractile SPH (Electrofluid) "
			  << "algorithm applied to C. elegans locomotion\n\n"
	 		  << "Usage: ./Release/Sibernetic [OPTION]\n\n"
			  << "    -no_g                      Run without graphics\n\n"
			  << "    -l_to                      Save simulation results to disk\n\n"
			  << "    -export_vtk                Save simulation results to VTK files\n\n"
			  << "    logstep=<value>            Set frequency of logging data into file in -l_to\n"
			  << "                               and -export_vtk modes by default it equals to 10\n\n"
			  << "    -l_from                    Load simulation results from disk\n\n"
			  << "    lpath=<value>              Indicates path where all buffers will be stored \n"
			  << "                               this option also works for -l_to and -l_from options\n\n"
			  << "    -test                      Run some physical tests\n\n"
			  << "    -f <filename>              Load configuration from file "
			  << "./configuration/<filename>\n\n"
			  << "    -f worm                    **Load Worm Body Simulation**\n\n"
			  << "    device=<device_type>       Trying to init OpenCL on device <type>\n"
			  << "                               it could be cpu or gpu default-ALL (it try to init "
			  << "most powerful available device)\n\n"
			  << "    timestep=<value>           Start simulation with time "
			  << "step = <value> in seconds\n\n"
			  << "    timelimit=<value>          Run simulation until <value> will be "
			  << "reached in seconds\n\n"
			  << "    leapfrog                   Run simulation using Leapfrog integration "
			  << "method for time integration\n\n"
			  << "    oclsourcepath=<value>      You can indicate path to you'r "
			  << "OpenCL program just using this option\n\n"
			  << "    -nrn <value>                Indicates that you plan run simulation with "
			  << "NEURON simulation = <value> \n"
			  << "                               value should be a file which "
			  << "can be run by NEURON simulator and \n"
			  << "                               also you should have installed neuron\n"
			  << "                               and sibernetic_neuron bridge\n\n"
			  << "    -c302                      Run worm model with c302 (use sibernetic_c302.py)\n\n"
			  << "    -help, -h, -?, --help      Print this information\n\n"
			  << "Full documentation at: <https://github.com/openworm/sibernetic>\n"
			  << "Please report any bugs/issues "
			  << "on: <https://github.com/openworm/sibernetic/issues>\n";
	return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
	int exitStatus;
    if (argc == 1) {
        std::cout << "Sibernetic: no arguments specified, run method executing\n";
        exitStatus = run(argc, argv);
    } else {
        bool graph = true;
        bool run_tests = false;

        for (int i = 1; i < argc; i++) {
            if (std::string("-help").compare(argv[i]) == 0 ||
				std::string("-?").compare(argv[i]) == 0 ||
				std::string("--help").compare(argv[i]) == 0 ||
				std::string("-h").compare(argv[i]) == 0) { // print usage information
				return usage();
            }
			if (std::string("-no_g").compare(argv[i]) == 0) // run without graphics
                graph = false;
            if (std::string("-l_to").compare(argv[i]) == 0) { // run load config to file mode
                std::cout << "-l_to flag: Sibernetic will save simulation results to disk\n";
                load_to = true;
            }
			if (std::string("-export_vtk").compare(argv[i]) == 0) {
				owVtkExport::isActive = true;
			}
            if (std::string("-l_from").compare(argv[i]) == 0) { // run load config from file mode
                graph = true;
                load_from_file = true;
            }
            if (std::string("-test").compare(argv[i]) == 0) { // run tests
                run_tests = true;
            }
        }
        if (run_tests) {
            test_energy_conservation(argc, argv);
        } else
        	exitStatus = run(argc, argv, graph);
    }
    return exitStatus;
}
