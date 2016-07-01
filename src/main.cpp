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


int main(int argc, char **argv) {
	int exitStatus;
    if (argc == 1) {
        std::cout << "Sibernetic: no arguments specified, run method executing\n";
        exitStatus = run(argc, argv);
    } else {
        bool graph = true;
        bool run_tests = false;
        std::string version = "0.0.4b";
        std::string helpFlag = "-help";
        std::string helpFlag2 = "-h";
        std::string helpFlag3 = "-?";
        std::string helpFlag4 = "--help";
        std::string noGraphicsFlag = "-no_g";
        std::string saveFlag = "-l_to";
        std::string loadFlag = "-l_from";
		std::string exportVtkFlag = "-export_vtk";
        std::string testFlag = "-test";
        std::string configFileFlag = "-f <filename>";
        std::string configFileWorm = "-f worm";
        std::string deviceTypeFlag = "device=<device_type>";
        std::string timeStepFlag = "timestep=<value>";
        std::string timeLimitFlag = "timelimit=<value>";
        std::string integrationMethodFlag = "leapfrog";
        std::string logStepFlag = "logstep=<value>";
        std::string loadFlagPath = "lpath=<value>";
        std::string oclSourcePath = "oclsourcepath=<value>";
        std::string nrnFlag = "nrn <value>";
        for (int i = 1; i < argc; i++) {
            if (helpFlag.compare(argv[i]) == 0 ||
            	helpFlag2.compare(argv[i]) == 0||
				helpFlag3.compare(argv[i]) == 0||
				helpFlag4.compare(argv[i]) == 0) { // print usage information
                std::cout << "\nSibernetic v" << version << "\n  This is a C++ implementation of the Contractile SPH (Electrofluid) algorithm applied to C. elegans locomotion\n\n";
                std::cout << "  Usage:  ./Release/Sibernetic [OPTION]\n\n";
                std::cout << "    " << noGraphicsFlag << "                      Run without graphics\n\n";
                std::cout << "    " << saveFlag << "                      Save simulation results to disk\n\n";
				std::cout << "    " << exportVtkFlag << "                Save simulation results to VTK files\n\n";
                std::cout << "        " << logStepFlag << "           Set frequency of logging data into file in -l_to and -export_vtk modes by default it equals to 10\n\n";
                std::cout << "    " << loadFlag << "                    Load simulation results from disk\n\n";
                std::cout << "        " << loadFlagPath << "                    Indicates path where all buffers will be stored this option also works for -l_to and -l_from options\n\n";
                std::cout << "    " << testFlag << "                      Run some tests\n\n";
                std::cout << "    " << configFileFlag << "              Load configuration from file ./configuration/<filename>\n\n";
                std::cout << "        " << configFileWorm << "                    **Load Worm Body Simulation**\n\n";
                std::cout << "    " << deviceTypeFlag << "       Trying to init OpenCL on device <type> it could be cpu or gpu default-ALL (it try to init most powerful available device)\n\n";
                std::cout << "    " << timeStepFlag << "           Start simulation with time step = <value> in seconds\n\n";
                std::cout << "    " << timeLimitFlag << "          Run simulation until <value> will be reached in seconds\n\n";
                std::cout << "    " << integrationMethodFlag << "                   Run simulation using Leapfrog integration method for time integration\n\n";
                std::cout << "    " << oclSourcePath << "      You can indicate path to you'r OpenCL program just using this option\n\n";
                std::cout << "    " << nrnFlag << "      Indicates that you plan run simulation with NEURON simulation = <value> value should be a file which can be run by NEURON simulator and also you should have installed neuron and sibernetic_neuron bridg\n\n";
                std::cout << "    " << helpFlag << "                      Print this information\n\n";
                std::cout << "  Please report any bugs/issues on: https://github.com/openworm/sibernetic/issues\n\n";
                return 0;
            }
            if (noGraphicsFlag.compare(argv[i]) == 0) // run without graphics
                graph = false;
            if (saveFlag.compare(argv[i]) == 0) { // run load config to file mode
                std::cout << saveFlag << " flag: Sibernetic will save simulation results to disk\n";
                load_to = true;
            }
			if (exportVtkFlag.compare(argv[i]) == 0) {
				owVtkExport::isActive = true;
			}
            if (loadFlag.compare(argv[i]) == 0) { // run load config from file mode
                graph = true;
                load_from_file = true;
            }
            if (testFlag.compare(argv[i]) == 0) { // run tests
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
