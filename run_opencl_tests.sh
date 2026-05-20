#!/bin/bash
set -ex

#This script will run a number of tests in OpenCL mode, directly on the Sibernetic executable. 

no_gui="-no_g"

# Falling cube 

cube_dir=simulations/test_opencl_cube

if [ ! -d "$cube_dir" ]; then
    mkdir -p "$cube_dir"
fi

if [[ "$@" == *"-show_cube"* ]]; then
    ./Release/Sibernetic -l_from lpath=${cube_dir}
else
    ./Release/Sibernetic -f demo1 -l_to lpath=${cube_dir} timelimit=0.01 logstep=10 device=ALL ${no_gui} -q 
fi
