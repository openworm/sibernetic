#!/bin/bash
set -ex

#This script will run a number of tests in OpenCL mode, directly on the Sibernetic executable. 

no_gui="-no_g"
export_vtk="-export_vtk"

# Falling cube 

cube_dir=simulations/test_opencl_cube

if [ ! -d "$cube_dir" ]; then
    mkdir -p "$cube_dir"
fi

if [[ "$@" == *"-show_cube"* ]]; then
    ./Release/Sibernetic -l_from lpath=${cube_dir}
else
    if [[ "$@" != *"-show"* ]]; then # i.e. not showing any other example
        rm -f ${cube_dir}/* # remove any previous simulation files
        ./Release/Sibernetic -f demo1 -l_to lpath=${cube_dir} timelimit=0.02 logstep=10 device=ALL ${no_gui} -q ${export_vtk}
        # Test reloading the simulation files with Python viewer
        python SiberneticReplay.py ${cube_dir}/position_buffer.txt -nogui
    fi
fi

# 2 elastic membranes example

memb_dir=simulations/test_opencl_membranes

if [ ! -d "$memb_dir" ]; then
    mkdir -p "$memb_dir"
fi

if [[ "$@" == *"-show_membranes"* ]]; then
    ./Release/Sibernetic -l_from lpath=${memb_dir}
else
    if [[ "$@" != *"-show"* ]]; then # i.e. not showing any other example
        rm -f ${memb_dir}/* # remove any previous simulation files
        ./Release/Sibernetic -f demo2 -l_to lpath=${memb_dir} timelimit=0.02 logstep=10 device=ALL ${no_gui} -q ${export_vtk}
        # Test reloading the simulation files with Python viewer
        python SiberneticReplay.py ${memb_dir}/position_buffer.txt -nogui
    fi
fi
