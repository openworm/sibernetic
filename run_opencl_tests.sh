#!/bin/bash
set -ex

# This script will run a number of tests in OpenCL mode, directly on the Sibernetic executable.

no_gui="-no_g"
export_vtk="-export_vtk"

export PYTHONPATH=$PYTHONPATH:.

# Arguments passed to this script, used below to decide which example(s) to show
script_args="$*"

# Run (or show a previous run of) a single Sibernetic configuration.
#   $1 - Sibernetic config to run (passed to -f)
#   $2 - timelimit for the simulation
#   $3 - timestep for the simulation
#   $4 - logstep for the simulation
run_test () {
    local config="$1"
    local timelimit="$2"
    local timestep="$3"
    local logstep="$4"

    sim_dir="simulations/test_opencl_${config}"
    show_flag="-show_${config}"

    if [ ! -d "$sim_dir" ]; then
        mkdir -p "$sim_dir"
    fi

    if [[ "$script_args" == *"$show_flag"* ]]; then
        ./Release/Sibernetic -l_from lpath=${sim_dir}
    else
        if [[ "$script_args" != *"-show"* ]]; then # i.e. not showing any other example
            rm -f ${sim_dir}/* # remove any previous simulation files
            ./Release/Sibernetic -f ${config} -l_to lpath=${sim_dir} timelimit=${timelimit} timestep=${timestep} logstep=${logstep} device=ALL ${no_gui} -q ${export_vtk}
            # Test reloading the simulation files with Python viewer
            python3 SiberneticReplay.py ${sim_dir}/position_buffer.txt -nogui
        fi
    fi
}

# Demo 1: Falling cube
run_test demo1 0.03 2e-5 20

# Demo 2: Elastic membranes
run_test demo2 0.02 2e-5 20

# Worm alone, half resolution (for faster testing)
run_test worm_alone_half_resolution 0.02 2e-5 20

# Worm swimming, half resolution (for faster testing)
run_test worm_swim_half_resolution 0.05 2e-5 50

# Worm crawling, half resolution (for faster testing)
run_test worm_crawl_half_resolution 0.05 2e-5 50

# Full worm
run_test worm 0.01 2e-5 20

# Crawling worm
run_test worm_crawling 0.01 2e-5 20

# Worm no water
run_test worm_no_water 0.1 2e-5 20
