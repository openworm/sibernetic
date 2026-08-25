#!/bin/bash
set -ex

# This script will run a number of tests in OpenCL mode, directly on the Sibernetic executable.

# Arguments passed to this script, used below to decide which example(s) to show
script_args=("$@")

no_gui="-no_g"

# Check if a given exact argument was passed to this script
has_arg () {
    local target="$1"
    local a
    for a in "${script_args[@]}"; do
        [[ "$a" == "$target" ]] && return 0
    done
    return 1
}

# Check if any argument passed to this script starts with the given prefix
has_arg_prefix () {
    local prefix="$1"
    local a
    for a in "${script_args[@]}"; do
        [[ "$a" == "$prefix"* ]] && return 0
    done
    return 1
}

if has_arg "-vtk"; then
    export_vtk="-export_vtk"
else
    export_vtk=""
fi

if has_arg "-q" || has_arg "-quick"; then
    quick_run=true
else
    quick_run=false
fi

export PYTHONPATH=$PYTHONPATH:.


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

    local config_basename="${config##*/}"
    sim_dir="simulations/test_opencl_${config_basename}"
    show_flag="-show_${config_basename}"

    if [ ! -d "$sim_dir" ]; then
        mkdir -p "$sim_dir"
    fi

    if has_arg "$show_flag"; then
        ./Release/Sibernetic -l_from lpath=${sim_dir}
    else
        if ! has_arg_prefix "-show"; then # i.e. not showing any other example
            rm -f ${sim_dir}/* # remove any previous simulation files
            ./Release/Sibernetic -f ${config} -l_to lpath=${sim_dir} timelimit=${timelimit} timestep=${timestep} logstep=${logstep} device=ALL ${no_gui} -q ${export_vtk}
            
            # check if there is an OMV test
            if [ -f "tests/.test.${config}.omt" ]; then

                if ! omv test -V tests/.test.${config}.omt; then
                    echo "OMV test failed: tests/.test.${config}.omt" >> omv_report
                else
                    echo "OMV test passed: tests/.test.${config}.omt" >> omv_report
                fi
            else
                echo "No OMV test found at tests/.test.${config}.omt, skipping that part of the test" >> omv_report
            fi
            
            # Test reloading the simulation files with Python viewer
            python3 SiberneticReplay.py ${sim_dir}/position_buffer.txt -nogui
        fi
    fi
}

# Demo 1: Falling cube
run_test demo1 0.03 2e-5 20

# Demo 2: Elastic membranes
run_test demo2 0.02 2e-5 20


if [ "$quick_run" = false ]; then

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

    # Test configuration: liquid cube falling into liquid
    run_test configuration/test/v_test_liquid 0.05 2.5e-5 20

    # Test configuration: one elastic particle falling under gravity, i.e. on a spring
    run_test configuration/test/one_spring_test 0.005 1e-5 1

    # Test configuration: falling liquid
    run_test configuration/test/test_energy 0.03 2e-5 20

fi
