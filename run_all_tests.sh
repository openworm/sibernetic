#!/bin/bash
set -ex

#This script will run a number of tests using sibernetic_c302.py & check the files produced

# Cube only
python3 sibernetic_c302.py -test -noc302 -duration 30 -dt 0.005 -configuration demo1 -simName test_cube -logstep 10 -q

# No c302
python3 sibernetic_c302.py -test -noc302 -duration 0.1 -simName test_noc302 

# No c302 - test logstep
python3 sibernetic_c302.py -test -noc302 -duration 0.054 -logstep 3 -simName test_noc302_logstep -q

# c302
python3 sibernetic_c302.py -test  -duration 1.1  -c302params C1 -simName test_c302 -q

# c302 + half_resolution
python3 sibernetic_c302.py -test  -duration 1 -c302params C0 -configuration worm_alone_half_resolution -simName test_c302_half_resolution -q

# c302 + TestMuscle 
python3 sibernetic_c302.py -test -duration 20 -c302params C0 -reference TargetMuscle -configuration worm_alone_half_resolution -logstep 500 -simName test_c302_half_resolution_target_muscle -q

if [[ ($# -eq 1) && ($1 == '-all') ]]; then

    # Run a simulation with the FW (forward locomotion) c302 configuration with C2 (cond based) cells
    python sibernetic_c302.py -test -duration 150 -dt 0.005 -dtNrn 0.05 -logstep 100 -device=CPU -configuration worm_crawl_half_resolution -reference FW -c302params C2 -datareader UpdatedSpreadsheetDataReader2 -simName test_C2_FW-q

    # Run a simulation with the FW (forward locomotion) c302 configuration with W2D (simple passive) cells
    ##python sibernetic_c302.py -test -duration 15.0 -dt 0.005 -dtNrn 0.05 -logstep 100 -device=CPU -configuration worm_crawl_half_resolution -reference FW -c302params W2D -datareader UpdatedSpreadsheetDataReader2

    omv all -V

fi


