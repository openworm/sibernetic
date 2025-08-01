#This script will run a number of tests using sibernetic_c302.py & check the 
# files produced
set -ex

# Ensure system site-packages are visible to embedded Python
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}/usr/local/lib/python3.12/dist-packages:$(pwd)"
export PATH="/usr/bin:$PATH"

# No c302
python3 sibernetic_c302.py -test -noc302 -duration 0.1

# No c302 - test logstep
python3 sibernetic_c302.py -test -noc302 -duration 0.054 -logstep 3

# c302
if command -v nrnivmodl >/dev/null 2>&1; then
    python3 sibernetic_c302.py -test  -duration 1.1  -c302params C1

    # c302 + half_resolution
    python3 sibernetic_c302.py -test  -duration 1  -c302params C0 -configuration worm_alone_half_resolution

    # c302 + TestMuscle
    python3 sibernetic_c302.py -test -duration 20 -c302params C0 -reference TargetMuscle -configuration worm_alone_half_resolution -logstep 500
else
    echo "Skipping c302 tests due to missing NEURON" >&2
fi

if [[ ($# -eq 1) && ($1 == '-all') ]]; then

    # Run a simulation with the FW (forward locomotion) c302 configuration with C2 (cond based) cells
    python sibernetic_c302.py -test -duration 15.0 -dt 0.005 -dtNrn 0.05 -logstep 100 -device=CPU -configuration worm_crawl_half_resolution -reference FW -c302params C2 -datareader UpdatedSpreadsheetDataReader2

    # Run a simulation with the FW (forward locomotion) c302 configuration with W2D (simple passive) cells
    ##python sibernetic_c302.py -test -duration 15.0 -dt 0.005 -dtNrn 0.05 -logstep 100 -device=CPU -configuration worm_crawl_half_resolution -reference FW -c302params W2D -datareader UpdatedSpreadsheetDataReader2

fi


