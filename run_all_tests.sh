#This script will run a number of tests using sibernetic_c302.py & check the 
# files produced
set -ex

# Ensure Python can locate packages installed via pip
PY_SITE=$(python3 - <<'EOF'
import site, sys
sys.stdout.write(site.getsitepackages()[0])
EOF
)
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}$PY_SITE:$(pwd)"
# Use the same Python that installed our packages
PY_BIN=$(command -v python3)
export PATH="$(dirname "$PY_BIN"):$PATH"

# Warn if pyneuroml is missing; the simulator can fall back to subprocess
python3 -c "import pyneuroml" 2>/dev/null || echo 'Warning: pyneuroml not found, using subprocess fallback'

# Cube only
python3 sibernetic_c302.py -test -noc302 -duration 30 -dt 0.005 -configuration demo1 -simName test_cube -logstep 10 -q

# No c302 but with the worm body 
python3 sibernetic_c302.py -test -noc302 -duration 0.1 -simName test_noc302 

# No c302 - test logstep
python3 sibernetic_c302.py -test -noc302 -duration 0.054 -logstep 3 -simName test_noc302_logstep -q

# c302
python3 sibernetic_c302.py -test  -duration 1.1  -c302params C1 -simName test_c302 -q

# c302 + half_resolution
python3 sibernetic_c302.py -test  -duration 1 -c302params C0 -configuration worm_alone_half_resolution -simName test_c302_half_resolution -q

# c302 + TestMuscle 
python3 sibernetic_c302.py -test -duration 20 -c302params C0 -reference TargetMuscle -configuration worm_alone_half_resolution -logstep 500 -simName test_c302_half_resolution_target_muscle -q

# c302 + C2 FW 
python sibernetic_c302.py -test -duration 15 -c302params C2 -reference FW -configuration worm_crawl_half_resolution -logstep 100

if [[ ($# -eq 1) && ($1 == '-all') ]]; then

    # Run a simulation with the FW (forward locomotion) c302 configuration with C2 (cond based) cells
    python sibernetic_c302.py -test -duration 150 -dt 0.005 -dtNrn 0.05 -logstep 100 -device=CPU -configuration worm_crawl_half_resolution -reference FW -c302params C2 -datareader UpdatedSpreadsheetDataReader2 -simName test_C2_FW -q

    # Run a simulation with the FW (forward locomotion) c302 configuration with W2D (simple passive) cells
    ##python sibernetic_c302.py -test -duration 15.0 -dt 0.005 -dtNrn 0.05 -logstep 100 -device=CPU -configuration worm_crawl_half_resolution -reference FW -c302params W2D -datareader UpdatedSpreadsheetDataReader2

    omv all -V

fi


# Run unit tests.
set +e
python3 -m pytest -q tests/test_energy.py tests/test_solver_logs.py
rc=$?
if [ "$rc" -ne 0 ] && [ "$rc" -ne 5 ]; then
    exit $rc
fi


