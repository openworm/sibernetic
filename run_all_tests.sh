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
export PATH="/usr/bin:$PATH"

# Warn if pyneuroml is missing; the simulator can fall back to subprocess
python3 -c "import pyneuroml" 2>/dev/null || echo 'Warning: pyneuroml not found, using subprocess fallback'

# No c302
python3 sibernetic_c302.py -test -noc302 -duration 0.1

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

# Run unit tests.  Skip engine comparison unless requested
RUN_ENGINE_TESTS=0 python3 -m pytest -q tests/test_pytorch_solver.py tests/test_energy.py
RUN_ENGINE_TESTS=1 python3 -m pytest -q tests/test_torch_backend.py


