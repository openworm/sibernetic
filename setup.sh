#!/bin/bash
set -ex

# Install system dependencies for building Sibernetic and the PyTorch solver
apt-get update
apt-get install -y python3-dev ocl-icd-opencl-dev libglu1-mesa-dev freeglut3-dev libglew-dev clinfo pocl-opencl-icd

# Install Python packages
pip install torch ruff pytest pyneuroml || echo "Warning: failed to install pyneuroml"

# Verify pyneuroml installed correctly if available
python3 - <<'EOF'
try:
    import pyneuroml
    print(pyneuroml.__version__)
except Exception:
    print('pyneuroml not installed, continuing')
EOF
