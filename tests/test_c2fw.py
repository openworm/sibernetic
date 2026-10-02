import sys
import os

from test_utils import restructure_output_for_omv

if os.getcwd().endswith("tests"):
    os.chdir("..")

sys.path.append(".")

from sibernetic_c302 import run

duration = 50
dt = 0.005  # Simulation time step

sim_dir, reportj = run(
    noc302=False,
    duration=duration,
    dt=dt,
    logstep=10,
    configuration="worm_alone_half_resolution",
    reference="FW",
    c302params="C2",
    simName="test_C2_FW_short",  # Explicitly set simulation name
)

print("TEST: Saved simulation to: %s" % sim_dir)

restructure_output_for_omv(sim_dir)
