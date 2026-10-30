import sys
import os

from test_utils import restructure_output_for_omv

if os.getcwd().endswith("tests"):
    os.chdir("..")

sys.path.append(".")

from sibernetic_c302 import run

duration = 10
dt = 0.005  # Simulation time step

sim_dir, reportj = run(
    noc302=True,
    duration=duration,
    dt=dt,
    logstep=1,
    configuration="worm_alone_half_resolution",
    simName="test_worm_alone_half",  # Explicitly set simulation name
)

print("TEST: Saved simulation to: %s" % sim_dir)

restructure_output_for_omv(sim_dir)
