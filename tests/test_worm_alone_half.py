import sys
import os
import numpy as np

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


def restructure_output_for_omv(sim_dir):
    """
    Restructure the output files for testing with OpenSourceBrain Model Validation (OMV).
    """

    report_file = os.path.join(sim_dir, "report.json")
    import json

    report_data = json.load(open(report_file, "r"))
    duration_s = (
        float(report_data["duration"].split()[0]) / 1000
    )  # Convert ms to seconds
    dt_s = float(report_data["dt"].split()[0]) / 1000  # Convert ms to seconds
    print(report_data)

    # Add times to the muscle activity buffer file
    src_file = os.path.join(sim_dir, "muscles_activity_buffer.txt")
    musc_dat = np.loadtxt(src_file)  # Ensure the file exists and is readable

    times = np.reshape(np.arange(0, duration_s, dt_s), (-1, 1))
    td = np.concatenate((times, musc_dat), axis=1)
    np.savetxt(
        os.path.join(sim_dir, "muscles_activity_buffer.dat"),
        td,
        fmt="%.6f",
        delimiter="\t",
    )


print("TEST: Saved simulation to: %s" % sim_dir)

restructure_output_for_omv(sim_dir)
