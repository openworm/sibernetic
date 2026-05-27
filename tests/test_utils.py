import os
import numpy as np


def restructure_output_for_omv(sim_dir):
    """
    Restructure the output files for testing with OpenSourceBrain Model Validation (OMV).
    """

    report_file = os.path.join(sim_dir, "report.json")
    print("Loading report file: %s" % report_file)
    import json

    report_data = json.load(open(report_file, "r"))
    duration_s = (
        float(report_data["duration"].split()[0]) / 1000
    )  # Convert ms to seconds
    dt_s = float(report_data["dt"].split()[0]) / 1000  # Convert ms to seconds
    logstep = int(report_data["logstep"])
    print(report_data)

    # Add times to the muscle activity buffer file
    src_file = os.path.join(sim_dir, "muscles_activity_buffer.txt")
    musc_dat = np.loadtxt(src_file)  # Ensure the file exists and is readable

    print("Loaded muscle activity buffer file: %s" % src_file)
    print("Muscle activity data shape: %s" % str(musc_dat.shape))

    times = np.reshape(np.arange(0, duration_s, dt_s * logstep), (-1, 1))
    print(f"Times: ({times[0]}, {times[-1]})")
    td = np.concatenate((times, musc_dat), axis=1)
    print(td)
    dat_file = os.path.join(sim_dir, "muscles_activity_buffer.dat")
    np.savetxt(
        dat_file,
        td,
        fmt="%.6f",
        delimiter="\t",
    )
    print("Saved restructured muscle activity data to: %s" % dat_file)


# main method
if __name__ == "__main__":
    restructure_output_for_omv("simulations/test_C2_FW")
