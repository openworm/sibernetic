import os
import sys
import numpy as np

import matplotlib.pyplot as plt

if os.getcwd().endswith("tests"):
    os.chdir("..")

sys.path.append(".")

from SibSimulation import SibSimulation


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


def plot_and_save_positions(sim_ref, sim_dir, show_plot=False, types_to_plot={}):

    print(
        f"Plotting positions for simulation: {sim_ref} from directory: {sim_dir}"
        + " and plotting "
        if show_plot
        else ""
    )

    if os.path.exists(os.path.join(sim_dir, "report.json")):
        sim = SibSimulation(report_file=os.path.join(sim_dir, "report.json"))
    elif os.path.exists(os.path.join(sim_dir, "position_buffer.txt")):
        sim = SibSimulation(position_file=os.path.join(sim_dir, "position_buffer.txt"))
    else:
        print(
            f"Files not found: {sim_dir}/position_buffer.txt or {sim_dir}/report.json"
        )
        return

    print("Plotting points: %s" % types_to_plot)

    for tp in types_to_plot:
        plt.figure()
        plt.title(f"Position of particles of type {tp} over time")
        plt.gcf().canvas.manager.set_window_title(
            f"Position of particles of type {tp} in: {sim_ref}"
        )
        for i in types_to_plot[tp]:
            times, positions = sim.get_particle_position(tp, i)
            if i == 0:
                # Save to file
                file_name = f"tests/positions_{sim_ref}_type_{tp}_index_{i}.dat"
                with open(file_name, "w") as f:
                    for t, pos in zip(times, positions):
                        f.write(f"{t}\t{pos[0]}\t{pos[1]}\t{pos[2]}\n")
                    f.write("\n")
                print(f"Saved positions for particle {i} of type {tp} to {file_name}")

            print(
                f"Retrieved positions for particle {i} of type {tp}: {len(positions)} points"
            )

            plt.plot(times, positions[:, 0], label=f"x ({i})", linestyle="-")
            plt.plot(times, positions[:, 1], label=f"y ({i})", linestyle=":")
            plt.plot(times, positions[:, 2], label=f"z ({i})", linestyle="--")

        plt.xlabel("Time")
        plt.ylabel("Position")
        plt.legend()

    if show_plot:
        plt.show()


# main method
if __name__ == "__main__":
    restructure_output_for_omv("simulations/test_C2_FW")
