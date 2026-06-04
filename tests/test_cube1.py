import sys
import os

import matplotlib.pyplot as plt


if os.getcwd().endswith("tests"):
    os.chdir("..")

sys.path.append(".")

from SibSimulation import SibSimulation


def test_cube1():
    sim_ref = "test_cube"
    sim_report = f"simulations/{sim_ref}/report.json"
    sim = SibSimulation(report_file=sim_report)

    sim_ref = "test_opencl_demo1"
    sim = SibSimulation(position_file=f"simulations/{sim_ref}/position_buffer.txt")

    all_3D_points = sim.all_3D_points

    print(f"Loaded {len(all_3D_points)} time points of types {all_3D_points[0].keys()}")

    types_to_plot = [1.1, 2.1]

    for tp in types_to_plot:
        plt.figure()
        plt.title(f"Position of particles of type {tp} over time")
        for i in [0, 60, 124]:
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

    if "-gui" in sys.argv:
        plt.show()


# main
if __name__ == "__main__":
    test_cube1()
